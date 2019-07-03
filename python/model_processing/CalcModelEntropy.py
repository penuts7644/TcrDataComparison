"""Calculates the entropy for between given IGoR model parameter and marginal
files. Kullback Leibler divergence is used to compare each model.

"""

import argparse
import itertools
import multiprocessing
import os

import numpy
import pandas

import _entropy
import _genmodel


def _parse_commandline():
    """Parses the commandline arguments given by user.

    Returns
    -------
    ArgumentParser
        Object containing parsed commandline arguments.

    """

    # Parse the commandline arguments given by user.
    parser = argparse.ArgumentParser(
        description="Calculates the entropy for between given IGoR model " \
                    "parameter and marginal files. Kullback Leibler " \
                    "divergence is used to compare each model.")
    parser.add_argument("-model", metavar=('<id>', '<parameters>', '<marginals>'),
                        type=str, required=True, action='append', nargs=3,
                        help="A unique model identifier (used in the output " \
                             "matrix data), IGoR parameters file and an IGoR " \
                             "marginals file.")
    parser.add_argument("--num-threads", type=int, nargs="?", default=1,
                        help="Number of threads in integer (default: %(default)s)")
    return parser.parse_args()


def compare_models(args):
    """Compares two IGoR models and returns the divergence value.

    Parameters
    ----------
    args : list
        A collection of arguments containing a list of model (id, parameters
        and marginals) permutations to be compared as well as the name of the
        event to process.

    Returns
    -------
    list
        A list of lists where each list contains the id for model 1, model 2
        and the computed divergence value of the two IGoR models.

    """
    # Process each model combination in the input data.
    df, kwargs = args
    type_list = kwargs['type_list']
    dkl_values = []
    for permutations in df:

        # Read in the generative models.
        model_1 = _genmodel.GenModel(permutations[0][1], permutations[0][2])
        model_2 = _genmodel.GenModel(permutations[1][1], permutations[1][2])

        # Calculate the entropy between the two models.
        entropy = _entropy.compute_models_dkl(
            model_1, model_2, type_list=type_list)
        dkl_values.append([permutations[0][0], permutations[1][0], entropy])
    return dkl_values


def multiprocess_dataframe(df, func, num_workers, **kwargs):
    """Applies multiprocessing on a pandas dataframe using given function.

    Parameters
    ----------
    df : list
        List like object to be split for multiple workers.
    func : Object
        A function object that the workers should apply.
    num_workers : int
        Integer specifying the number of workers (threads) to create.
    **kwargs : dict
        The remaining arguments to be given to the input function.

    Returns
    -------
    list
        Contains the results from each of the workers.

    """
    # Check out available worker count and adjust accordingly.
    if len(df) < num_workers:
        num_workers = len(df)

    # Divide the array into chucks for the workers.
    pool = multiprocessing.Pool(processes=num_workers)
    result = pool.map(func, [(d, kwargs)
                             for d in numpy.array_split(df, num_workers)])
    pool.close()
    return result


def main():
    """Main function when script ran through command-line."""

    # Parse arguments and make model permutations.
    args = _parse_commandline()
    model_permutations = list(itertools.permutations(args.model, 2))
    model_identifiers = [model[0] for model in args.model]
    dkl_comparisons = {
        'D_trim_5': [r'Deletion_D_gene_Five_prime_.*'],
        'dinuc_markov_VD': [r'DinucMarkov_VD_genes_Undefined_side_.*'],
        'J_trim_5': [r'Deletion_J_gene_Five_prime_.*'],
        'DJ_insert_length': [r'Insertion_DJ_gene_Undefined_side_.*'],
        'D_trim_3': [r'Deletion_D_gene_Three_prime_.*'],
        'J': [r'GeneChoice_J_gene_Undefined_side_.*'],
        'VD_insert_length': [r'Insertion_VD_genes_Undefined_side_.*'],
        'V': [r'GeneChoice_V_gene_Undefined_side_.*'],
        'dinuc_markov_DJ': [r'DinucMarkov_DJ_gene_Undefined_side_.*'],
        'V_trim_3': [r'Deletion_V_gene_Three_prime_.*'],
        'D': [r'GeneChoice_D_gene_Undefined_side_.*'],
        'total': None
    }

    # For each event type (and combined model) calculate DKL.
    for name, event_types in dkl_comparisons.iteritems():

        # Process the dataframe in various threads and append results.
        results = multiprocess_dataframe(
            df=model_permutations, func=compare_models,
            num_workers=args.num_threads, type_list=event_types)
        results = [i for j in results for i in j]

        # Create the zero-d divergence matrix and fill it up with data.
        # NOTE - READ ROW x COLUMN
        matrix_df = pandas.DataFrame(
            0, columns=model_identifiers, index=model_identifiers, dtype=float)
        for value in results:
            matrix_df[value[1]][value[0]] = value[2]

        # Writes the new dataframe to a CSV file.
        output_filename = os.path.join(
            os.getcwd(), 'calc_model_entropy_{}.tsv'.format(name))
        pandas.DataFrame.to_csv(matrix_df, path_or_buf=output_filename,
                                sep="\t", na_rep="NA", index=True)
        print("Written '{}' file".format(output_filename))


if __name__ == "__main__":
    main()
