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


def compare_models(df):
    """Compares two IGoR models and returns the divergence value.

    Parameters
    ----------
    df : list
        A collection of arguments containing a list of model (id, parameters
        and marginals) combinations to be compared.

    Returns
    -------
    list
        A list of lists where each list contains the id for model 1, model 2
        and the computed divergence value of the two IGoR models.

    """
    # Process each model combination in the input data.
    dkl_values = []
    for combination in df:

        # Read in the generative models.
        model_1 = _genmodel.GenModel(combination[0][1], combination[0][2])
        model_2 = _genmodel.GenModel(combination[1][1], combination[1][2])

        # Calculate the entropy between the two models.
        entropy = _entropy.compute_models_dkl(model_1, model_2)
        dkl_values.append([combination[0][0], combination[1][0], entropy])
    return dkl_values


def multiprocess_dataframe(df, func, num_workers):
    """Applies multiprocessing on a pandas dataframe using given function.

    Parameters
    ----------
    df : list
        List like object to be split for multiple workers.
    func : Object
        A function object that the workers should apply.
    num_workers : int
        Integer specifying the number of workers (threads) to create.

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
    result = pool.map(func, [d for d in numpy.array_split(df, num_workers)])
    pool.close()
    return result


def main():
    """Main function when script ran through command-line."""

    # Parse arguments and make model combinations.
    args = _parse_commandline()
    model_combinations = list(itertools.combinations(args.model, 2))
    model_identifiers = [model[0] for model in args.model]

    # Process the dataframe in various threads and append results.
    results = multiprocess_dataframe(
        df=model_combinations, func=compare_models, num_workers=args.num_threads)
    results = [i for j in results for i in j]

    # Create the zerod divergence matrix and fill it up with data.
    matrix_df = pandas.DataFrame(0, columns=model_identifiers, index=model_identifiers, dtype=float)
    for value in results:
        matrix_df[value[0]][value[1]] = value[2]
        matrix_df[value[1]][value[0]] = value[2]

    # Writes the new dataframe to a CSV file.
    output_filename_1 = os.path.join(os.getcwd(), 'calc_model_entropy.tsv')
    pandas.DataFrame.to_csv(matrix_df, path_or_buf=output_filename_1,
                            sep="\t", na_rep="NA", index=True)
    print("Written '{}' file".format(output_filename_1))


if __name__ == "__main__":
    main()
