"""Checks and locates any stop codons in the generated nucleotide sequences.

"""

import argparse
import itertools
import multiprocessing
import os
import re
import sys

import pandas
import numpy


def _parse_commandline():
    """Parses the commandline arguments given by user.

    Returns
    -------
    ArgumentParser
        Object containing parsed commandline arguments.

    """

    # Parse the commandline arguments given by user.
    parser = argparse.ArgumentParser(description="Locates any stop codons in "
                                                 "the given sequences.")
    parser.add_argument("file", metavar="F", type=str,
                        help="File containing sequences separated by a ';'")
    parser.add_argument("column", metavar="C", type=str,
                        help="The column name in the file containing the "
                             "nucleotide of amino acid sequences")
    parser.add_argument("--num-threads", type=int, nargs="?", default=1,
                        help="Number of threads in integer (default: %(default)s)")
    return parser.parse_args()


def locate_stop_codons_sequence(seq):
    """Locates the stop codons for a sequence string.

    Parameters
    ----------
    seq : str
        Sequence in string format.

    Returns
    -------
    list
        Contains tuples with matches as following (start position, end
        position, codon string).

    """

    # Locate the codons with pattern in the sequence string
    stop_codon_pattern = re.compile("[UT]AA|[UT]AG|[UT]GA")
    matches = [(m.start(0), m.end(0), m.group(0))
               for m in re.finditer(stop_codon_pattern, seq)]
    return matches


def locate_stop_codons(args):
    """Locates the stop codons for each sequence in a dataframe.

    Parameters
    ----------
    args : list
        A collection of arguments containing dataframe, seq_col.

    Returns
    -------
    list
        Containing lists with matches and index value for the row.

    """
    df, kwargs = args
    data_list = list()
    for i, row in df.iterrows():
        matches = locate_stop_codons_sequence(seq=row[kwargs["seq_col"]])
        data_list.append([i, matches])
    return data_list


def multiprocess_dataframe(df, func, num_workers, **kwargs):
    """Applies multiprocessing on a pandas dataframe using given function.

    Parameters
    ----------
    df : DataFrame
        Pandas dateframe to be split for multiple workers.
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
    return list(result)


def main():
    """Main function when script ran through command-line."""

    # TODO: Some CSV files are not consistent with the separator ';' or ','.

    # Parse arguments and read in file.
    args = _parse_commandline()
    df = pandas.read_csv(args.file, sep=";", index_col="seq_index", header=0)
    output_filename = args.file.split("/")[-1].split(".")[0] \
        + "_stop_codon_indices.csv"

    # Produce error if column not exists.
    if args.column not in df.columns:
        print("Error: given column name does not exit in dataframe")
        sys.exit(0)

    # Locate stop codons and create new dataframe with results.
    matches = itertools.chain(
        *multiprocess_dataframe(df=df, func=locate_stop_codons,
                                num_workers=args.num_threads,
                                seq_col=args.column))
    new_df = pandas.DataFrame.from_records(matches,
                                           columns=[df.index.name,
                                                    "sequence_stop_codons"])

    # Write resulting dataframe to CSV file.
    pandas.DataFrame.to_csv(new_df, path_or_buf=os.path.join(
        os.getcwd(), output_filename), index=False, sep=";")
    print("Written '{}' file to '{}' directory".format(output_filename,
                                                       os.getcwd()))


if __name__ == "__main__":
    main()
