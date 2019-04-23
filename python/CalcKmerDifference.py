"""Calculates the kmer usage between simulated (IGoR) and experimental (e.g.
blood sample) data in CSV format.

"""

import argparse
import multiprocessing
from collections import Counter
import glob
import os
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
    parser = argparse.ArgumentParser(description="Calculates and compares the "
                                                 "Kmer usage between both "
                                                 "experimental and simulated "
                                                 "datasets.")
    parser.add_argument("experimental", metavar="E", type=str,
                        help="Experimental ';' separated file")
    parser.add_argument("simulated", metavar="S", type=str,
                        help="Simulated ';' separated IGoR file")
    parser.add_argument("--kmer-length", type=int, nargs="?", default=3,
                        help="Length of the Kmer's as integer (default: %(default)s)")
    parser.add_argument("--num-threads", type=int, nargs="?", default=1,
                        help="Number of threads in integer (default: %(default)s)")
    return parser.parse_args()


def count_kmers_sequence(seq, kmer_len):
    """Counts the occurrences of each Kmer within a given sequence.

    Parameters
    ----------
    seq : string
        Either nucleotide or amino acid sequence string.
    kmer_len : int
        The length of the Kmer's.

    Returns
    -------
    Counter
        Object containing the counts of each Kmer within the sequence.

    """
    return Counter(seq[start:start + kmer_len]
                   for start in xrange(len(seq) - kmer_len))


def count_kmers(args):
    """Counts the combined occurrences of each Kmer for pandas.DataFrame rows.

    Parameters
    ----------
    args : list
        A collection of arguments containing dataframe, column_name and Kmer
        length.

    Returns
    -------
    Counter
        Object containing the combined counts of each Kmer from multiple
        sequences.

    """
    df, kwargs = args
    kmer_counts = Counter()
    for _, row in df.iterrows():
        kmer_counts.update(count_kmers_sequence(seq=row[kwargs["col_name"]],
                                                kmer_len=kwargs["kmer_len"]))
    return kmer_counts


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
    pool = multiprocessing.Pool(processes=num_workers)
    result = pool.map(func, [(d, kwargs)
                             for d in numpy.array_split(df, num_workers)])
    pool.close()
    return list(result)


def calc_count_difference(c1, c2):
    """Calculates the difference between two counters.

    Parameters
    ----------
    c1 : Counter
        The first Counter object with counts.
    c2 : Counter
        The second Counter object with counts.

    Returns
    -------
    dict
        Dictionary containing the count differences between the two counters.

    """
    diff_dict = {}
    for kmer, count in c1.iteritems():
        if kmer in c2:
            diff_dict[kmer] = abs(count - c2[kmer])
        else:
            diff_dict[kmer] = count
    for kmer, count in c2.iteritems():
        if kmer not in c1:
            diff_dict[kmer] = count
    return diff_dict


def main():
    """Main function when script ran through command-line."""

    # Parse arguments and read in given files.
    args = _parse_commandline()
    exp_df = pandas.read_csv(args.experimental, sep=";", header=0)
    sim_df = pandas.read_csv(args.simulated, sep=";", index_col="seq_index",
                             header=0)

    # Create counters to store all Kmer counts and process both dataframes.
    exp_counts, sim_counts = Counter(), Counter()
    for i in multiprocess_dataframe(df=exp_df, func=count_kmers,
                                    num_workers=args.num_threads,
                                    col_name="N..Seq..CDR3.1",
                                    kmer_len=args.kmer_length):
        exp_counts.update(i)
    for i in multiprocess_dataframe(df=sim_df,
                                    func=count_kmers,
                                    num_workers=args.num_threads,
                                    col_name="nt_CDR3",
                                    kmer_len=args.kmer_length):
        sim_counts.update(i)

    # Calculate the count differences and create dataframe.
    df_data, count = dict(), 0
    for k, v in calc_count_difference(exp_counts, sim_counts).iteritems():
        df_data[count] = [k, exp_counts[k], sim_counts[k], v]
        count += 1
    count_df = pandas.DataFrame.from_dict(
        data=df_data, orient="index",
        columns=["kmer",
                 os.path.basename(os.path.splitext(args.experimental)[0]),
                 os.path.basename(os.path.splitext(args.simulated)[0]),
                 "diff_count"])
    count_df.index.name = "index"

    # Writes the new dataframe to a CSV file.
    pandas.DataFrame.to_csv(count_df, path_or_buf=os.path.join(
        os.getcwd(), "kmer_count.csv"), index=True, sep=";")
    print("Written '{}' file to '{}' directory".format("kmer_count.csv",
                                                       os.getcwd()))


if __name__ == "__main__":
    main()
