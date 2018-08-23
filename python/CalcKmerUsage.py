"""Calculates the kmer usage between simulated (IGoR) and experimental (e.g.
blood sample) data in CSV format.

"""

import argparse
import multiprocessing
import numpy
from collections import Counter
import glob
import os
import sys

import pandas


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
                        help="Experimental ';' separated file (required)")
    parser.add_argument("simulated", metavar="S", type=str,
                        help="Simulated ';' separated IGoR file (required)")
    parser.add_argument("--kmer-length", type=int, nargs="?", default=3,
                        help="Length of the Kmer's as integer (default: 3)")
    parser.add_argument("--num-threads", type=int, nargs="?", default=1,
                        help="Number of threads in integer (default: 1)")
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

    print(exp_counts, sim_counts)


if __name__ == "__main__":
    main()
