"""Checks and locates any stop codons in the generated nt/CDR3 sequences.

Usage
-----
    First argument points to the file containing the sequences (either
    nt_sequence or nt_CDR3 column).
"""

import os
import re
import sys

import pandas


def main():
    """Main function when script ran through command-line."""
    # Locate all stop codons in the given file and return the file index,
    # codon start and end index as well as the codon itself for each of them.
    print("Given sequences file path: {}".format(sys.argv[1]))
    df = pandas.read_csv(sys.argv[1], sep=";", index_col="seq_index", header=0)
    output_filename = sys.argv[1].split("/")[-1].split(".")[0] \
        + "_stop_codon_indices.csv"
    data_list = list()
    for i, row in df.iterrows():
        stop_codon_pattern = re.compile("[UT]AA|[UT]AG|[UT]GA")
        matches = [(m.start(0), m.end(0), m.group(0))
                   for m in re.finditer(stop_codon_pattern,
                                        row["nt_sequence"])]
        data_list.append([i, matches])
    new_df = pandas.DataFrame.from_records(data_list,
                                           columns=[df.index.name,
                                                    "sequence_stop_codons"])
    pandas.DataFrame.to_csv(new_df, path_or_buf=os.path.join(
        os.getcwd(), output_filename), index=False, sep=";")
    print("Written '{}' file to '{}' directory".format(output_filename,
                                                       os.getcwd()))


if __name__ == "__main__":
    main()
