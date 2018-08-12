"""Merges the output CSV files from one IGoR run into one file.

Usage
-----
    One argument that specifies the directory from which to locate and combine
    all the CSV files.
"""

import glob
import os
import sys

import pandas


def main():
    """Main function when script ran through command-line."""
    # Locate all CSV files in a given directory, merges them together and
    # writes the new dataframe to a csv file.
    print("Given directory path: {}".format(os.path.abspath(sys.argv[1])))
    combined_df = pandas.DataFrame()
    for csv_file in glob.glob(os.path.abspath(sys.argv[1]) + "/*.csv"):

        # TODO: Some CSV files are not consistent with the separator ';' or ','

        df = pandas.read_csv(csv_file, sep=";", index_col="seq_index",
                             header=0)
        combined_df = pandas.concat([combined_df, df], axis=1, sort=False)
    pandas.DataFrame.to_csv(combined_df, path_or_buf=os.path.join(
        os.getcwd(), "generated_merged.csv"), index=True, sep=";")
    print("Written '{}' file to '{}' directory".format("generated_merged.csv",
                                                       os.getcwd()))


if __name__ == "__main__":
    main()
