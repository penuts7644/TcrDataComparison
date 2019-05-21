"""Merges CSV files together.

"""
import argparse
import glob
import os

import pandas


def _parse_commandline():
    """Parses the commandline arguments given by user.

    Returns
    -------
    ArgumentParser
        Object containing parsed commandline arguments.

    """

    # Parse the commandline arguments given by user.
    parser = argparse.ArgumentParser(description="Merges output CSV file in "
                                                 "a directory to one file.")
    parser.add_argument("directory", metavar="D", type=str,
                        help="Directory containing input CSV files")
    return parser.parse_args()


def main():
    """Main function when script ran through command-line."""

    # Parse arguments from user, locate all CSV files in the given directory
    # and merge them together.
    args = _parse_commandline()
    combined_df = pandas.DataFrame()
    for csv_file in glob.glob(os.path.abspath(args.directory) + "/*.csv"):
        df = pandas.read_csv(csv_file, sep=";", index_col=0, header=0)
        combined_df = pandas.concat([combined_df, df], axis=1, sort=False)

    # Writes the new dataframe to a CSV file.
    pandas.DataFrame.to_csv(combined_df, path_or_buf=os.path.join(
        os.getcwd(), "generated_merged.csv"), index=True, sep=";")
    print("Written '{}' file to '{}' directory".format("generated_merged.csv",
                                                       os.getcwd()))


if __name__ == "__main__":
    main()
