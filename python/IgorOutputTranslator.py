"""Translate IGoR output realization data to human readable format with pygor.

"""
import argparse
import csv
import os
import sys

import pandas

from pygor.counters.bestscenarios import bestscenarios


def _parse_commandline():
    """Parses the commandline arguments given by user.

    Returns
    -------
    ArgumentParser
        Object containing parsed commandline arguments.

    """

    # Parse the commandline arguments given by user.
    parser = argparse.ArgumentParser(description="Translates IGoR output "
                                                 "realization data to human "
                                                 "readable format with pygor.")
    parser.add_argument("realizations", metavar="R", type=str,
                        help="Realizations data ';' separated file")
    parser.add_argument("model_parms", metavar="M", type=str,
                        help="IGoR Model parameters text file used for generating "
                             "the realizations file")
    return parser.parse_args()


def main():
    """Main function when script ran through command-line."""

    # Convert the realizations csv output file to sensible output data using
    # the corresponding model's parameters text file.
    args = _parse_commandline()
    df = bestscenarios.read_bestscenarios_values(
        scenarios_file=args.realizations, model_parms_file=args.model_parms)
    output_filename = args.realizations.split("/")[-1].split(".")[0] \
        + "_translated.csv"

    # Write resulting dataframe to CSV file.
    pandas.DataFrame.to_csv(df, path_or_buf=os.path.join(
        os.getcwd(), output_filename), index=False, sep=";")
    print("Written '{}' file to '{}' directory".format(output_filename,
                                                       os.getcwd()))


if __name__ == "__main__":
    main()
