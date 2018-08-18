"""Translate IGoR output realization data to human readable format.

Usage
-----
    First argument points to the realizations file.
    Second argument points to the model parameters file that was used for the
    construction of the realizations file.
"""

import os
import sys

import pandas

from pygor.counters.bestscenarios import bestscenarios


def main():
    """Main function when script ran through command-line."""
    # Convert the realizations csv output file to sensible output data using
    # the corresponding model its parms text file.
    print("Given realizations file path: {}".format(sys.argv[1]))
    print("Given model parameters file path: {}".format(sys.argv[2]))
    df = bestscenarios.read_bestscenarios_values(scenarios_file=sys.argv[1],
                                                 model_parms_file=sys.argv[2])
    output_filename = sys.argv[1].split("/")[-1].split(".")[0] \
        + "_translated.csv"
    pandas.DataFrame.to_csv(df, path_or_buf=os.path.join(
        os.getcwd(), output_filename), index=False, sep=";")
    print("Written '{}' file to '{}' directory".format(output_filename,
                                                       os.getcwd()))


if __name__ == "__main__":
    main()
