#!/usr/bin/env bash

# Usage docstring.
USAGE="
$(basename "$0") [-h] [-t num_threads] [-s name] [-c name] [-I csv_file] [-A] [-E num_scenarios] [-G num_sequences]

The following options are allowed:
    -h    help        show this help text
    -t    threads     number of threads to use in IGoR                                     (default: 1)
    -s    species     species name                                                         (default: human)
    -c    chain       chain name                                                           (default: alpha)
    -I    index       index the given CSV file with sequences
    -A    align       aligns all segments (V, D and J) against indexed sequences           (requires: -s, -c)
    -E    evaluate    evaluate the created alignment and produce scenarios                 (requires: -s, -c)
    -G    generate    generate number of error free sequences (+ CDR3) using IGoR model    (requires: -s, -c)

NOTE: Options -A, -E and -G require options -s and -c to be specified, otherwise default options are
      used in script. LOWERCASE options need to be specified before using UPPERCASE options.
"

# Reset getopts.
OPTIND=1

# Set some general variables.
WD_PATH=$PWD
THREADS=1
SEQS_PATH=""
SPECIES="human"
CHAIN="alpha"
GENERATE=0
SCENARIOS=0
IGOR="igor -set_wd ${WD_PATH}"


# Show usage if no arguments given, else parse arguments.
if [[ $# == 0 ]]; then
    echo "$USAGE"
    exit
else
    while getopts ":ht:s:c:I:AE:G:" opt; do
        case "$opt" in
            h)
                echo "$USAGE" >&2
                exit
                ;;
            t)
                THREADS=$OPTARG >&2
                ;;
            s)
                SPECIES=$OPTARG >&2
                ;;
            c)
                CHAIN=$OPTARG >&2
                ;;
            I)
                SEQS_PATH=$OPTARG >&2
                eval "${IGOR} -read_seqs ${SEQS_PATH}" # Read the given sequences.
                ;;
            A)
                SPECIES_CHAIN="-species ${SPECIES} -chain ${CHAIN}"
                eval "${IGOR} ${SPECIES_CHAIN} -align --all" # Align indexed sequences.
                ;;
            E)
                SPECIES_CHAIN="-species ${SPECIES} -chain ${CHAIN}"
                SCENARIOS=$OPTARG >&2
                eval "${IGOR} ${SPECIES_CHAIN} -evaluate -output --Pgen --scenarios ${SCENARIOS}" # Evaluate and generate output data.
                ;;
            G)
                SPECIES_CHAIN="-species ${SPECIES} -chain ${CHAIN}"
                GENERATE=$OPTARG >&2
                eval "${IGOR} ${SPECIES_CHAIN} -generate ${GENERATE} --CDR3 --noerr" # Generate sequences.
                ;;
            \?)
                echo "invalid option: -$OPTARG" >&2
                echo "$USAGE" >&2
                exit 1
                ;;
            :)
                echo "option -$OPTARG requires an argument." >&2
                echo "$USAGE" >&2
                exit 1
                ;;
            *)
                echo "unimplemented option: -$OPTARG" >&2;
                exit 1
                ;;
        esac
    done
fi

shift $((OPTIND - 1))