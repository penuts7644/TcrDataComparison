#!/bin/bash

# Specify the initial subset IDs
subsetIDArray=(
    'all'
    '100000'
    '50000'
    '10000'
    '5000'
    '1000'
    '500'
)

# Make initial output file
echo -e "id\tmodel\ttype\tcount" > "process_5_files/sequence_counts.tsv"

# For each subset ID process data files
NUMBER_OF_SUBSETS=`ls "process_5_files" | wc -l`
NUMBER_OF_SUBSETS=`expr ${NUMBER_OF_SUBSETS} - 2`
for i in $(seq 0 $NUMBER_OF_SUBSETS)
do

    # For each data extract count the number of sequences used for training
    NUMBER_OF_FILES=`ls "process_5_files/${subsetIDArray[${i}]}/data_extract_and_models" | wc -l`
    NUMBER_OF_FILES=`expr ${NUMBER_OF_FILES} - 1`
    for j in $(seq 0 $NUMBER_OF_FILES)
    do
        echo -e "${subsetIDArray[${i}]}\t${j}\tproductive\t`tail -n +2 "process_5_files/${subsetIDArray[${i}]}/data_extract_and_models/immuno_probs_${j}/converted_full_length_productive.tsv" | wc -l`" >> "process_5_files/sequence_counts.tsv"
        echo -e "${subsetIDArray[${i}]}\t${j}\tunproductive\t`tail -n +2 "process_5_files/${subsetIDArray[${i}]}/data_extract_and_models/immuno_probs_${j}/converted_full_length_unproductive.tsv" | wc -l`" >> "process_5_files/sequence_counts.tsv"
        echo -e "${subsetIDArray[${i}]}\t${j}\tall\t`tail -n +2 "process_5_files/${subsetIDArray[${i}]}/data_extract_and_models/immuno_probs_${j}/converted_full_length.tsv" | wc -l`" >> "process_5_files/sequence_counts.tsv"
        echo -e "${subsetIDArray[${i}]}\t${j}\tevaluated\t`tail -n +2 "process_5_files/${subsetIDArray[${i}]}/data_extract_and_models/immuno_probs_${j}/converted_CDR3.tsv" | wc -l`" >> "process_5_files/sequence_counts.tsv"
    done

    # Add combined data extract count
    echo -e "${subsetIDArray[${i}]}\tcombined\tproductive\t`tail -n +2 "process_5_files/${subsetIDArray[${i}]}/combined_data_extract_and_models/converted_full_length_productive.tsv" | wc -l`" >> "process_5_files/sequence_counts.tsv"
    echo -e "${subsetIDArray[${i}]}\tcombined\tunproductive\t`tail -n +2 "process_5_files/${subsetIDArray[${i}]}/combined_data_extract_and_models/converted_full_length_unproductive.tsv" | wc -l`" >> "process_5_files/sequence_counts.tsv"
    echo -e "${subsetIDArray[${i}]}\tcombined\tall\t`tail -n +2 "process_5_files/${subsetIDArray[${i}]}/combined_data_extract_and_models/converted_full_length.tsv" | wc -l`" >> "process_5_files/sequence_counts.tsv"
    echo -e "${subsetIDArray[${i}]}\tcombined\tevaluated\t`tail -n +2 "process_5_files/${subsetIDArray[${i}]}/combined_data_extract_and_models/converted_CDR3.tsv" | wc -l`" >> "process_5_files/sequence_counts.tsv"
done
