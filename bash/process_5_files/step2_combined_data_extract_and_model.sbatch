#!/bin/bash

#SBATCH --account=nn9603k

#SBATCH --job-name=ImmunoProbs

#SBATCH --time=1-00:00:00

#SBATCH --partition=bigmem

#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=32

#SBATCH --mem-per-cpu=14G

# Set up job environment, load modules and setup python
module restore system
module use .local/easybuild/modules/all
module load IGoR/1.3.0-GCC-7.3.0-2.30
module load Python/2.7.15-intel-2018b
# Make sure to have ImmunoProbs installed for python with:
# pip install --user immuno-probs

# Set safer defaults for bash
set -o errexit
set -o nounset

SUBSET_ID='all'
# SUBSET_ID='100000'
# SUBSET_ID='50000'
# SUBSET_ID='10000'
# SUBSET_ID='5000'
# SUBSET_ID='1000'
# SUBSET_ID='500'

# Copy over files to work dir
cd ${SCRATCH}
cp "${SLURM_SUBMIT_DIR}/TcrDataComparison/immuno_probs_config.ini" .
cp -r "${SLURM_SUBMIT_DIR}/human_TRB" .

# Mark the output dir for automatic copying to $SLURM_SUBMIT_DIR afterwards
savefile process_5_files

# Create the job dir for the output files
mkdir process_5_files
cd process_5_files
mkdir ${SUBSET_ID}
cd ${SUBSET_ID}
mkdir combined_data_extract_and_models
cd combined_data_extract_and_models

# Combine the data extracts from the previous models
NUMBER_OF_FILES=`ls "${SLURM_SUBMIT_DIR}/process_5_files/${SUBSET_ID}/data_extract_and_models" | wc -l`
NUMBER_OF_FILES=`expr ${NUMBER_OF_FILES} - 1`
for i in $(seq 0 $NUMBER_OF_FILES)
do
    if [ $i == 0 ]; then
        cp "${SLURM_SUBMIT_DIR}/process_5_files/${SUBSET_ID}/data_extract_and_models/immuno_probs_${i}/converted_full_length_productive.tsv" .
        cp "${SLURM_SUBMIT_DIR}/process_5_files/${SUBSET_ID}/data_extract_and_models/immuno_probs_${i}/converted_full_length_unproductive.tsv" .
        cp "${SLURM_SUBMIT_DIR}/process_5_files/${SUBSET_ID}/data_extract_and_models/immuno_probs_${i}/converted_full_length.tsv" .
        cp "${SLURM_SUBMIT_DIR}/process_5_files/${SUBSET_ID}/data_extract_and_models/immuno_probs_${i}/converted_CDR3.tsv" .
    else
        tail -n +2 -q "${SLURM_SUBMIT_DIR}/process_5_files/${SUBSET_ID}/data_extract_and_models/immuno_probs_${i}/converted_full_length_productive.tsv" >> converted_full_length_productive.tsv
        tail -n +2 -q "${SLURM_SUBMIT_DIR}/process_5_files/${SUBSET_ID}/data_extract_and_models/immuno_probs_${i}/converted_full_length_unproductive.tsv" >> converted_full_length_unproductive.tsv
        tail -n +2 -q "${SLURM_SUBMIT_DIR}/process_5_files/${SUBSET_ID}/data_extract_and_models/immuno_probs_${i}/converted_full_length.tsv" >> converted_full_length.tsv
        tail -n +2 -q "${SLURM_SUBMIT_DIR}/process_5_files/${SUBSET_ID}/data_extract_and_models/immuno_probs_${i}/converted_CDR3.tsv" >> converted_CDR3.tsv
    fi
done

# Create a model (params and marginals) for each data extract
mkdir model
cd model
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'unproductive' -config-file '../../../../immuno_probs_config.ini' build -ref V '../../../../human_TRB/TRBV.fasta' -ref D '../../../../human_TRB/TRBD.fasta' -ref J '../../../../human_TRB/TRBJ.fasta' -seqs '../converted_full_length_unproductive.tsv' -n-iter 10 -type 'beta'
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'productive' -config-file '../../../../immuno_probs_config.ini' build -ref V '../../../../human_TRB/TRBV.fasta' -ref D '../../../../human_TRB/TRBD.fasta' -ref J '../../../../human_TRB/TRBJ.fasta' -seqs '../converted_full_length_productive.tsv' -n-iter 10 -type 'beta'
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'all' -config-file '../../../../immuno_probs_config.ini' build -ref V '../../../../human_TRB/TRBV.fasta' -ref D '../../../../human_TRB/TRBD.fasta' -ref J '../../../../human_TRB/TRBJ.fasta' -seqs '../converted_full_length.tsv' -n-iter 10 -type 'beta'
cd ../

# Evaluate the CDR3 sequences with each model
mkdir evaluate
cd evaluate
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'pgen_estimate_unproductive_CDR3' -config-file '../../../../immuno_probs_config.ini' evaluate -custom-model '../model/unproductive_params.txt' '../model/unproductive_marginals.txt' -seqs '../converted_CDR3.tsv' -type 'beta' -cdr3 -anchor V '../../../../human_TRB/V_gene_CDR3_anchors.tsv' -anchor J '../../../../human_TRB/J_gene_CDR3_anchors.tsv'
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'pgen_estimate_productive_CDR3' -config-file '../../../../immuno_probs_config.ini' evaluate -custom-model '../model/productive_params.txt' '../model/productive_marginals.txt' -seqs '../converted_CDR3.tsv' -type 'beta' -cdr3 -anchor V '../../../../human_TRB/V_gene_CDR3_anchors.tsv' -anchor J '../../../../human_TRB/J_gene_CDR3_anchors.tsv'
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'pgen_estimate_all_CDR3' -config-file '../../../../immuno_probs_config.ini' evaluate -custom-model '../model/all_params.txt' '../model/all_marginals.txt' -seqs '../converted_CDR3.tsv' -type 'beta' -cdr3 -anchor V '../../../../human_TRB/V_gene_CDR3_anchors.tsv' -anchor J '../../../../human_TRB/J_gene_CDR3_anchors.tsv'
cd ../

# Exit succesfully
exit 0
