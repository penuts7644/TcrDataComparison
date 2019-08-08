#!/bin/bash

#SBATCH --account=nn9603k

#SBATCH --job-name=ImmunoProbs

#SBATCH --time=0-01:00:00

#SBATCH --partition=bigmem

#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=32

#SBATCH --mem-per-cpu=8G

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

# Copy over files to work dir
cd ${SCRATCH}
cp "${SLURM_SUBMIT_DIR}/TcrDataComparison/immuno_probs_config.ini" .
cp "${SLURM_SUBMIT_DIR}/process_5_files/${SUBSET_ID}/combined_data_extract_and_models/sequence_data_extractor_CDR3.tsv" .
cp -r "${SLURM_SUBMIT_DIR}/human_TRB" .

# Mark the output dir for automatic copying to $SLURM_SUBMIT_DIR afterwards
savefile process_5_files

# Calculate entropy for all the models
mkdir process_5_files
cd process_5_files
mkdir ${SUBSET_ID}
cd ${SUBSET_ID}
mkdir evaluated_CDR3
cd evaluated_CDR3
cp -r "${SLURM_SUBMIT_DIR}/process_5_files/${SUBSET_ID}/combined_data_extract_and_models/evaluate" .
mv evaluate evaluate_combined

# Evaluate using the build-in igor model
mkdir evaluate_igor
cd evaluate_igor
cp -r "${SLURM_SUBMIT_DIR}/igor_human_TCRB_model" .
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'pgen_estimate_CDR3' -config-file '../../../../immuno_probs_config.ini' evaluate -custom-model 'igor_human_TCRB_model/model_params.txt' 'igor_human_TCRB_model/model_marginals.txt' -seqs '../../../../sequence_data_extractor_CDR3.tsv' -type 'beta' -cdr3 -anchor V '../../../../human_TRB/V_gene_CDR3_anchors.tsv' -anchor J '../../../../human_TRB/J_gene_CDR3_anchors.tsv'
rm -r igor_human_TCRB_model
cd ../

# Exit succesfully
exit 0
