#!/bin/bash

#SBATCH --account=nn9603k
#SBATCH --job-name=ImmunoProbs
#SBATCH --array=0-4
#SBATCH --time=0-06:00:00
#SBATCH --nodes=5 --ntasks-per-node=1 --cpus-per-task=32


###########################
# Load necessary modules
###########################
module restore system
module use .local/easybuild/modules/all
module load IGoR/1.3.0-GCC-7.3.0-2.30
module load Python/2.7.15-intel-2018b
set -o errexit
set -o nounset


###########################
# Specify job variables
###########################
filesArray=(
    "${SLURM_SUBMIT_DIR}/BruskoT1D/PB/MC/Brusko_20803_TCRB.tsv"
    "${SLURM_SUBMIT_DIR}/BruskoT1D/PB/MC/Brusko_4276_TCRB.tsv"
    "${SLURM_SUBMIT_DIR}/BruskoT1D/PB/MC/Brusko_24036_TCRB.tsv"
    "${SLURM_SUBMIT_DIR}/BruskoT1D/PB/MC/Brusko_24273_TCRB.tsv"
    "${SLURM_SUBMIT_DIR}/BruskoT1D/PB/MC/Brusko_8469_TCRB.tsv"
)

# Copy over files to work dir
cd ${SCRATCH}
cp "${SLURM_SUBMIT_DIR}/TcrDataComparison/immuno_probs_config.ini" .
cp -r "${SLURM_SUBMIT_DIR}/human_TRB" .
cp ${filesArray[${SLURM_ARRAY_TASK_ID}]} .
BASENAME=`basename ${filesArray[${SLURM_ARRAY_TASK_ID}]}`

# Mark the output dir for automatic copying to $SLURM_SUBMIT_DIR afterwards
savefile claim_1

# Create the job dir for the output files
mkdir claim_1
cd claim_1
mkdir models
cd models
mkdir subject_${SLURM_ARRAY_TASK_ID}
cd subject_${SLURM_ARRAY_TASK_ID}


###########################
# Build model files
###########################
# Extract the data files from the Brusko data file
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'converted' -config-file '../../../immuno_probs_config.ini' convert -ref V '../../../human_TRB/TRBV.fasta' -ref J '../../../human_TRB/TRBJ.fasta' -seqs "../../../${BASENAME}"

# Create a model (params and marginals) for each data extract
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'unproductive' -config-file '../../../immuno_probs_config.ini' build -ref V '../../../human_TRB/TRBV.fasta' -ref D '../../../human_TRB/TRBD.fasta' -ref J '../../../human_TRB/TRBJ.fasta' -seqs 'converted_full_length_unproductive.tsv' -n-iter 10 -type 'beta'
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'productive' -config-file '../../../immuno_probs_config.ini' build -ref V '../../../human_TRB/TRBV.fasta' -ref D '../../../human_TRB/TRBD.fasta' -ref J '../../../human_TRB/TRBJ.fasta' -seqs 'converted_full_length_productive.tsv' -n-iter 10 -type 'beta'
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'all' -config-file '../../../immuno_probs_config.ini' build -ref V '../../../human_TRB/TRBV.fasta' -ref D '../../../human_TRB/TRBD.fasta' -ref J '../../../human_TRB/TRBJ.fasta' -seqs 'converted_full_length.tsv' -n-iter 10 -type 'beta'


# Exit succesfully
exit 0
