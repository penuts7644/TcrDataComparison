#!/bin/bash

#SBATCH --account=nn9603k
#SBATCH --job-name=ImmunoProbs
#SBATCH --array=0-4
#SBATCH --time=0-04:00:00
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
    "${SLURM_SUBMIT_DIR}/Emerson2017/reduced_15e4/HIP05832.tsv"
    "${SLURM_SUBMIT_DIR}/Emerson2017/reduced_15e4/HIP09051.tsv"
    "${SLURM_SUBMIT_DIR}/Emerson2017/reduced_15e4/HIP13809.tsv"
    "${SLURM_SUBMIT_DIR}/Emerson2017/reduced_15e4/HIP08816.tsv"
    "${SLURM_SUBMIT_DIR}/Emerson2017/reduced_15e4/HIP08346.tsv"
)

# Copy over files to work dir
cd ${SCRATCH}
cp "${SLURM_SUBMIT_DIR}/TcrDataComparison/immuno_probs_config.ini" .
cp -r "${SLURM_SUBMIT_DIR}/human_TRB" .
cp ${filesArray[${SLURM_ARRAY_TASK_ID}]} .
BASENAME=`basename ${filesArray[${SLURM_ARRAY_TASK_ID}]}`

# Mark the output dir for automatic copying to $SLURM_SUBMIT_DIR afterwards
savefile claim_2

# Create the job dir for the output files
mkdir claim_2
cd claim_2
mkdir models
cd models
mkdir emerson
cd emerson
mkdir subject_${SLURM_ARRAY_TASK_ID}
cd subject_${SLURM_ARRAY_TASK_ID}


###########################
# Build model files
###########################
# Extract the data files from the data file
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'converted' -config-file '../../../../immuno_probs_config.ini' convert -ref V '../../../../human_TRB/TRBV.fasta' -ref J '../../../../human_TRB/TRBJ.fasta' -seqs "../../../../${BASENAME}"

# Create a model (params and marginals) for each data extract
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'all' -config-file '../../../../immuno_probs_config.ini' build -ref V '../../../../human_TRB/TRBV.fasta' -ref D '../../../../human_TRB/TRBD.fasta' -ref J '../../../../human_TRB/TRBJ.fasta' -seqs 'converted_full_length.tsv' -n-iter 10 -type 'beta'


# Exit succesfully
exit 0
