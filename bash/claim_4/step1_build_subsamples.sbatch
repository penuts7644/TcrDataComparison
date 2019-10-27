#!/bin/bash

#SBATCH --account=nn9603k
#SBATCH --job-name=ImmunoProbs
#SBATCH --array=0-4
#SBATCH --time=1-00:00:00
#SBATCH --qos=preproc
#SBATCH --ntasks-per-node=1 --cpus-per-task=32


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

subsampleArrayText=(
    'all'
    '50000'
    '10000'
    '5000'
    '1000'
    '500'
    '100'
)

subsampleArrayNumber=(
    0
    50000
    10000
    5000
    1000
    500
    100
)

# Copy over files to work dir
cd ${SCRATCH}
cp "${SLURM_SUBMIT_DIR}/TcrDataComparison/immuno_probs_config.ini" .
cp -r "${SLURM_SUBMIT_DIR}/human_TRB" .
cp ${filesArray[${SLURM_ARRAY_TASK_ID}]} .
BASENAME=`basename ${filesArray[${SLURM_ARRAY_TASK_ID}]}`

# Mark the output dir for automatic copying to $SLURM_SUBMIT_DIR afterwards
savefile claim_4

# Create the job dir for the output files
mkdir claim_4
cd claim_4
mkdir models
cd models
mkdir subject_${SLURM_ARRAY_TASK_ID}
cd subject_${SLURM_ARRAY_TASK_ID}


###########################
# Build model files
###########################
# Make models for each subsample size
for i in "${!subsampleArrayText[@]}"
do
    mkdir "${subsampleArrayText[$i]}"
    cd "${subsampleArrayText[$i]}"

    # Extract the data files from data file and make model
    immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'converted' -config-file '../../../../immuno_probs_config.ini' convert -ref V '../../../../human_TRB/TRBV.fasta' -ref J '../../../../human_TRB/TRBJ.fasta' -seqs "../../../../${BASENAME}" -n-random "${subsampleArrayNumber[$i]}"

    # Create a model (params and marginals)
    immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'all' -config-file '../../../../immuno_probs_config.ini' build -ref V '../../../../human_TRB/TRBV.fasta' -ref D '../../../../human_TRB/TRBD.fasta' -ref J '../../../../human_TRB/TRBJ.fasta' -seqs 'converted_full_length.tsv' -n-iter 10 -type 'beta'
    cd ../
done


# Exit succesfully
exit 0
