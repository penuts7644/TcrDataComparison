#!/bin/bash

#SBATCH --account=nn9603k
#SBATCH --job-name=ImmunoProbs
#SBATCH --time=1-00:00:00
#SBATCH --partition=bigmem
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=32
#SBATCH --mem-per-cpu=14G


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
    "${SLURM_SUBMIT_DIR}/Peakman2018/reduced_15e4/T1D9_CM.tsv"
    "${SLURM_SUBMIT_DIR}/Peakman2018/reduced_15e4/T1D10_CM.tsv"
    "${SLURM_SUBMIT_DIR}/Peakman2018/reduced_15e4/T1D11_CM.tsv"
    "${SLURM_SUBMIT_DIR}/Peakman2018/reduced_15e4/HD10_CM.tsv"
    "${SLURM_SUBMIT_DIR}/Peakman2018/reduced_15e4/HD14_CM.tsv"
)

cd ${SCRATCH}
cp "${SLURM_SUBMIT_DIR}/TcrDataComparison/immuno_probs_config.ini" .
cp -r "${SLURM_SUBMIT_DIR}/human_TRB" .

# Mark the output dir for automatic copying to $SLURM_SUBMIT_DIR afterwards
savefile claim_3_and_4

# Create the job dir for the output files
mkdir claim_3_and_4
cd claim_3_and_4
mkdir models
cd models
mkdir peakman
cd peakman


###########################
# Build model files
###########################
# Combine the data extracts from the previous models
for i in "${!filesArray[@]}"
do
    if [ $i == 0 ]; then
        cp "${filesArray[$i]}" tmp_combined_dataset.tsv
    else
        tail -n +2 -q "${filesArray[$i]}" >> tmp_combined_dataset.tsv
    fi
done

# Extract the data files from the peakman data file
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'converted' -config-file '../../../immuno_probs_config.ini' convert -ref V '../../../human_TRB/TRBV.fasta' -ref J '../../../human_TRB/TRBJ.fasta' -seqs tmp_combined_dataset.tsv
rm tmp_combined_dataset.tsv

# Create a model (params and marginals) for each data extract
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'all' -config-file '../../../immuno_probs_config.ini' build -ref V '../../../human_TRB/TRBV.fasta' -ref D '../../../human_TRB/TRBD.fasta' -ref J '../../../human_TRB/TRBJ.fasta' -seqs 'converted_full_length.tsv' -n-iter 10 -type 'beta'


# Exit succesfully
exit 0
