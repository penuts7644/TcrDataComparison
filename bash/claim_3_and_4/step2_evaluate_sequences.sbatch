#!/bin/bash

#SBATCH --account=nn9603k
#SBATCH --job-name=ImmunoProbs
#SBATCH --array=0-2
#SBATCH --time=1-00:00:00
#SBATCH --partition=bigmem
#SBATCH --nodes=3 --ntasks-per-node=1 --cpus-per-task=32
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
    "${SLURM_SUBMIT_DIR}/BruskoT1D/subset/Brusko_20803_TCRB.tsv"
    "${SLURM_SUBMIT_DIR}/BruskoT1D/subset/Brusko_4276_TCRB.tsv"
    "${SLURM_SUBMIT_DIR}/BruskoT1D/subset/Brusko_24036_TCRB.tsv"
    "${SLURM_SUBMIT_DIR}/BruskoT1D/subset/Brusko_24273_TCRB.tsv"
    "${SLURM_SUBMIT_DIR}/BruskoT1D/subset/Brusko_8469_TCRB.tsv"
    "${SLURM_SUBMIT_DIR}/Emerson2017/HIP05595.tsv"
    "${SLURM_SUBMIT_DIR}/Emerson2017/HIP14071.tsv"
    "${SLURM_SUBMIT_DIR}/Emerson2017/HIP13865.tsv"
    "${SLURM_SUBMIT_DIR}/Emerson2017/HIP10376.tsv"
    "${SLURM_SUBMIT_DIR}/Emerson2017/HIP13810.tsv"
    "${SLURM_SUBMIT_DIR}/Peakman2018/reduced_18e4/T1D9_TN.tsv"
    "${SLURM_SUBMIT_DIR}/Peakman2018/reduced_18e4/T1D10_TN.tsv"
    "${SLURM_SUBMIT_DIR}/Peakman2018/reduced_18e4/T1D11_TN.tsv"
    "${SLURM_SUBMIT_DIR}/Peakman2018/reduced_18e4/HD8_TN.tsv"
    "${SLURM_SUBMIT_DIR}/Peakman2018/reduced_18e4/HD12_TN.tsv"
)

projectArray=(
    'brusko'
    'emerson'
    'peakman'
)

# Copy over files to work dir
cd ${SCRATCH}
cp "${SLURM_SUBMIT_DIR}/TcrDataComparison/immuno_probs_config.ini" .
cp -r "${SLURM_SUBMIT_DIR}/human_TRB" .

# Mark the output dir for automatic copying to $SLURM_SUBMIT_DIR afterwards
savefile claim_3

# Create the job dir for the output files
mkdir claim_3
cd claim_3
mkdir evaluations
cd evaluations
mkdir ${projectArray[$SLURM_ARRAY_TASK_ID]}
cd ${projectArray[$SLURM_ARRAY_TASK_ID]}


###########################
# Evaluate combined dataset
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

# Extract the data files from the Brusko data file
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'converted' -config-file '../../../immuno_probs_config.ini' convert -ref V '../../../human_TRB/TRBV.fasta' -ref J '../../../human_TRB/TRBJ.fasta' -seqs tmp_combined_dataset.tsv
rm tmp_combined_dataset.tsv

# In a loop copy over necessary model files and evaluate
cp -r "${SLURM_SUBMIT_DIR}/claim_3/models/${projectArray[$SLURM_ARRAY_TASK_ID]}" model
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'pgen_estimate_all_CDR3' -config-file '../../../immuno_probs_config.ini' evaluate -custom-model 'model/all_params.txt' 'model/all_marginals.txt' -seqs 'converted_CDR3.tsv' -type 'beta' -cdr3 -anchor V '../../../human_TRB/V_gene_CDR3_anchors.tsv' -anchor J '../../../human_TRB/J_gene_CDR3_anchors.tsv'
rm -r model


# Exit succesfully
exit 0
