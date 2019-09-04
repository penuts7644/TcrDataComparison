#!/bin/bash

#SBATCH --account=nn9603k
#SBATCH --job-name=ImmunoProbs
#SBATCH --array=0-4
#SBATCH --time=1-00:00:00
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
    "${SLURM_SUBMIT_DIR}/BruskoT1D/subset/Brusko_20803_TCRB.tsv"
    "${SLURM_SUBMIT_DIR}/BruskoT1D/subset/Brusko_4276_TCRB.tsv"
    "${SLURM_SUBMIT_DIR}/BruskoT1D/subset/Brusko_24036_TCRB.tsv"
    "${SLURM_SUBMIT_DIR}/BruskoT1D/subset/Brusko_24273_TCRB.tsv"
    "${SLURM_SUBMIT_DIR}/BruskoT1D/subset/Brusko_8469_TCRB.tsv"
)

subsampleArrayText=(
    'all'
    '100000'
    '50000'
    '10000'
    '5000'
    '1000'
    '500'
)

# Copy over files to work dir
cd ${SCRATCH}
cp "${SLURM_SUBMIT_DIR}/TcrDataComparison/immuno_probs_config.ini" .
cp -r "${SLURM_SUBMIT_DIR}/human_TRB" .

# Mark the output dir for automatic copying to $SLURM_SUBMIT_DIR afterwards
savefile claim_5

# Create the job dir for the output files
mkdir claim_5
cd claim_5
mkdir evaluations
cd evaluations
mkdir subject_${SLURM_ARRAY_TASK_ID}
cd subject_${SLURM_ARRAY_TASK_ID}


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

# Make models for each subsample size
for i in "${!subsampleArrayText[@]}"
do
    mkdir "${subsampleArrayText[$i]}"
    cd "${subsampleArrayText[$i]}"

    # Extract the data files from the Brusko data file
    immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'converted' -config-file '../../../../immuno_probs_config.ini' convert -ref V '../../../../human_TRB/TRBV.fasta' -ref J '../../../../human_TRB/TRBJ.fasta' -seqs '../tmp_combined_dataset.tsv'

    # Evaluate the CDR3 sequences with each model
    cp -r "${SLURM_SUBMIT_DIR}/claim_5/models/subject_${SLURM_ARRAY_TASK_ID}/${subsampleArrayText[$i]}" model
    immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'pgen_estimate_all_CDR3' -config-file '../../../../immuno_probs_config.ini' evaluate -custom-model 'model/all_params.txt' 'model/all_marginals.txt' -seqs 'converted_CDR3.tsv' -type 'beta' -cdr3 -anchor V '../../../../human_TRB/V_gene_CDR3_anchors.tsv' -anchor J '../../../../human_TRB/J_gene_CDR3_anchors.tsv'
    rm -r model
    cd ../
done
rm tmp_combined_dataset.tsv


# Exit succesfully
exit 0
