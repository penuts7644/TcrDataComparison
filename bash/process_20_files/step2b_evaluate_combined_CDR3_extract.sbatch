#!/bin/bash

#SBATCH --account=nn9603k

#SBATCH --job-name=ImmunoProbs

#SBATCH --array=0-19

#SBATCH --time=0-06:00:00

#SBATCH --nodes=20 --ntasks-per-node=1 --cpus-per-task=32

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

# Copy over files to work dir
cd ${SCRATCH}
cp "${SLURM_SUBMIT_DIR}/TcrDataComparison/immuno_probs_config.ini" .
cp "${SLURM_SUBMIT_DIR}/20files_combined_data_extract_and_models/sequence_data_extractor_CDR3.tsv" .
cp -r "${SLURM_SUBMIT_DIR}/human_TRB" .

# Mark the output dir for automatic copying to $SLURM_SUBMIT_DIR afterwards
savefile 20files_evaluated_CDR3

# Evaluate the sequences for given model in array
mkdir 20files_evaluated_CDR3
cd 20files_evaluated_CDR3
mkdir "evaluate_${SLURM_ARRAY_TASK_ID}"
cd "evaluate_${SLURM_ARRAY_TASK_ID}"
cp -r "${SLURM_SUBMIT_DIR}/20files_data_extract_and_models/immuno_probs_${SLURM_ARRAY_TASK_ID}/model" .
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'pgen_estimate_unproductive_CDR3' -config-file '../../immuno_probs_config.ini' evaluate -custom-model 'model/unproductive_params.txt' 'model/unproductive_marginals.txt' -seqs '../../sequence_data_extractor_CDR3.tsv' -type 'beta' -cdr3 -anchor V '../../human_TRB/V_gene_CDR3_anchors.tsv' -anchor J '../../human_TRB/J_gene_CDR3_anchors.tsv'
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'pgen_estimate_productive_CDR3' -config-file '../../immuno_probs_config.ini' evaluate -custom-model 'model/productive_params.txt' 'model/productive_marginals.txt' -seqs '../../sequence_data_extractor_CDR3.tsv' -type 'beta' -cdr3 -anchor V '../../human_TRB/V_gene_CDR3_anchors.tsv' -anchor J '../../human_TRB/J_gene_CDR3_anchors.tsv'
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'pgen_estimate_all_CDR3' -config-file '../../immuno_probs_config.ini' evaluate -custom-model 'model/all_params.txt' 'model/all_marginals.txt' -seqs '../../sequence_data_extractor_CDR3.tsv' -type 'beta' -cdr3 -anchor V '../../human_TRB/V_gene_CDR3_anchors.tsv' -anchor J '../../human_TRB/J_gene_CDR3_anchors.tsv'
rm -r model
cd ../

# Exit succesfully
exit 0
