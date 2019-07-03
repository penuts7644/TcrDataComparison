#!/bin/bash

#SBATCH --account=nn9603k

#SBATCH --job-name=ImmunoProbs

#SBATCH --time=1-00:00:00

#SBATCH --qos=preproc

#SBATCH --ntasks-per-node=1 --cpus-per-task=32

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
cp "${SLURM_SUBMIT_DIR}/20files_data_extract_and_models/immuno_probs_combined/sequence_data_extractor_CDR3.tsv" .
cp -r "${SLURM_SUBMIT_DIR}/human_TRB" .

# Mark the output dir for automatic copying to $SLURM_SUBMIT_DIR afterwards
savefile 20files_evaluated_CDR3

# Calculate entropy for all the models
mkdir 20files_evaluated_CDR3
cd 20files_evaluated_CDR3
cp -r "${SLURM_SUBMIT_DIR}/20files_data_extract_and_models/immuno_probs_combined/evaluate" .
mv evaluate evaluate_combined

# Evaluate using the build-in igor model
mkdir evaluate_igor
cd evaluate_igor
cp -r "${SLURM_SUBMIT_DIR}/igor_human_TCRB_model" .
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'pgen_estimate_CDR3' evaluate-seqs -custom-model 'igor_human_TCRB_model/model_params.txt' 'igor_human_TCRB_model/model_marginals.txt' -seqs '../../sequence_data_extractor_CDR3.tsv' -type 'beta' -cdr3 -anchor V '../../human_TRB/V_gene_CDR3_anchors.tsv' -anchor J '../../human_TRB/J_gene_CDR3_anchors.tsv' &> 'unproductive_log.txt'
rm -r igor_human_TCRB_model
cd ../

# In a loop copy over necessary model files and create command string
NUMBER_OF_FILES=`ls | wc -l`
NUMBER_OF_FILES= expr ${NUMBER_OF_FILES} - 2
for i in {0..${NUMBER_OF_FILES}}
do
    mkdir "evaluate_${i}"
    cd "evaluate_${i}"
    cp -r "${SLURM_SUBMIT_DIR}/20files_data_extract_and_models/immuno_probs_${i}/model" .
    immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'pgen_estimate_unproductive_CDR3' evaluate-seqs -custom-model 'model/unproductive_params.txt' 'model/unproductive_marginals.txt' -seqs '../../sequence_data_extractor_CDR3.tsv' -type 'beta' -cdr3 -anchor V '../../human_TRB/V_gene_CDR3_anchors.tsv' -anchor J '../../human_TRB/J_gene_CDR3_anchors.tsv' &> 'unproductive_log.txt'
    immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'pgen_estimate_productive_CDR3' evaluate-seqs -custom-model 'model/productive_params.txt' 'model/productive_marginals.txt' -seqs '../../sequence_data_extractor_CDR3.tsv' -type 'beta' -cdr3 -anchor V '../../human_TRB/V_gene_CDR3_anchors.tsv' -anchor J '../../human_TRB/J_gene_CDR3_anchors.tsv' &> 'productive_log.txt'
    immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'pgen_estimate_all_CDR3' evaluate-seqs -custom-model 'model/all_params.txt' 'model/all_marginals.txt' -seqs '../../sequence_data_extractor_CDR3.tsv' -type 'beta' -cdr3 -anchor V '../../human_TRB/V_gene_CDR3_anchors.tsv' -anchor J '../../human_TRB/J_gene_CDR3_anchors.tsv' &> 'all_log.txt'
    rm -r model
    cd ../
done

# Exit succesfully
exit 0
