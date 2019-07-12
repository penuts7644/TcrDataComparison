#!/bin/bash

#SBATCH --account=nn9603k

#SBATCH --job-name=ImmunoProbs

#SBATCH --array=0-1585%200

#SBATCH --time=0-04:00:00

#SBATCH --nodes=1586 --ntasks-per-node=1 --cpus-per-task=32

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

filesArray=(
    "${SLURM_SUBMIT_DIR}/BruskoT1D/PB/MC/*"
)

# Copy over files to work dir
cd ${SCRATCH}
cp "${SLURM_SUBMIT_DIR}/TcrDataComparison/immuno_probs_config.ini" .
cp "${SLURM_SUBMIT_DIR}/TcrDataComparison/python/SequenceDataExtractor.py" .
cp -r "${SLURM_SUBMIT_DIR}/human_TRB" .
cp ${filesArray[${SLURM_ARRAY_TASK_ID}]} .

# Mark the output dir for automatic copying to $SLURM_SUBMIT_DIR afterwards
savefile all_data_extract_and_models

# Create the job dir for the output files
mkdir all_data_extract_and_models
cd all_data_extract_and_models
mkdir immuno_probs_${SLURM_ARRAY_TASK_ID}
cd immuno_probs_${SLURM_ARRAY_TASK_ID}

# Extract the data files from the Brusko data file
FILE_BASENAME=`basename ${filesArray[${SLURM_ARRAY_TASK_ID}]}`
python ../../SequenceDataExtractor.py --num-threads ${OMP_NUM_THREADS} "../../${FILE_BASENAME}" '../../human_TRB/TRBV.fasta' '../../human_TRB/TRBJ.fasta' '\t' &> 'sequence_extract_log.txt'

# Create a model (params and marginals) for each data extract
mkdir model
cd model
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'unproductive' -config-file '../../../immuno_probs_config.ini' build-igor-model -ref V '../../../human_TRB/TRBV.fasta' -ref D '../../../human_TRB/TRBD.fasta' -ref J '../../../human_TRB/TRBJ.fasta' -seqs '../sequence_data_extractor_unproductive.tsv' -n-iter 10 -type 'beta' &> 'unproductive_log.txt'
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'productive' -config-file '../../../immuno_probs_config.ini' build-igor-model -ref V '../../../human_TRB/TRBV.fasta' -ref D '../../../human_TRB/TRBD.fasta' -ref J '../../../human_TRB/TRBJ.fasta' -seqs '../sequence_data_extractor_productive.tsv' -n-iter 10 -type 'beta' &> 'productive_log.txt'
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'all' -config-file '../../../immuno_probs_config.ini' build-igor-model -ref V '../../../human_TRB/TRBV.fasta' -ref D '../../../human_TRB/TRBD.fasta' -ref J '../../../human_TRB/TRBJ.fasta' -seqs '../sequence_data_extractor_all.tsv' -n-iter 10 -type 'beta' &> 'all_log.txt'
cd ../

# Evaluate the CDR3 sequences with each model
mkdir evaluate
cd evaluate
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'pgen_estimate_unproductive_CDR3' -config-file '../../../immuno_probs_config.ini' evaluate-seqs -custom-model '../model/unproductive_params.txt' '../model/unproductive_marginals.txt' -seqs '../sequence_data_extractor_CDR3.tsv' -type 'beta' -cdr3 -anchor V '../../../human_TRB/V_gene_CDR3_anchors.tsv' -anchor J '../../../human_TRB/J_gene_CDR3_anchors.tsv' &> 'unproductive_log.txt'
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'pgen_estimate_productive_CDR3' -config-file '../../../immuno_probs_config.ini' evaluate-seqs -custom-model '../model/productive_params.txt' '../model/productive_marginals.txt' -seqs '../sequence_data_extractor_CDR3.tsv' -type 'beta' -cdr3 -anchor V '../../../human_TRB/V_gene_CDR3_anchors.tsv' -anchor J '../../../human_TRB/J_gene_CDR3_anchors.tsv' &> 'productive_log.txt'
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'pgen_estimate_all_CDR3' -config-file '../../../immuno_probs_config.ini' evaluate-seqs -custom-model '../model/all_params.txt' '../model/all_marginals.txt' -seqs '../sequence_data_extractor_CDR3.tsv' -type 'beta' -cdr3 -anchor V '../../../human_TRB/V_gene_CDR3_anchors.tsv' -anchor J '../../../human_TRB/J_gene_CDR3_anchors.tsv' &> 'all_log.txt'
cd ../

# Exit succesfully
exit 0
