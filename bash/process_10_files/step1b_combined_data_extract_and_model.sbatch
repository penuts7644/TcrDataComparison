#!/bin/bash

#SBATCH --account=nn9603k

#SBATCH --job-name=ImmunoProbs

#SBATCH --time=7-00:00:00

#SBATCH --partition=bigmem

#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=32

#SBATCH --mem-per-cpu=16G

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
cp "${SLURM_SUBMIT_DIR}/TcrDataComparison/python/SequenceDataExtractor.py" .
cp -r "${SLURM_SUBMIT_DIR}/human_TRB" .
mkdir data
cp "${SLURM_SUBMIT_DIR}/BruskoT1D/PB/MC/Brusko_13524_TCRB.tsv" data/.
cp "${SLURM_SUBMIT_DIR}/BruskoT1D/PB/MC/Brusko_9374_TCRB.tsv" data/.
cp "${SLURM_SUBMIT_DIR}/BruskoT1D/PB/MC/Brusko_22340_TCRB.tsv" data/.
cp "${SLURM_SUBMIT_DIR}/BruskoT1D/PB/MC/Brusko_29173_TCRB.tsv" data/.
cp "${SLURM_SUBMIT_DIR}/BruskoT1D/PB/MC/Brusko_7122_TCRB.tsv" data/.
cp "${SLURM_SUBMIT_DIR}/BruskoT1D/PB/MC/Brusko_17842_TCRB.tsv" data/.
cp "${SLURM_SUBMIT_DIR}/BruskoT1D/PB/MC/Brusko_26184_TCRB.tsv" data/.
cp "${SLURM_SUBMIT_DIR}/BruskoT1D/PB/MC/Brusko_4162_TCRB.tsv" data/.
cp "${SLURM_SUBMIT_DIR}/BruskoT1D/PB/MC/Brusko_10571_TCRB.tsv" data/.
cp "${SLURM_SUBMIT_DIR}/BruskoT1D/PB/MC/Brusko_26852_TCRB.tsv" data/.

# Mark the output dir for automatic copying to $SLURM_SUBMIT_DIR afterwards
savefile 10files_combined_data_extract_and_models

# Create the job dir for the output files
mkdir 10files_combined_data_extract_and_models
cd 10files_combined_data_extract_and_models

# Extract the data files from the Brusko data file
python ../SequenceDataExtractor.py --num-threads ${OMP_NUM_THREADS} "../data" '../human_TRB/TRBV.fasta' '../human_TRB/TRBJ.fasta' '\t' &> 'sequence_extract_log.txt'

# Create a model (params and marginals) for each data extract
mkdir model
cd model
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'unproductive' build-igor-model -ref V '../../human_TRB/TRBV.fasta' -ref D '../../human_TRB/TRBD.fasta' -ref J '../../human_TRB/TRBJ.fasta' -seqs '../sequence_data_extractor_unproductive.tsv' -n-iter 10 -type 'beta' &> 'unproductive_log.txt'
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'productive' build-igor-model -ref V '../../human_TRB/TRBV.fasta' -ref D '../../human_TRB/TRBD.fasta' -ref J '../../human_TRB/TRBJ.fasta' -seqs '../sequence_data_extractor_productive.tsv' -n-iter 10 -type 'beta' &> 'productive_log.txt'
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'all' build-igor-model -ref V '../../human_TRB/TRBV.fasta' -ref D '../../human_TRB/TRBD.fasta' -ref J '../../human_TRB/TRBJ.fasta' -seqs '../sequence_data_extractor_all.tsv' -n-iter 10 -type 'beta' &> 'all_log.txt'
cd ../

# Evaluate the CDR3 sequences with each model
mkdir evaluate
cd evaluate
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'pgen_estimate_unproductive_CDR3' evaluate-seqs -custom-model '../model/unproductive_params.txt' '../model/unproductive_marginals.txt' -seqs '../sequence_data_extractor_CDR3.tsv' -type 'beta' -cdr3 -anchor V '../../human_TRB/V_gene_CDR3_anchors.tsv' -anchor J '../../human_TRB/J_gene_CDR3_anchors.tsv' &> 'unproductive_log.txt'
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'pgen_estimate_productive_CDR3' evaluate-seqs -custom-model '../model/productive_params.txt' '../model/productive_marginals.txt' -seqs '../sequence_data_extractor_CDR3.tsv' -type 'beta' -cdr3 -anchor V '../../human_TRB/V_gene_CDR3_anchors.tsv' -anchor J '../../human_TRB/J_gene_CDR3_anchors.tsv' &> 'productive_log.txt'
immuno-probs -threads ${OMP_NUM_THREADS} -out-name 'pgen_estimate_all_CDR3' evaluate-seqs -custom-model '../model/all_params.txt' '../model/all_marginals.txt' -seqs '../sequence_data_extractor_CDR3.tsv' -type 'beta' -cdr3 -anchor V '../../human_TRB/V_gene_CDR3_anchors.tsv' -anchor J '../../human_TRB/J_gene_CDR3_anchors.tsv' &> 'all_log.txt'
cd ../

# Exit succesfully
exit 0
