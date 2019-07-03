#!/bin/bash

#SBATCH --account=nn9603k

#SBATCH --job-name=ImmunoProbs

#SBATCH --time=0-02:00:00

#SBATCH --qos=preproc

#SBATCH --ntasks-per-node=1 --cpus-per-task=32

# Set up job environment, load modules and setup python
module restore system
module load Python/2.7.15-intel-2018b

# Set safer defaults for bash
set -o errexit
set -o nounset


# Copy over files to work dir
cd ${SCRATCH}
cp -r "${SLURM_SUBMIT_DIR}/TcrDataComparison/python/model_processing" .
cp -r "${SLURM_SUBMIT_DIR}/20files_data_extract_and_models/immuno_probs_combined/model" .
mv model model_combined
cp -r "${SLURM_SUBMIT_DIR}/igor_human_TCRB_model" .
mv igor_human_TCRB_model model_igor

# In a loop copy over necessary model files and create command string (already include the combined and igor models)
UNPRODUCTIVE_MODELS="-model igor_default ../model_igor/model_params.txt ../model_igor/model_marginals.txt -model combined_unproductive ../model_combined/unproductive_params.txt ../model_combined/unproductive_marginals.txt"
PRODUCTIVE_MODELS="-model igor_default ../model_igor/model_params.txt ../model_igor/model_marginals.txt -model combined_productive ../model_combined/productive_params.txt ../model_combined/productive_marginals.txt"
ALL_MODELS="-model igor_default ../model_igor/model_params.txt ../model_igor/model_marginals.txt -model combined_all ../model_combined/all_params.txt ../model_combined/all_marginals.txt"
NUMBER_OF_FILES=`ls | wc -l`
NUMBER_OF_FILES= expr ${NUMBER_OF_FILES} - 2
for i in {0..${NUMBER_OF_FILES}}
do
    cp -r "${SLURM_SUBMIT_DIR}/20files_data_extract_and_models/immuno_probs_${i}/model" .
    mv model "model_${i}"
    UNPRODUCTIVE_MODELS+=" -model ${i}_unproductive ../model_${i}/unproductive_params.txt ../model_${i}/unproductive_marginals.txt"
    PRODUCTIVE_MODELS+=" -model ${i}_productive ../model_${i}/productive_params.txt ../model_${i}/productive_marginals.txt"
    ALL_MODELS+=" -model ${i}_all ../model_${i}/all_params.txt ../model_${i}/all_marginals.txt"
done

# Mark the output dir for automatic copying to $SLURM_SUBMIT_DIR afterwards
savefile 20files_model_entropies

# Calculate entropy for all the models
mkdir 20files_model_entropies
cd 20files_model_entropies
python ../model_processing/CalcModelEntropy.py --num-threads ${OMP_NUM_THREADS} ${UNPRODUCTIVE_MODELS} &> 'unproductive_model_entropy_log.txt'
mv calc_model_entropy.tsv unproductive_calc_model_entropy.tsv
python ../model_processing/CalcModelEntropy.py --num-threads ${OMP_NUM_THREADS} ${PRODUCTIVE_MODELS} &> 'productive_model_entropy_log.txt'
mv calc_model_entropy.tsv productive_calc_model_entropy.tsv
python ../model_processing/CalcModelEntropy.py --num-threads ${OMP_NUM_THREADS} ${ALL_MODELS} &> 'all_model_entropy_log.txt'
mv calc_model_entropy.tsv all_calc_model_entropy.tsv
python ../model_processing/CalcModelEntropy.py --num-threads ${OMP_NUM_THREADS} ${UNPRODUCTIVE_MODELS} ${PRODUCTIVE_MODELS} ${ALL_MODELS} &> 'model_entropy_log.txt'
cd ../

# Exit succesfully
exit 0