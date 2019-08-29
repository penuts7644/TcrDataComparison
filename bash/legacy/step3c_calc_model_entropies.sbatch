#!/bin/bash

#SBATCH --account=nn9603k

#SBATCH --job-name=ImmunoProbs

#SBATCH --time=0-00:30:00

#SBATCH --qos=preproc

#SBATCH --ntasks-per-node=1 --cpus-per-task=32

# Set up job environment, load modules and setup python
module restore system
module load Python/2.7.15-intel-2018b

# Set safer defaults for bash
set -o errexit
set -o nounset

# Change this value to '100000', '50000', '10000', '5000', '1000' or '500'
SUBSET_ID='all'

# Copy over files to work dir
cd ${SCRATCH}
cp -r "${SLURM_SUBMIT_DIR}/TcrDataComparison/python/model_processing" .
cp -r "${SLURM_SUBMIT_DIR}/process_5_files/${SUBSET_ID}/combined_data_extract_and_models/model" .
mv model model_combined

# In a loop copy over necessary model files and create command string (already include the combined models)
UNPRODUCTIVE_MODELS="-model combined ../../../../model_combined/unproductive_params.txt ../../../../model_combined/unproductive_marginals.txt"
PRODUCTIVE_MODELS="-model combined ../../../../model_combined/productive_params.txt ../../../../model_combined/productive_marginals.txt"
ALL_MODELS="-model combined ../../../../model_combined/all_params.txt ../../../../model_combined/all_marginals.txt"
NUMBER_OF_FILES=`ls "${SLURM_SUBMIT_DIR}/process_5_files/${SUBSET_ID}/data_extract_and_models" | wc -l`
NUMBER_OF_FILES=`expr ${NUMBER_OF_FILES} - 1`
for i in $(seq 0 $NUMBER_OF_FILES)
do
    cp -r "${SLURM_SUBMIT_DIR}/process_5_files/${SUBSET_ID}/data_extract_and_models/immuno_probs_${i}/model" .
    mv model "model_${i}"
    UNPRODUCTIVE_MODELS+=" -model ${i} ../../../../model_${i}/unproductive_params.txt ../../../../model_${i}/unproductive_marginals.txt"
    PRODUCTIVE_MODELS+=" -model ${i} ../../../../model_${i}/productive_params.txt ../../../../model_${i}/productive_marginals.txt"
    ALL_MODELS+=" -model ${i} ../../../../model_${i}/all_params.txt ../../../../model_${i}/all_marginals.txt"
done

# Mark the output dir for automatic copying to $SLURM_SUBMIT_DIR afterwards
savefile process_5_files

# Calculate entropy for all the models
mkdir process_5_files
cd process_5_files
mkdir ${SUBSET_ID}
cd ${SUBSET_ID}
mkdir model_entropies
cd model_entropies
mkdir unproductive
cd unproductive
python ../../../../model_processing/CalcModelEntropy.py --num-threads ${OMP_NUM_THREADS} ${UNPRODUCTIVE_MODELS}
cd ../
mkdir productive
cd productive
python ../../../../model_processing/CalcModelEntropy.py --num-threads ${OMP_NUM_THREADS} ${PRODUCTIVE_MODELS}
cd ../
mkdir all
cd all
python ../../../../model_processing/CalcModelEntropy.py --num-threads ${OMP_NUM_THREADS} ${ALL_MODELS}
cd ../../

# Exit succesfully
exit 0
