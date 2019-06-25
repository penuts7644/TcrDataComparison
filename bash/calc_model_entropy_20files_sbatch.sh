#!/bin/bash

#SBATCH --account=nn9603k

#SBATCH --job-name=ImmunoProbs

#SBATCH --time=0-00:20:00

#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=4

# Set up job environment, load modules and setup python
module restore system
module load Python/2.7.15-intel-2018b

# Set safer defaults for bash
set -o errexit
set -o nounset


# Copy over files to work dir
cd ${SCRATCH}
cp -r ${SLURM_SUBMIT_DIR}/TcrDataComparison/python/model_processing .

# In a loop copy over necessary model files and create command string
MODELS_VARIABLE=""
for i in {0..19}
do
    mkdir "model_${i}"
    cp "${SLURM_SUBMIT_DIR}/immuno_probs_20files_models/immuno_probs_${i}/model/." "model_${i}/."
    MODELS_VARIABLE+=" -model ${i}_unproductive model_${i}/unproductive_params.txt model_${i}/unproductive_marginals.txt"
    MODELS_VARIABLE+=" -model ${i}_productive model_${i}/productive_params.txt model_${i}/productive_marginals.txt"
    MODELS_VARIABLE+=" -model ${i}_all model_${i}/all_params.txt model_${i}/all_marginals.txt"
done

# Mark the output dir for automatic copying to $SLURM_SUBMIT_DIR afterwards
savefile model_processing/calc_model_entropy.tsv model_processing/model_entropy_log.txt

# Calculate entropy for all the models
cd model_processing
python CalcModelEntropy.py --num-threads ${OMP_NUM_THREADS} ${MODELS_VARIABLE} &> 'model_entropy_log.txt'
cd ../

# Exit succesfully
exit 0
