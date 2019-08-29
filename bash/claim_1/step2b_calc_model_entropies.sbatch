#!/bin/bash

#SBATCH --account=nn9603k
#SBATCH --job-name=ImmunoProbs
#SBATCH --time=0-00:30:00
#SBATCH --qos=preproc
#SBATCH --ntasks-per-node=1 --cpus-per-task=32


###########################
# Load necessary modules
###########################
module restore system
module load Python/2.7.15-intel-2018b
set -o errexit
set -o nounset


###########################
# Specify job variables
###########################
cd ${SCRATCH}
cp -r "${SLURM_SUBMIT_DIR}/TcrDataComparison/python/model_processing" .

# In a loop copy over necessary model files and create command string (already include the combined models)
UNPRODUCTIVE_MODELS=""
PRODUCTIVE_MODELS=""
ALL_MODELS=""
NUMBER_OF_FILES=`ls "${SLURM_SUBMIT_DIR}/claim_1/models" | wc -l`
NUMBER_OF_FILES=`expr ${NUMBER_OF_FILES} - 1`
for i in $(seq 0 $NUMBER_OF_FILES)
do
    cp -r "${SLURM_SUBMIT_DIR}/claim_1/models/subject_${i}/model" .
    mv model "model_${i}"
    UNPRODUCTIVE_MODELS+="-model ${i} ../../../model_${i}/unproductive_params.txt ../../../model_${i}/unproductive_marginals.txt "
    PRODUCTIVE_MODELS+="-model ${i} ../../../model_${i}/productive_params.txt ../../../model_${i}/productive_marginals.txt "
    ALL_MODELS+="-model ${i} ../../../model_${i}/all_params.txt ../../../model_${i}/all_marginals.txt "
done

# Mark the output dir for automatic copying to $SLURM_SUBMIT_DIR afterwards
savefile claim_1

# Create the job dir for the output files
mkdir claim_1
cd claim_1
mkdir entropies
cd entropies


###########################
# Calculate model entropies
###########################
mkdir unproductive
cd unproductive
python ../../../model_processing/CalcModelEntropy.py --num-threads ${OMP_NUM_THREADS} ${UNPRODUCTIVE_MODELS}
cd ../
mkdir productive
cd productive
python ../../../model_processing/CalcModelEntropy.py --num-threads ${OMP_NUM_THREADS} ${PRODUCTIVE_MODELS}
cd ../
mkdir all
cd all
python ../../../model_processing/CalcModelEntropy.py --num-threads ${OMP_NUM_THREADS} ${ALL_MODELS}
cd ../


# Exit succesfully
exit 0
