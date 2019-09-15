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
NUMBER_OF_FILES=`ls "${SLURM_SUBMIT_DIR}/claim_2/models/emerson" | wc -l`
NUMBER_OF_FILES=`expr ${NUMBER_OF_FILES} - 1`
for i in $(seq 0 $NUMBER_OF_FILES)
do
    cp -r "${SLURM_SUBMIT_DIR}/claim_2/models/emerson/subject_${i}" .
    ALL_MODELS+="-model ${i} ../../../subject_${i}/all_params.txt ../../../subject_${i}/all_marginals.txt "
done

# Mark the output dir for automatic copying to $SLURM_SUBMIT_DIR afterwards
savefile claim_2

# Create the job dir for the output files
mkdir claim_2
cd claim_2
mkdir entropies
cd entropies
mkdir emerson
cd emerson


###########################
# Calculate model entropies
###########################
python ../../../model_processing/CalcModelEntropy.py --num-threads ${OMP_NUM_THREADS} ${ALL_MODELS}


# Exit succesfully
exit 0
