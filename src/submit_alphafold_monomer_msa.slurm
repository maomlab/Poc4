#!/bin/bash
#SBATCH --job-name=af2_monomer_msa
#SBATCH --output=intermediate_data/parafold/logs/%x_%j.out
#SBATCH --error=intermediate_data/parafold/logs/%x_%j.err


if [ -z ${INPUT_PATH+x} ]; then
    echo "ERROR: The \${INPUT_PATH} variable is not set"
    echo "ERROR: Please pass it into the slum command like this:"
    echo "ERROR:    sbatch ... --export=INPUT_PATH='<input_path>',OUTPUT_PATH='<output_path>' ..."
fi

if [ -z ${OUTPUT_PATH+x} ]; then
    echo "ERROR: The \${OUTPUT_PATH} variable is not set"
    echo "ERROR: Please pass it into the slum command like this:"
    echo "ERROR:    sbatch ... --export=INPUT_PATH='<input_path>',OUTPUT_PATH='<output_path>' ..."
fi



#sbatch \
#  --account=tromeara99 \
#  --partition=standard \
#  --array=1-3 \
#  --ntasks-per-node=1 \
#  --cpus-per-task=8 \
#  --time=10:00:00 \
#  --mem-per-cpu=6GB \
#  --output=intermediate_data/parafold/Frame2Seq_fixinterface_T1.0/logs/%x_%j.out \
#  --error=intermediate_data/parafold/Frame2Seq_fixinterface_T1.0/logs/%x_%j.err \
#  --export=\
#    INPUT_PATH='intermediate_data/parafold/Frame2Seq_fixinterface_T1.0/input',\
#    OUTPUT_PATH='intermediate_data/parafold/Frame2Seq_fixinterface_T1.0/output' \
#  src/submit_alphafold_monomer_msa.slurm

# This script is run once for each SLURM_ARRAY_TASK_ID,
# it looks up the structure id in the intermediate_data/parafold/msa.todo
# and submits it to run_alphafold.sh with the "-f" flag to only generate
# the features.pkl e.g. only compute the multiple sequence alignments
# This is done first because it is largely CPU bound while the structure
# prediction is done second and is largely GPU bound.


module purge
module restore cauris
module list
# e.g.
# 1) openjdk/18.0.1.1   2) gcc/13.2.0   3) image-libraries/220318   4) R/4.4.0   5) cuda/11.3.0   6) cudnn/11.3-v8.2.1

export PATH=${HOME}/miniconda3/envs/parafold/bin/:${PATH}:

# install parafold using the precise instructions https://github.com/Zuricho/ParallelFold
${HOME}/miniconda/bin/activate parafold

CONDA_ENV="parafold"
PARAFOLD_PATH="${HOME}/turbo/opt/ParallelFold"
STRUCTURE_ID=$(cat ${OUTPUT_PATH}/msa.todo | sed "${SLURM_ARRAY_TASK_ID}q;d")
INPUT_FASTA=${INPUT_PATH}/${STRUCTURE_ID}/input.fasta
OUTPUT_STRUCTURE_PATH=${OUTPUT_PATH}/${STRUCTURE_ID}


echo "currentw working directory: $(pwd)"
echo "CONDA_ENV: ${CONDA_ENV}"
echo "PARAFOLD_PATH: ${PARAFOLD_PATH}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "INTPUT_PATH: ${INPUT_PATH}"
echo "OUTPUT_PATH: ${OUTPUT_PATH}"
echo "STRUCTURE_ID: ${STRUCTURE_ID}"
echo "INPUT_FASTA: ${INPUT_FASTA}"
echo "OUTPUT_STRUCTURE_PATH: ${OUTPUT_STRUCTURE_PATH}"

echo "Creating output path..."
mkdir -p ${OUTPUT_STRUCTURE_PATH}


${PARAFOLD_PATH}/run_alphafold.sh \
  -d ${PARAFOLD_PATH}/data \
  -i ${INPUT_FASTA} \
  -o ${OUTPUT_STRUCTURE_PATH} \
  -p monomer \
  -f

# Check if the job is complete and write that to the output folder
if [ -f ${OUTPUT_STRUCTURE_PATH}/input/features.pkl ] || \
   [ -f ${OUTPUT_STRUCTURE_PATH}/input/msas/bfd_uniref_hits.a3m ] || \
   [ -f ${OUTPUT_STRUCTURE_PATH}/input/msas/mgnify_hits.sto     ] || \
   [ -f ${OUTPUT_STRUCTURE_PATH}/input/msas/pdb_hits.sto        ] || \
   [ -f ${OUTPUT_STRUCTURE_PATH}/input/msas/uniref90_hits.sto   ]; then
   echo "DONE"
   echo -e "job_id\n$SLURM_JOB_ID" > "${OUTPUT_PATH}/msa.done"
fi


   
