#!/bin/bash
#SBATCH --array=0-13
#SBATCH -q training
#SBATCH -p gpu
#SBATCH --time=13-00:00:00
#SBATCH --gres=gpu:1
#SBATCH --mem=100G
#SBATCH --output=./slurm-log/slurm-%A_%a.out


source ~/miniconda3/etc/profile.d/conda.sh
conda activate /projects/li-lab/agarwa/conda_envs/cube


base_dir="/projects/li-lab/agarwa/CUBE/DeepTact/dataset/DeepTact_tmp.1/TrainingData"

cell=$1
num_rep=$2
python DeepTact_0.py ${base_dir}/${cell} P-E ${num_rep}
# eval_cell=$3
# python DeepTact_0.py ${base_dir}/${cell} P-E ${num_rep} ${eval_cell}

##

# sbatch scriptGPU.sh -a 0-<NUM_ENSEMBL-1> <CELL> <NUM_REP>