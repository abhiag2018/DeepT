#!/bin/bash
#SBATCH --array=0-3
#SBATCH -q training
#SBATCH -p gpu
#SBATCH --gres gpu:4
#SBATCH --time=13-00:00:00
#SBATCH --mem-per-gpu=24G
#SBATCH --output=./slurm-log/slurm-%A_%a.out


source ~/miniconda3/etc/profile.d/conda.sh
conda activate /projects/li-lab/agarwa/conda_envs/cube

# python ../DeepTACT-master/DeepTACT.py /projects/li-lab/agarwa/CUBE/DeepTact/dataset/DeepTact_tmp.1/train_Mon P-E 1 1

base_dir="/projects/li-lab/agarwa/CUBE/DeepTact/dataset/DeepTact_tmp.1/TrainingData"

cell=$1
num_rep=$2
# THEANO_FLAGS=device=cuda python ../DeepTACT-master/DeepTACT.py ${base_dir}/${cell} P-E ${num_rep} 2
python DeepTact_0.py ${base_dir}/${cell} P-E ${num_rep} $SLURM_ARRAY_TASK_COUNT $SLURM_ARRAY_TASK_ID
# python testGPU.py



# sbatch scriptGPU.sh <CELL> <NUM_REP>