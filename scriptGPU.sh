#!/bin/bash
#SBATCH -q training
#SBATCH -p gpu
#SBATCH --time=13-00:00:00
#SBATCH --gres=gpu:v100:1
#SBATCH --mem=100G



source ~/miniconda3/etc/profile.d/conda.sh
conda activate /projects/li-lab/agarwa/conda_envs/cube

PYTHONPATH=NN-feature-prep/bin:$PYTHONPATH
base_dir="/projects/li-lab/agarwa/CUBE/DeepTact/code/storeDir"

job=$1
cell=$2
num_rep=$3
# python preprocessing/bin/DeepTact_0.py ${base_dir}/${cell} P-E ${num_rep}
eval_cell=$4
python DeepTact_0.py $job ${base_dir}/${cell} P-E ${num_rep} ${eval_cell}

##

# sbatch scriptGPU.sh test tB 3 tB
# sbatch scriptGPU.sh -a 0-<NUM_ENSEMBL-1> <CELL> <NUM_REP>
