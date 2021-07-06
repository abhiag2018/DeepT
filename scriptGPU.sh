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

job=$1; shift
cell=$1; shift
num_rep=$1; shift
eval_cell=$1; shift
python DeepTact_0.py $job ${base_dir}/${cell} P-E ${num_rep} ${eval_cell} "$@"
