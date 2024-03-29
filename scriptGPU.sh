#!/bin/bash
#SBATCH -q training
#SBATCH -p gpu
#SBATCH --time=13-00:00:00
#SBATCH --gres=gpu:v100:1
#SBATCH --mem=100G



source ~/miniconda3/etc/profile.d/conda.sh
conda activate /projects/li-lab/agarwa/conda_envs/cube

PYTHONPATH=NN-feature-prep/bin:$PYTHONPATH
# base_dir="/projects/li-lab/agarwa/CUBE/DeepTact/code/storeDir"
base_dir="/projects/li-lab/agarwa/CUBE/DeepTact/code/storeDir/data-all"

job=$1; shift
cell=$1; shift
num_rep=$1; shift
eval_cell=$1; shift

if [ $job == "train" ] && [ "$#" -ne 2 ]; then
    echo "2 arguments needed : train <celltype> <numRep>"
fi

if [ $job == "split" ] && [ "$#" -ne 2 ]; then
    echo "2 arguments needed : split <celltype> <numRep>"
fi

if [ $job == "test" ] && [ "$#" -lt 5 ]; then
    echo ">=5 arguments needed : test <celltype> <numRep> <evalCell> <appendStr> <bootstrap1> <bootstrap2>.."
fi

python DeepTact_0.py $job ${base_dir}/${cell} P-E ${num_rep} ${eval_cell} "$@"
