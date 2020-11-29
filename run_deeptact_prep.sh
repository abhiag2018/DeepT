#!/bin/bash
#SBATCH --array=0-0
#SBATCH --output=./slurm-log/slurm-%A_%a.out
#SBATCH --error=./slurm-log/slurm-%A.out
#SBATCH --qos batch 
#SBATCH --time 3-00:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

###
nextflow run main.sh -resume