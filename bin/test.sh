#!/bin/bash
#SBATCH --array=0-0
#SBATCH --time=6:00:00
#SBATCH --mem=50G
#SBATCH --output=./slurm-log/slurm-%A_%a.out
#SBATCH --ntasks=128

python multiprocessing_copy.py