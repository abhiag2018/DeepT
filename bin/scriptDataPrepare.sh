#!/bin/bash
#SBATCH --array=0-0
#SBATCH --output=./slurm-log/slurm-%A_%a.out
#SBATCH --error=./slurm-log/slurm-%A.out
#SBATCH --qos batch 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time 12:00:00


#####
source ~/miniconda3/etc/profile.d/conda.sh
conda activate dt1

python ../DeepTACT-master/DataPrepare.py "/fastscratch/agarwa/DeepTact_tmp/TrainingData/$1" P-E