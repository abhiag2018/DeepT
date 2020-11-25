#!/bin/bash
#SBATCH --array=0-0
#SBATCH --output=./slurm-log/slurm-%A_%a.out
#SBATCH --error=./slurm-log/slurm-%A.out
#SBATCH --cpus-per-task=32
#SBATCH --qos batch 
#SBATCH --mem=500G
#SBATCH --time 3-00:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

nfoptions=$1

nextflow run atacseq/main.nf --genome GRCh37 --profile docker --input design.csv $nfoptions