#!/bin/bash
#SBATCH --array=0-0
#SBATCH --time=01:00:00
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output=./slurm-log/slurm-%A_%a.out

#####
source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube


bowtie2-build --threads 16 /projects/li-lab/agarwa/CUBE/DeepTact/dataset/hg19.fa /projects/li-lab/agarwa/CUBE/DeepTact/dataset/hg19.fa.index