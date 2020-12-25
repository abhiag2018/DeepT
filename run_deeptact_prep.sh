#!/bin/bash
#SBATCH --array=0-0
#SBATCH --output=./slurm-log/slurm-%A_%a.out
#SBATCH --error=./slurm-log/slurm-%A.out
#SBATCH --mem=100G
#SBATCH --qos batch 
#SBATCH --time 3-00:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

###
NXF_OPTS='-Xms16g -Xmx64g'

resumeDir=$1
nextflow run preprocessing/main.nf -profile slurm -w "/fastscratch/agarwa/nf-tmp/work" -resume $resumeDir
# nextflow -C preprocessing/nextflow-test.config run preprocessing/main.nf -profile slurm -w "/fastscratch/agarwa/nf-tmp/work"

# nextflow run preprocessing/main_1.nf -profile slurm -w "/fastscratch/agarwa/nf-tmp/work" -resume $resumeDir