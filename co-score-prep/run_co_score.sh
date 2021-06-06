#!/bin/bash
#SBATCH --mem=10G
#SBATCH --qos batch 
#SBATCH --time 6:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

###
NXF_OPTS='-Xms16g -Xmx64g'

resumeDir=$1
nextflow  -c nextflow.config \
	-c  params.config \
	-c ../pr-enh-prep/params.config \
	run main.nf \
	-profile slurm -w "/fastscratch/agarwa/nf-tmp/work" -with-timeline \
	-resume $resumeDir 

