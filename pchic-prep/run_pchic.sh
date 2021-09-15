#!/bin/bash
#SBATCH --mem=10G
#SBATCH --qos batch 
#SBATCH --time 3-00:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

###
NXF_OPTS='-Xms16g -Xmx64g'

basedir="/projects/li-lab/agarwa/CUBE/DeepTact/code/pchic-prep"
resumeDir=$1
nextflow -c $basedir/nextflow.config \
	-c params.config \
	-c $basedir/../pr-enh-prep/params.config \
	-c $basedir/../genome-seq-prep/params.config \
	run $basedir/main.nf -profile slurm \
	-w "/fastscratch/agarwa/nf-tmp/work" -with-timeline -resume ${resumeDir} "$@"
	# --species mm 
	# --dev
