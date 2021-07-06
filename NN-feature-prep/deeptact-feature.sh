#!/bin/bash
#SBATCH --mem=1G
#SBATCH --qos batch 
#SBATCH --time 3-00:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

###
NXF_OPTS='-Xms16g -Xmx64g'

resumeDir=$1; shift
# nextflow -C NN-feature-prep/nextflow.config run NN-feature-prep/main.nf -profile slurm -w "/fastscratch/agarwa/nf-tmp/work" -with-timeline -resume $resumeDir
# nextflow -C DeepTact-input-preprocessing/nextflow-test.config run preprocessing/main.nf -profile slurm -w "/fastscratch/agarwa/nf-tmp/work"

# gzip -d storeDir/*.gz
nextflow -c nextflow.config \
	-c  configs/co-score-prep.config \
	-c  configs/genome-seq-prep.config \
	-c  configs/pchic-prep.config \
	-c  configs/pr-enh-prep.config \
	run main.v1.nf -profile slurm \
	-w "/fastscratch/agarwa/nf-tmp/work" -with-timeline \
	-resume $resumeDir "$@" \
	## --dev
# gzip storeDir/*.npz 
