#!/bin/bash
#SBATCH --mem=10G
#SBATCH --qos batch 
#SBATCH --time 1-00:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

###
NXF_OPTS='-Xms16g -Xmx64g'

resumeDir=$1
# nextflow -C NN-feature-prep/nextflow.config run NN-feature-prep/main.nf -profile slurm -w "/fastscratch/agarwa/nf-tmp/work" -with-timeline -resume $resumeDir
# nextflow -C DeepTact-input-preprocessing/nextflow-test.config run preprocessing/main.nf -profile slurm -w "/fastscratch/agarwa/nf-tmp/work"

# gzip -d storeDir/*.gz
nextflow -c NN-feature-prep/nextflow-test.config \
	-c  NN-feature-prep/configs/co-score-prep.config \
	-c  NN-feature-prep/configs/genome-seq-prep.config \
	-c  NN-feature-prep/configs/pchic-prep-test.config \
	-c  NN-feature-prep/configs/pr-enh-prep.config \
	run NN-feature-prep/main.v1.nf -profile slurm \
	--dev \
	-w "/fastscratch/agarwa/nf-tmp/work" -with-timeline \
	-resume $resumeDir
# gzip storeDir/*.npz 
