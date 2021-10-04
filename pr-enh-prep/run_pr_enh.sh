#!/bin/bash
#SBATCH --mem=10G
#SBATCH --qos batch 
#SBATCH --time 6:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

###

basedir="/projects/li-lab/agarwa/CUBE/DeepTact/code/pr-enh-prep"

species=$1; shift	
resumeID=$1; shift	
nextflow run $basedir/main.nf \
	-profile slurm \
	-w "/fastscratch/agarwa/nf-tmp/work" -with-timeline \
	--species $species \
	-resume $resumeID "$@"
# --dev \
