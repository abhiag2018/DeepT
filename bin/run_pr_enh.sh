#!/bin/bash
#SBATCH --mem=10G
#SBATCH --qos batch 
#SBATCH --time 6:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

###

#basedir="/projects/li-lab/agarwa/CUBE/DeepTact/code"

configDir=$1; shift
species=$1; shift	
resumeID=$1; shift	

basedir=`dirname $configDir` 
ln -sf $configDir $basedir/configs

nextflow run $basedir/pr-enh-prep/main.nf \
	-profile slurm \
	-w "/fastscratch/agarwa/nf-tmp/work" -with-timeline \
	--species $species \
	-resume $resumeID "$@"
# --dev \
