#!/bin/bash
#SBATCH --mem=10G
#SBATCH --qos batch 
#SBATCH --time 6:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

###
NXF_OPTS='-Xms16g -Xmx64g'



SCRIPT_DIR="/projects/li-lab/agarwa/CUBE/DeepTact/code/co-score-prep"

# resumeDir=$1; shift
# bamInput=$1; shift
nextflow  -c $SCRIPT_DIR/nextflow.config \
	-c  $SCRIPT_DIR/params.config \
	-c $SCRIPT_DIR/../pr-enh-prep/params.config \
	run $SCRIPT_DIR/main.nf \
	-profile slurm -w "/fastscratch/agarwa/nf-tmp/work" -with-timeline \
	-resume "$@"
	# --bamInput $bamInput

sbatch co-score-prep --bamInput Mon.csv
