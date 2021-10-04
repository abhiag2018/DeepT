#!/bin/bash
#SBATCH --mem=10G
#SBATCH --qos batch 
#SBATCH --time 6:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

###
NXF_OPTS='-Xms16g -Xmx64g'



basedir="/projects/li-lab/agarwa/CUBE/DeepTact/code/co-score-prep"

species=$1; shift	
bamInput=$1; shift
resumeID=$1; shift	
nextflow run $basedir/main.nf \
	-profile slurm -w "/fastscratch/agarwa/nf-tmp/work" -with-timeline \
	--species $species \
	--bamInput $bamInput \
	-resume $resumeID "$@"

#sbatch run_co_score.sh [species]hg [bamInput]tB.csv
