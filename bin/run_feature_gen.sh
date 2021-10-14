#!/bin/bash
#SBATCH --mem=200G
#SBATCH --qos batch 
#SBATCH --time 0-12:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

###
NXF_OPTS='-Xms16g -Xmx64g'

if [ "$#" -lt 1 ]; then
    echo "at least 1 argument <dtype (val/all/train/test)> <resumeID>"
    echo "$#"
    exit 1
fi

configDir=$1; shift
species=$1; shift	
dtype=$1; shift
# resumeID=$1; shift

basedir=`dirname $configDir`
ln -sf $configDir $basedir/configs

nextflow run $basedir/NN-feature-prep/main.nf \
	-profile slurm \
	-w "/fastscratch/agarwa/nf-tmp/work" -with-timeline \
	--save_dir `pwd` \
	--species $species \
	--dtype $dtype "$@"
	# -resume ${resumeID}
	# --splitParts 10
	## --splitByChr \
	## --dev

