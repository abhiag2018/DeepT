#!/bin/bash
#SBATCH --mem=200G
#SBATCH --qos batch 
#SBATCH --time 1-00:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

###
NXF_OPTS='-Xms16g -Xmx64g'

if [ "$#" -lt 1 ]; then
    echo "at least 1 argument <dtype (val/all/train/test)> <resumeID>"
    echo "$#"
    exit 1
fi

species=$1; shift	
dtype=$1; shift
resumeID=$1; shift
# nextflow -C NN-feature-prep/nextflow.config run NN-feature-prep/main.nf -profile slurm -w "/fastscratch/agarwa/nf-tmp/work" -with-timeline -resume $resumeID
# nextflow -C DeepTact-input-preprocessing/nextflow-test.config run preprocessing/main.nf -profile slurm -w "/fastscratch/agarwa/nf-tmp/work"

basedir="/projects/li-lab/agarwa/CUBE/DeepTact/code/NN-feature-prep"

nextflow run $basedir/main.nf \
	-profile slurm \
	-w "/fastscratch/agarwa/nf-tmp/work" -with-timeline \
	--species $species \
	--dtype $dtype \
	-resume ${resumeID} "$@"
	# --splitParts 10
	## --splitByChr \
	## --dev

