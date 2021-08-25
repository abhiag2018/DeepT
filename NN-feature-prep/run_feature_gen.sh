#!/bin/bash
#SBATCH --mem=200G
#SBATCH --qos batch 
#SBATCH --time 3-00:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

###
NXF_OPTS='-Xms16g -Xmx64g'

if [ "$#" -lt 1 ]; then
    echo "at least 1 argument <dtype (val/all/train/test)> <resumeID>"
    echo "$#"
    exit 1
fi

dtype=$1; shift
resumeID=$1; shift
# nextflow -C NN-feature-prep/nextflow.config run NN-feature-prep/main.nf -profile slurm -w "/fastscratch/agarwa/nf-tmp/work" -with-timeline -resume $resumeID
# nextflow -C DeepTact-input-preprocessing/nextflow-test.config run preprocessing/main.nf -profile slurm -w "/fastscratch/agarwa/nf-tmp/work"

STORE_DIR="/projects/li-lab/agarwa/CUBE/DeepTact/code/NN-feature-prep"

# gzip -d storeDir/*.gz
nextflow -c $STORE_DIR/nextflow.config \
	-c $STORE_DIR/configs/co-score-prep.config \
	-c $STORE_DIR/configs/genome-seq-prep.config \
	-c $STORE_DIR/configs/pchic-prep.config \
	-c $STORE_DIR/configs/pr-enh-prep.config \
	run $STORE_DIR/main.nf -profile slurm \
	-w "/fastscratch/agarwa/nf-tmp/work" -with-timeline \
	--dtype $dtype \
	-resume $resumeID "$@" \
	## --dev

