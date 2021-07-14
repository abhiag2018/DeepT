#!/bin/bash
#SBATCH --mem=200G
#SBATCH --qos batch 
#SBATCH --time 3-00:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

###
NXF_OPTS='-Xms16g -Xmx64g'

resumeDir=$1; shift
# nextflow -C NN-feature-prep/nextflow.config run NN-feature-prep/main.nf -profile slurm -w "/fastscratch/agarwa/nf-tmp/work" -with-timeline -resume $resumeDir
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
	-resume $resumeDir "$@" \
	## --dev
# gzip storeDir/*.npz 


# --pchic_data pchic_data-val.csv --dataType data_val
# --pchic_data pchic_data-test.csv --dataType data_test
# --pchic_data pchic_data-train.csv --dataType data_train
