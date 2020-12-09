# Intro

Preprocess/generate input features and training labels for DeepTact. 

To run:
	sbatch run_deeptact_prep.sh <resume-session-id>
The session-id can be found by running command : *nextflow log*


# Input Params and Files

The input file path is specified in the nextflow config file parameters. The input .bam files for DNase/ATAC-seq input are specified separately in a .csv specified in $params.coScore_data


# Preprocessing Steps

	1. Preprocess regulatory element lists (PREPROCESS_PROMOTER, PREPROCESS_ENHANCER)
	2 Construct chromatin openness score (CO score) profiles for input reg elements (from step 1):
		1. Generate .bam.bai index files for ATAC/DNase input bam files (GEN_BAM_INDEX)
		2. Construct CO score from .bam files  (CHROM_OPENN_SCORE)
		3. Normalize CO score and output CO score for reg elements (CHROM_OPENN_SCORE_PROFILE_ENHANCER, CHROM_OPENN_SCORE_PROFILE_PROMOTER)
	3. Generate genomic sequence for reg element lists from fasta files (GENOMIC_SEQUENCE_PROMOTER, GENOMIC_SEQUENCE_ENHANCER)
	4. Training Label Preparation
		1. split hic input for faster preprocessing into *$params.hic_split_process* parts (SPLIT_HIC)
		2. match hic input pairs to reg elements from step 1 (MATCH_HIC_TO_REG_ELEMENTS)
		3. combine matched files (COMBINE_MATCHED_HIC)
		4. combine promoters which correspond to the same gene i.e. in the promoter list there are a lot of promoters that are essentially the same; combine them (COMBINE_PROMOTERS)
		5. generate interactions for negative training data (GEN_NEG_LABEL)
		6. augment training data by a factor of *$params.hic_augment_factor* and output training .csv file (GEN_AUGMENTED_LABEL) 
	5. split hic augmented file into *$params.hic_split_combine* parts (SPLIT_HIC_AUG)
	6. Combine steps 1,2,3, and 4 ie. generate features for reg elemetns from step 5:
		1. combine CO score profile data :
			i COMBINE_PCHIC_CO_SCORE
			ii COMBINE_PCHIC_CO_SCORE_ENH, COMBINE_PCHIC_CO_SCORE_PR
			iii COMBINE_CO_SCORE_REPS
		2. combine DNA-seq data :
			i. COMBINE_PCHIC_DNA_SEQ
			ii. COMBINE_PCHIC_OUT_ENHANCER, COMBINE_PCHIC_DNA_SEQ_PROMOTER
		3. combine all the features into one file and break up by data points.
			i. SAVE_ENH_DNA_SEQ, SAVE_PR_DNA_SEQ, SAVE_ENH_CO_SCORE, SAVE_PR_CO_SCORE
			ii. SEPARATE_DATA


