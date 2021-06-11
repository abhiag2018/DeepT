#!/usr/bin/env bash
gzip -df $enhancer_hic_CO_gz $promoter_hic_CO_gz $enhancer_hic_DNAseq_gz $promoter_hic_DNAseq_gz

python ${projectDir}/templates/separate_data.py $hic_aug $enhancer_hic_DNAseq $promoter_hic_DNAseq $enhancer_hic_CO $promoter_hic_CO $params.sepdata_split


gzip input_feature_part_*.h5
rm $enhancer_hic_CO $promoter_hic_CO $enhancer_hic_DNAseq $promoter_hic_DNAseq