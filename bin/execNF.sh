#!/bin/bash

# a/i. 
run_pr_enh.sh /projects/li-lab/agarwa/CUBE/DeepTact/code/configs_hg38_mm10/ hg

# a/ii. 
run_genome_seq.sh /projects/li-lab/agarwa/CUBE/DeepTact/code/configs_hg38_mm10/ hg

# b.
run_co_score.sh /projects/li-lab/agarwa/CUBE/DeepTact/code/configs_hg38_mm10/ hg 

# a/iii
run_pchic.sh /projects/li-lab/agarwa/CUBE/DeepTact/code/configs_hg38_mm10/ hg

# a/iv.
run_feature_gen.sh /projects/li-lab/agarwa/CUBE/DeepTact/code/configs_hg38_mm10/ hg --splitByChr

