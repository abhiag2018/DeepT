"""
- specify input files needed for pre-processing the input needed for DeepTact.
- specify the pre-processing, input feature parameters 
- only called through preprocessing.py
"""
import colored_traceback.always
import os
import re, glob
from collections import OrderedDict

import pandas as pd

import preptools as pt


baseDataDir = "/projects/li-lab/agarwa/CUBE/DeepTact/dataset" #dataset base directory
AUG_LEN = 1000
AUG_STEP = 50

# directory to store most of the intermediate processed files
# tmpBaseDir = "/fastscratch/agarwa/DeepTact_tmp.1"
tmpBaseDir = '/projects/li-lab/agarwa/CUBE/DeepTact/dataset/DeepTact_tmp.2'
codeTmpDir = f"{tmpBaseDir}/tmp_data" #directory for short lived (temporary) files 
bgWindow = int(1e6) # background length for regulatory elements = 1MB (both promoters and enhancers)

hg19 = f"{baseDataDir}/hg19.fa" #DNA-Sequence input : fasta file

## promoter parameters
promoter = {'headers': ['chrom', 'txStart', 'txEnd', 'name', 'score', 'strand',
               'cdsStart', 'cdsEnd', '--', 'exonCount', 'exonStarts',
               'exonEnds'], # headers for all_field.bed transcriptID file downloaded from ENCODE. '--' is a column that does not exist in all_field.bed
    'window': 1000,  # length of DNA sequence & DNase-Seq data
    'allfield-bed': f"{baseDataDir}/hg19_promoter_allFields.bed" # all_field.bed transcriptID file downloaded from ENCODE
}
promoter['ext-window'] = promoter['window'] + AUG_LEN
promoter['aug_step'] = AUG_STEP

    
## enhancer parameters
enhancer = {'headers': ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
           'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes',
           'blockStarts'], # headers for .bed file downloaded from FANTOM5
    'window': 2000, # length of DNA sequence & DNase-Seq data
    'allfield-bed': f"{baseDataDir}/enhancers_fantom5/human_permissive_enhancers_phase_1_and_2.bed" # bed file downloaded from FANTOM5
    # 'allfield-bed': f"{baseDataDir}/enhancers_EnhancerAtlas2/human_permissive_enhancers_phase_1_and_2.bed" # bed file downloaded from FANTOM5
    # 'allfield-bed': f"{baseDataDir}/enhancers_fantom5/F5.hg38.enhancers.bed" # bed file downloaded from FANTOM5
    }
enhancer['ext-window'] = enhancer['window'] + AUG_LEN
enhancer['aug_step'] = AUG_STEP

## DNase input parameters
LIMITBAM = None #limit bamfiles pre-processing to make development set; set to None to process everything
# bamDir = f"{baseDataDir}/Dnase-Seq/cellTypes" # .bam files for DNase-Seq input
# bamfilesInit = [fn for fn in glob.glob(f"{bamDir}/**/*.bam",recursive=True) if re.match("[^.]*\.bam$",fn)][:LIMITBAM] # list the .bam files to process inside the main directory
bamDir = f"{baseDataDir}/ATAC-Seq/nCD4_nature2018" # .bam files for DNase-Seq input
# bamfilesInit = [fn for fn in glob.glob(f"{bamDir}/*.bam",recursive=True) if re.match("[^.]*\.bam$",fn)][:LIMITBAM] # list the .bam files to process inside the main directory
bamfilesInit = [fn for fn in glob.glob(f"{bamDir}/*.bam",recursive=True)][:LIMITBAM] # list the .bam files to process inside the main directory

DnaseCells = OrderedDict()  #cell type info for the .bam files
DnaseCells['nCD4'] = [['naiveCD4_R1.mLb.clN.sorted.bam'], ['naiveCD4_R2.mLb.clN.sorted.bam'], ['naiveCD4_R3.mLb.clN.sorted.bam']]
# DnaseCells['nCD4'] = [['ENCFF138MTV.bam'], ['ENCFF947LCB.bam']]
# DnaseCells['tCD4'] = [['ENCFF145YPS.bam'] ,['ENCFF938VVP.bam'] ,['ENCFF308TOC.bam']]
# DnaseCells['tCD8'] = [['ENCFF054ZTY.bam'], ['ENCFF736QAD.bam'],['ENCFF790RMQ.bam']]
# DnaseCells['tB'] = [['ENCFF203BEH.bam'] ,['ENCFF469VSO.bam'] ,['ENCFF504HES.bam', 'ENCFF552WWJ.bam']]
# DnaseCells['Mon'] = [['ENCFF295OEK.bam', 'ENCFF780FKA.bam'],['ENCFF175LYA.bam'] ,['ENCFF227VKS.bam'] ,['ENCFF418HYD.bam']] 
# DnaseCells['FoeT'] = [['ENCFF840GSK.bam'] ,['ENCFF487IUY.bam'] ,['ENCFF315TUQ.bam'] ,['ENCFF127FDA.bam'] ,['ENCFF719KCS.bam'] ,['ENCFF149PAO.bam'] ,['ENCFF224XWI.bam'] ,['ENCFF347QEH.bam'] ,['ENCFF397UIC.bam']]


## PCHiC Training Data

hicTSV = f"{baseDataDir}/Javierre_ref_18/DATA_S1/PCHiC_peak_matrix_cutoff5.tsv" # PCHi-C input file

# gtf = f"{baseDataDir}/gencode.v19.annotation.gtf"
gtf = f"{baseDataDir}/gencode.v3c.annotation.GRCh37.gtf" #ref release. 09/2009 #gene annotation file to match Transcript IDs to gene symbols

## Running Parameters

clearRun = False # setting this to True deletes everything in preprocessing.py
reRun = False #used to check whether to re-run computation or just use the existing intermediate files when available
























