"""
- specify intermediate and final preprocessing output file names and path
- do elementary pre-processing
"""
import colored_traceback.always
import os, glob, shutil, sys, itertools
import pandas as pd
import preptools as pt
import argparse, itertools
from parameters import baseDataDir, tmpBaseDir, codeTmpDir, bgWindow, promoter, enhancer, bamfilesInit, bamDir, clearRun, reRun, hg19, hicTSV, gtf, DnaseCells


## promoter parameters
promoter['bed-path']=  f"{tmpBaseDir}/hg19_promoter_TSS.bed"  #main bed file containing element with the appropriate length 1k/2k/..
promoter['bg-path']= f"{tmpBaseDir}/hg19_promoter_BG.bed" #bed file defining the background for element with the appropriate length 1MB..

promoter['fa-out']=f"{tmpBaseDir}/promoter.fa" # fasta file for promoters of length 1k/2k/..
promoter['dna-out'] = f"{tmpBaseDir}/promoterDNA.npz" # DNA-Sequence output .npz file


## enhancer parameters
enhancer['bed-path']= f"{tmpBaseDir}/human_enhancers_window.bed"
enhancer['bg-path']= f"{tmpBaseDir}/hg19_enhancer_BG.bed"

enhancer['fa-out']=f"{tmpBaseDir}/enhancer.fa"
enhancer['dna-out'] = f"{tmpBaseDir}/enhancerDNA.npz"


## DNase Input
dnaseTmpDir = lambda bamf: f"""{tmpBaseDir}/{os.path.basename(bamf).split(".")[0]}""" # bam processing output 

all_chrom = [str(i) for i in range(1,23)]+list('XY')
chromLen_GRCh37v13 = {'chr1':249250621, 'chr2':243199373, 'chr3':198022430, 'chr4':191154276, 'chr5':180915260, 'chr6':171115067, 'chr7':159138663, 'chr8':146364022, 'chr9':141213431, 'chr10':135534747, 'chr11':135006516, 'chr12':133851895, 'chr13':115169878, 'chr14':107349540, 'chr15':102531392, 'chr16':90354753, 'chr17':81195210, 'chr18':78077248, 'chr19':59128983, 'chr20':63025520, 'chr21':48129895, 'chr22':51304566, 'chrX':155270560, 'chrY':59373566}
chromAppend = lambda chrom:f"REF_chr{chrom}" # default for bamtools split
bamfiles = list(itertools.chain.from_iterable(glob.glob(f"""{os.path.dirname(bamf)}/*.{chromAppend('*')}.bam""") for bamf in bamfilesInit))


## Hi-C Input
tmpBaseDir_hic = f"{tmpBaseDir}/TrainingData"
HiC_Match = lambda cell:f"{tmpBaseDir_hic}/pchicMatch_{cell}.pkl"
HiC_GroupMatch = lambda cell:f"{tmpBaseDir_hic}/pchicGroupMatch_{cell}.pkl"
HiC_UniqueMatchPE = lambda cell:f"{tmpBaseDir_hic}/pchicUniqueMatchPE_{cell}.pkl"
HiC_UniqueMatchEP = lambda cell:f"{tmpBaseDir_hic}/pchicUniqueMatchEP_{cell}.pkl"
HiC_TrainingPos = lambda cell:f"{tmpBaseDir_hic}/pchicTrainingPos_{cell}.csv"
HiC_Training = lambda cell:f"{tmpBaseDir_hic}/pchicTraining_{cell}.csv"

HiCParts = [f"{tmpBaseDir_hic}/pchic_{i}.csv" for i in range(10)]

def elem_preprocessing(elem, bg_window, func):
    """ 
    element pre-processing function : from the downloaded files for the elements --promoters and enhancers-- creates 
        the background and main bed file with the appropriate window sizes.
        The output path is specified in the variables *element*['bed-path'] and *element*['bg-path'] wher *element* = [promoter, enhancer]
    """
    bed_path = elem['bed-path']
    headers = elem['headers']
    window = elem['window']

    bg_path = elem['bg-path']
    allfield_bed = elem['allfield-bed']

    print(f"generating main bed file: {os.path.basename(bed_path)}", end=".. ",flush=True)
    func( allfield_bed, bed_path, headers, window=window)
    print(".")

    print(f"generating background bed file: {os.path.basename(bg_path)}", end=".. ",flush=True)
    func( allfield_bed, bg_path, headers, window=bg_window)
    print(".")

    return 0


def split_bam(args, bamfiles):
    args.nTasks = min(args.nTasks,len(bamfiles))

    pt.distribute_task(task_list = bamfiles, 
        nTasks = int(args.nTasks), file_index=int(args.file_index), 
        func= pt.split_bam, num_tasks= len(bamfiles),
        dry_run=False)
    return 0


if __name__=="__main__":
    args = pt.process_inputArgs(input_parse=sys.argv[1:])
    # args = pt.process_inputArgs(input_parse=['--file_index',0,'--nTasks',52,'--taskType','pWin'])
 
    assert args.taskType in ['prepBam', 'prepPr' , 'prepEnh' , 'prepHiC', 'splitBam' ]

    clean_run = clearRun

    if args.taskType=="splitBam":
        from parameters import bamfilesInit
        split_bam(args, bamfilesInit)

    if args.taskType == "prepBam":
        for bam_file in bamfilesInit:
            pt.makedirs(dnaseTmpDir(bam_file),exist_ok=not clean_run)
    pt.makedirs(tmpBaseDir_hic,exist_ok=not clean_run)
    

    if args.taskType=="prepPr":
        elem_preprocessing(promoter, bgWindow, pt.process_promoter_bed)

    if args.taskType=="prepEnh":
        elem_preprocessing(enhancer, bgWindow, pt.process_enhancer_bed)

    if args.taskType=="prepHiC":
        pt.splitCSV(hicTSV, HiCParts, readArgs= {'delimiter':'\t'}, writeArgs= {'index':False,'sep':'\t'})



