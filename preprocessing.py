"""
- specify intermediate and final preprocessing output file names and path
- do elementary pre-processing
"""
import colored_traceback.always
import os, glob, shutil, sys, itertools
import pandas as pd
import preptools as pt
import argparse, itertools
from parameters import baseDataDir, tmpBaseDir, codeTmpDir, bgWindow, promoter, enhancer, bamfilesInit, bamDir, intersectOptions, clearRun, reRun, hg19, hicTSV, gtf


baseBaseName = lambda f: f.split("/")[-1].split(".")[0] # extract basename from file path (no extension)
midBaseName = lambda f: f.split("/")[-1].split(".")[1:-1]  # extract middle names list from file path (no extension)


NameSys = lambda bam_file,name,append,ext: f"""{baseBaseName(bam_file) if bam_file else "."}/{name}{"."+".".join(append) if append else ""}.{ext}"""
# rel path and file name format for .bam file intersection with promoter/enhancer.
# input : 
#     bam_file: specifies the rel directory
#     name: main name of file (specified whether .bam file is intersected with background length or main .bed window)
#     append: append to specify split file chromosome or window index
#     ext: extension

nameWinBed = 'win' # main name for .bed file split into 1k/2k/.. windows. (Used later for intersection with .bam files for DNase input to DeepTact)
nameIntersectBG = 'intersect_BG' # main name for .bam and background .bed file intersection output
nameIntersectWin = 'intersect_win' # main name for .bam and window .bed file intersection output
nameDnaseOut = 'DNase' # main name for .npz output for DNase-Seq data
nameDNAout = 'DNA_Seq' # main name for .npz output for DNA-Sequence data

## promoter parameters
promoter['tmp-dir']= f"{tmpBaseDir}/tmp_promoter.1" # temporary data directory for promoter

promoter['bed-path']=  f"{promoter['tmp-dir']}/hg19_promoter_TSS.bed"  #main bed file containing element with the appropriate length 1k/2k/..
promoter['bg-path']= f"{promoter['tmp-dir']}/hg19_promoter_BG.bed" #bed file defining the background for element with the appropriate length 1MB..

promoter['intersect-tasklist'] = f"{promoter['tmp-dir']}/TaskList.npz"  # tasklist for pt.distribute_task(..) function to run 'bedtools intersect'
promoter['combine-bed-tasklist'] = f"{promoter['tmp-dir']}/CombineBedTasklist.npz"  # tasklist for pt.distribute_task(..) function to combine output of chromosome split 'bedtools intersect' output

promoter['winBed-path']=  lambda X:f"""{promoter['tmp-dir']}/{NameSys('',nameWinBed,[str(X)],'bed')}""" # .bed filepath after splitting promoter's main .bed file into windows for intersection with .bam
promoter['bg-intersect-out'] =  lambda bam_file: f"""{promoter['tmp-dir']}/{NameSys(bam_file,nameIntersectBG,midBaseName(bam_file),"bed")}""" # bedtools intersect output filepath for intersection of .bam file with element's background .bed file
promoter['intersect-out']=  lambda bam_file,X:f"""{promoter['tmp-dir']}/{NameSys(bam_file,nameIntersectWin,[str(X)]+midBaseName(bam_file),"bed")}""" # bedtools intersect output filepath for intersection of .bam file with element's window .bed file. X : window index
promoter['dnase-out']=  lambda bam_file:f"""{promoter['tmp-dir']}/{NameSys(bam_file,nameDnaseOut,[],"npz")}""" # DNase-Seq output .npz file

promoter['fa-out']=f"{promoter['tmp-dir']}/promoter.fa" # fasta file for promoters of length 1k/2k/..
promoter['dna-out'] = f"{promoter['tmp-dir']}/{nameDNAout}.npz" # DNA-Sequence output .npz file

## enhancer parameters
enhancer['tmp-dir']= f"{tmpBaseDir}/tmp_enhancer.1"

enhancer['bed-path']= f"{enhancer['tmp-dir']}/human_enhancers_window.bed"
enhancer['bg-path']= f"{enhancer['tmp-dir']}/hg19_enhancer_BG.bed"

enhancer['intersect-tasklist'] = f"{enhancer['tmp-dir']}/TaskList.npz"
enhancer['combine-bed-tasklist'] = f"{enhancer['tmp-dir']}/CombineBedTasklist.npz"

enhancer['winBed-path']=  lambda X:f"""{enhancer['tmp-dir']}/{NameSys('',nameWinBed,[str(X)],'bed')}"""
enhancer['bg-intersect-out'] =  lambda bam_file: f"""{enhancer['tmp-dir']}/{NameSys(bam_file,nameIntersectBG,midBaseName(bam_file),"bed")}"""
enhancer['intersect-out']=  lambda bam_file,X:f"""{enhancer['tmp-dir']}/{NameSys(bam_file,nameIntersectWin,[str(X)]+midBaseName(bam_file),"bed")}"""
enhancer['dnase-out']=  lambda bam_file:f"""{enhancer['tmp-dir']}/{NameSys(bam_file,nameDnaseOut,[],"npz")}"""

enhancer['fa-out']=f"{enhancer['tmp-dir']}/enhancer.fa"
enhancer['dna-out'] = f"{enhancer['tmp-dir']}/{nameDNAout}.npz"


## DNase Input

all_chrom = list(range(1,23))+list('XY')
chromAppend = lambda chrom:f"REF_chr{chrom}" # default for bamtools split
bamfiles = list(itertools.chain.from_iterable(glob.glob(f"""{os.path.dirname(bamf)}/*.{chromAppend('*')}.bam""") for bamf in bamfilesInit))


## Hi-C Input
HiC_Match = lambda cell:f"{tmpBaseDir}/pchicMatch_{cell}.pkl"
HiC_GroupMatch = lambda cell:f"{tmpBaseDir}/pchicGroupMatch_{cell}.pkl"
HiC_UniqueMatchPE = lambda cell:f"{tmpBaseDir}/pchicUniqueMatchPE_{cell}.pkl"
HiC_UniqueMatchEP = lambda cell:f"{tmpBaseDir}/pchicUniqueMatchEP_{cell}.pkl"
HiC_TrainingPos = lambda cell:f"{tmpBaseDir}/pchicTrainingPos_{cell}.csv"
HiC_Training = lambda cell:f"{tmpBaseDir}/pchicTraining_{cell}.csv"

HiCParts = [f"{tmpBaseDir}/pchic_{i}.csv" for i in range(10)]

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
 
    assert args.taskType in ["prep","prepPr","prepEnh","splitBam","prepHiC"]

    clean_run = clearRun

    if args.taskType=="splitBam":
        from parameters import bamfilesInit
        split_bam(args, bamfilesInit)

    if args.taskType=="prep" or args.taskType=="prepPr":
        pt.makedirs(promoter['tmp-dir'],exist_ok=not clean_run)
        for bam_file in bamfilesInit:
            pt.makedirs(os.path.dirname(promoter['bg-intersect-out'](bam_file)),exist_ok=not clean_run)
        elem_preprocessing(promoter, bgWindow, pt.process_promoter_bed)

    if args.taskType=="prep" or args.taskType=="prepEnh":
        pt.makedirs(enhancer['tmp-dir'],exist_ok=not clean_run)
        for bam_file in bamfilesInit:
            pt.makedirs(os.path.dirname(enhancer['bg-intersect-out'](bam_file)),exist_ok=not clean_run)
        elem_preprocessing(enhancer, bgWindow, pt.process_enhancer_bed)

    if args.taskType=="prep" or args.taskType=="prepHiC":
        pt.splitCSV(hicTSV, HiCParts, readArgs= {'delimiter':'\t'}, writeArgs= {'index':False,'sep':'\t'})



