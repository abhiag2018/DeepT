import colored_traceback.always
import os, glob, shutil, sys, itertools
import pandas as pd
import preptools as pt
import argparse, itertools
from parameters import baseDataDir, bgWindow, promoter, enhancer, bamfilesInit, bamDir, intersectOptions, clearRun, reRun, hg19


def makedirs(dirpath,exist_ok = False):
    try:
        os.makedirs(dirpath,exist_ok=exist_ok)
    except FileExistsError:
        shutil.rmtree(dirpath)
        os.makedirs(dirpath,exist_ok=exist_ok)
    return 0

baseBaseName = lambda f: f.split("/")[-1].split(".")[0]
midBaseName = lambda f: f.split("/")[-1].split(".")[1:-1]

NameSys = lambda bam_file,name,append,ext: f"""{baseBaseName(bam_file) if bam_file else "."}/{name}{"."+".".join(append) if append else ""}.{ext}"""


nameWinBed = 'win'
nameIntersectBG = 'intersect_BG'
nameIntersectWin = 'intersect_win'
nameDnaseOut = 'DNase'

tmpBaseDir = "/fastscratch/agarwa/DeepTact_tmp"

## promoter parameters
promoter['bed-path']=  f"{baseDataDir}/hg19_promoter_TSS.bed"
promoter['bg-path']= f"{baseDataDir}/hg19_promoter_BG.bed"
promoter['tmp-dir']= f"{tmpBaseDir}/tmp_promoter"

promoter['intersect-tasklist'] = f"{promoter['tmp-dir']}/prTaskList.npz"

promoter['winBed-path']=  lambda X:f"""{promoter['tmp-dir']}/{NameSys('',nameWinBed,[str(X)],'bed')}"""
promoter['bg-intersect-out'] =  lambda bam_file: f"""{promoter['tmp-dir']}/{NameSys(bam_file,nameIntersectBG,midBaseName(bam_file),"bed")}"""
promoter['intersect-out']=  lambda bam_file,X:f"""{promoter['tmp-dir']}/{NameSys(bam_file,nameIntersectWin,[str(X)]+midBaseName(bam_file),"bed")}"""
promoter['dnase-out']=  lambda bam_file:f"""{promoter['tmp-dir']}/{NameSys(bam_file,nameDnaseOut,[],"npz")}"""

promoter['fa-out']=f"{baseDataDir}/promoter.fa"

## enhancer parameters
enhancer['bed-path']= f"{baseDataDir}/enhancers_fantom5/human_enhancers_window.bed"
enhancer['bg-path']= f"{baseDataDir}/hg19_enhancer_BG.bed"
enhancer['tmp-dir']= f"{tmpBaseDir}/tmp_enhancer"

enhancer['intersect-tasklist'] = f"{enhancer['tmp-dir']}/enhTaskList.npz"

enhancer['winBed-path']=  lambda X:f"""{enhancer['tmp-dir']}/{NameSys('',nameWinBed,[str(X)],'bed')}"""
enhancer['bg-intersect-out'] =  lambda bam_file: f"""{enhancer['tmp-dir']}/{NameSys(bam_file,nameIntersectBG,midBaseName(bam_file),"bed")}"""
enhancer['intersect-out']=  lambda bam_file,X:f"""{enhancer['tmp-dir']}/{NameSys(bam_file,nameIntersectWin,[str(X)]+midBaseName(bam_file),"bed")}"""
enhancer['dnase-out']=  lambda bam_file:f"""{enhancer['tmp-dir']}/{NameSys(bam_file,nameDnaseOut,[],"npz")}"""

enhancer['fa-out']=f"{baseDataDir}/enhancer.fa"


## DNase Input


all_chrom = list(range(1,23))+list('XY')
chromAppend = lambda chrom:f"REF_chr{chrom}" # default for bamtools split
bamfiles = list(itertools.chain.from_iterable(glob.glob(f"""{os.path.dirname(bamf)}/*.{chromAppend('*')}.bam""") for bamf in bamfilesInit))

def elem_preprocessing(elem, bg_window, func):
    """ element pre-processing"""
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
 
    assert args.taskType in ["prep","prepPr","prepEnh","splitBam"]

    clean_run = clearRun

    if args.taskType=="splitBam":
        from parameters import bamfilesInit
        split_bam(args, bamfilesInit)

    if args.taskType=="prep" or args.taskType=="prepPr":
        makedirs(promoter['tmp-dir'],exist_ok=not clean_run)
        for bam_file in bamfilesInit:
            makedirs(os.path.dirname(promoter['bg-intersect-out'](bam_file)),exist_ok=not clean_run)
        # elem_preprocessing(promoter, bgWindow, pt.process_promoter_bed)

    if args.taskType=="prep" or args.taskType=="prepEnh":
        makedirs(enhancer['tmp-dir'],exist_ok=not clean_run)
        for bam_file in bamfilesInit:
            makedirs(os.path.dirname(enhancer['bg-intersect-out'](bam_file)),exist_ok=not clean_run)
        # elem_preprocessing(enhancer, bgWindow, pt.process_enhancer_bed)




