import colored_traceback.always
import os, glob, shutil, sys
import pandas as pd
import preptools as pt
import argparse, itertools
from parameters import baseDataDir, bgWindow, promoter, enhancer, bamDir, intersectOptions, clearRun, reRun, hg19


def makedirs(dirpath,exist_ok = False):
    try:
        os.makedirs(dirpath,exist_ok=exist_ok)
    except FileExistsError:
        shutil.rmtree(dirpath)
        os.makedirs(dirpath,exist_ok=exist_ok)
    return 0

basename = lambda f: f.split("/")[-1].split(".")[0]+ext

def tmp_bam_dir(tmp_basedir,bam_file):
    bamf_name = pt.replace_suffix(os.path.basename(bam_file),ext='')
    tmpBAMdir = f"{tmp_basedir}/{bamf_name}"
    os.makedirs(tmpBAMdir,exist_ok=True)
    return tmpBAMdir

## promoter parameters
promoter['bed-path']=  "hg19_promoter_TSS.bed"
promoter['bg-path']= "hg19_promoter_BG.bed"
promoter['tmp_dir']= "tmp_promoter"

promoter['winBed-path']=  lambda X:f"{promoter['tmp_dir']}/pr_win_{X}.bed"
promoter['bg-intersect-out']=  lambda bam_file: f"{tmp_bam_dir(promoter['tmp_dir'],bam_file)}/intersect_pr_BG.bed"
promoter['intersect-out']=  lambda bam_file,X:f"{tmp_bam_dir(promoter['tmp_dir'],bam_file)}/intersect_pr_{X}.bed"
promoter['dnase-out']=  lambda bam_file:f"{tmp_bam_dir(promoter['tmp_dir'],bam_file)}/DNase_pr.npz"

promoter['fa-out']="promoter.fa"

## enhancer parameters
enhancer['bed-path']= "enhancers_fantom5/human_enhancers_window.bed"
enhancer['bg-path']= "hg19_enhancer_BG.bed"
enhancer['tmp_dir']= "tmp_enhancer"

enhancer['winBed-path']= lambda X:f"{enhancer['tmp_dir']}/enh_win_{X}.bed"
enhancer['bg-intersect-out']= lambda bam_file: f"{tmp_bam_dir(enhancer['tmp_dir'],bam_file)}/intersect_enh_BG.bed"
enhancer['intersect-out']= lambda bam_file,X:f"{tmp_bam_dir(enhancer['tmp_dir'],bam_file)}/intersect_enh_{X}.bed"
enhancer['dnase-out']= lambda bam_file:f"{tmp_bam_dir(enhancer['tmp_dir'],bam_file)}/DNase_enh.npz"

enhancer['fa-out']="enhancer.fa"


## DNase Input

bamfiles = glob.glob(f"{bamDir}/**/*.REF_chr*.bam",recursive=True)

def elem_preprocessing(basedir, elem, bg_window, func):
    """ element pre-processing"""
    bed_path = f"{basedir}/{elem['bed-path']}"
    headers = elem['headers']
    window = elem['window']

    bg_path = f"{basedir}/{elem['bg-path']}"
    allfield_bed = f"{basedir}/{elem['allfield-bed']}"

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
        makedirs(f"{baseDataDir}/{promoter['tmp_dir']}",exist_ok=not clean_run)
        for bam_file in bamfiles:
            makedirs(os.path.dirname(f"{baseDataDir}/{promoter['bg-intersect-out'](bam_file)}"),exist_ok=not clean_run)
        elem_preprocessing(baseDataDir, promoter, bgWindow, pt.process_promoter_bed)

    if args.taskType=="prep" or args.taskType=="prepEnh":
        makedirs(f"{baseDataDir}/{enhancer['tmp_dir']}",exist_ok=not clean_run)
        for bam_file in bamfiles:
            makedirs(os.path.dirname(f"{baseDataDir}/{enhancer['bg-intersect-out'](bam_file)}"),exist_ok=not clean_run)
        elem_preprocessing(baseDataDir, enhancer, bgWindow, pt.process_enhancer_bed)




