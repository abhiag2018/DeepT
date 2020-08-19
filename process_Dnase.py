import colored_traceback.always
import glob
import sys, os, shutil, glob, re
import argparse, itertools
import pandas as pd
import numpy as np
import math
import preptools as pt
import functools



LIMIT = None

# generate element tasklist
def intersect_tasklist(elem,bamfiles, tasktype='all', recompute=True):
    numTasks = 0
    output_path = elem['intersect-tasklist']
    site_window = elem['siteWindow']
    num_windows = elem['window']-site_window+1

    assert tasktype in ['all','bg','win']
    bg_intersect_bam = []
    if tasktype == 'all' or tasktype=='bg':
        bg_intersect_bam = [("intersect-bg",-1,bam_file) for bam_file in bamfiles if recompute or not os.path.exists(elem['bg-intersect-out'](bam_file))]
        numTasks += len(bg_intersect_bam)

    win_intersect_bam = []
    if tasktype == 'all' or tasktype=='win':
        win_intersect_bam = [("intersect-win",i,bam_file) for bam_file in bamfiles for i in range(num_windows) if recompute or not os.path.exists(elem['intersect-out'](bam_file, i))]
        numTasks += len(win_intersect_bam)

    tasklist = bg_intersect_bam+win_intersect_bam

    numTasks = LIMIT if LIMIT else numTasks
    tasklist = tasklist[:numTasks]
    np.savez(output_path,tasklist=tasklist,numTasks=numTasks)
    return tasklist,numTasks


def generate_winBed(args, reRun, elem):
    bed_path = elem['bed-path']
    headers = elem['headers']
    winBedPath = lambda idx: elem['winBed-path'](idx)
    window =  elem['window']
    site_window = elem['siteWindow']

    def preprocess_elem(winIdx):
        outputf = winBedPath(winIdx)
        if reRun or not os.path.exists(outputf):
            pt.generate_elem_window_bed(bed_path, headers, winIdx, site_window, outputf)
        return 0

    numTasks = window-site_window+1
    numTasks = LIMIT if LIMIT else numTasks

    args.nTasks = min(args.nTasks,numTasks)

    pt.distribute_task(task_list = range(0,numTasks), 
        nTasks = int(args.nTasks), file_index=int(args.file_index), 
        func=preprocess_elem, num_tasks=numTasks,
        dry_run=False)
    return 0

def intersect_elem(args, elem, bamfiles, intersectOptions, tasktype='all', recompute=True):
    data = np.load(elem['intersect-tasklist'],allow_pickle=True)
    elem_tasks = data['tasklist']
    elem_numTasks = data['numTasks']

    def process_elem(task):
        task_type, winIdx, bam_file = task
        assert task_type in ["intersect-win","intersect-bg"]
        if task_type == "intersect-win":
            inputBed = elem['winBed-path'](winIdx)
            outBed = elem['intersect-out'](bam_file, winIdx)
        else:
            inputBed = elem['bg-path']
            outBed = elem['bg-intersect-out'](bam_file)
        pt.process_peaks(a=inputBed,b=bam_file, out_path = outBed,options = intersectOptions, recompute=recompute)
        # else:
            # print("not executing..")
        return 0

    args.nTasks = min(args.nTasks,elem_numTasks)

    output =  pt.distribute_task(task_list = elem_tasks, 
        nTasks = int(args.nTasks), file_index=int(args.file_index), 
        func= process_elem, num_tasks=elem_numTasks,
        dry_run=False)
    return output


def combineBed_tasks(elem, bgNameSys, winNameSys, all_chrom, recompute=True):
    namesys = lambda bamDir, index, appendList:winNameSys(bamDir, index, appendList) if index is not None else bgNameSys(bamDir, appendList)

    elemBaseDir = elem['tmp-dir']
    numWindows = elem['window']-elem['siteWindow']+1

    bamDirList = [os.path.dirname(f) for f in glob.glob(f"{elemBaseDir}/*" + os.path.sep)]
    tasklist = [(bamDir,None) for bamDir in bamDirList  if recompute or not os.path.exists(namesys(bamDir,None,[])) ]
    tasklist += [(bamDir,n) for n in range(numWindows) for bamDir in bamDirList if recompute or not os.path.exists(namesys(bamDir,n,[])) ]

    #check all chromosomes present
    for bamDir, index in tasklist:
        chromExist = [os.path.exists(namesys(bamDir, index, [c])) for c in all_chrom]
        if not functools.reduce(lambda x, y: x and y, chromExist):
            print(functools.reduce(lambda x, y: x and y, chromExist),chromExist)
            raise NameError(f"all chromosomes not found for {bamDir}; index:{index}")

    np.savez(elem['combine-bed-tasklist'],tasklist=tasklist)
    tasklist = tasklist[:LIMIT] if LIMIT else tasklist

    return tasklist

def combine_bed_output(args, elem, bgNameSys, winNameSys, all_chrom, recompute=True):
    data = np.load(elem['combine-bed-tasklist'],allow_pickle=True)
    tasklist = data['tasklist']

    namesys = lambda bamDir, index, appendList:winNameSys(bamDir, index, appendList) if index is not None  else bgNameSys(bamDir, appendList)

    headers = elem['headers']

    def combine_files(bamDir,index, recompute=False):
        namefunc = lambda X: namesys(bamDir, index, X)
        combine_input = [namefunc([c]) for c in all_chrom]
        outputf = namefunc([])
        if not recompute and os.path.exists(outputf):
            return outputf
        pt.combine_intersectBed(combine_input, headers, outputf, remove=False)
        return outputf


    args.nTasks = min(args.nTasks,len(tasklist))

    output =  pt.distribute_task(task_list = tasklist, 
        nTasks = int(args.nTasks), file_index=int(args.file_index), 
        func= lambda x:combine_files(*x, recompute=recompute), num_tasks=len(tasklist),
        dry_run=False)
    return output


def elem_postprocessing(args, reRun, elem, bg_window, bamfiles):
    bed_path = elem['bed-path']
    headers = elem['headers']
    window = elem['window']
    site_window = elem['siteWindow']

    bg_inter = elem['bg-intersect-out']
    out_p = elem['dnase-out']

    intersect_winfiles = elem['intersect-out']

    elem_bed = pd.read_csv(bed_path,delimiter="\t",names=headers)
    out_dim = (elem_bed.shape[0],window-site_window+1)

    def post_processing(bam_file):
        outputf = out_p(bam_file)
        if reRun or not os.path.exists(outputf):
            pt.post_process(bg_inter(bam_file), bg_window, headers, lambda X:intersect_winfiles(bam_file,X), out_dim, outputf)
        return 0

    args.nTasks = min(args.nTasks,len(bamfiles))

    pt.distribute_task(task_list = bamfiles, 
        nTasks = int(args.nTasks), file_index=int(args.file_index), 
        func= post_processing, num_tasks= len(bamfiles),
        dry_run=False)
    return 0



if __name__=="__main__":
    """
    outputs openness score computed from DNase-Seq data for elements.
    The data is sorted according to ['chrom','name'] in pt.post_process
    """
    import preprocessing as prep

    args = pt.process_inputArgs(input_parse=sys.argv[1:])
    # args = pt.process_inputArgs(input_parse=['--file_index','0','--nTasks','52','--taskType','pWin'])


    _taskTypes = ["Win",'TaskList',"Intersect",'cTaskList', "Combine", "PostProcess"]
    taskTypes = [t+_t for t in ['','p','e'] for _t in _taskTypes]
    assert args.taskType in taskTypes

    p = prep.promoter
    e = prep.enhancer
    bgwin = prep.bgWindow

    bamfsInit = prep.bamfilesInit
    bamfs = prep.bamfiles
    interesctOpt = prep.intersectOptions

    all_chrom = [prep.chromAppend(c) for c in prep.all_chrom ]
    bgNameSys = lambda bamDir, appendList:f"{bamDir}/{prep.NameSys('',prep.nameIntersectBG, appendList,'bed')}"
    winNameSys = lambda bamDir, index, appendList:f"{bamDir}/{prep.NameSys('',prep.nameIntersectWin, [str(index)]+appendList,'bed')}"

    reRun = prep.reRun

    if args.taskType=="pWin" or args.taskType=="Win":
        generate_winBed(args, reRun, p)
    if args.taskType=="eWin" or args.taskType=="Win":
        generate_winBed(args, reRun, e)

    if args.taskType=="TaskList" or  args.taskType=="pTaskList":
        intersect_tasklist(p, bamfs, tasktype='all', recompute=reRun)
    if args.taskType=="TaskList" or  args.taskType=="eTaskList":
        intersect_tasklist(e, bamfs, tasktype='all', recompute=reRun)

    if args.taskType=="pIntersect" or args.taskType=="Intersect":
        intersect_elem(args, p, bamfs, interesctOpt, tasktype='all', recompute=reRun)
    if args.taskType=="eIntersect" or args.taskType=="Intersect":
        intersect_elem(args, e, bamfs, interesctOpt, tasktype='all', recompute=reRun)

    if args.taskType=="cTaskList" or  args.taskType=="pcTaskList":
        combineBed_tasks(p, bgNameSys, winNameSys, all_chrom, recompute=reRun)
    if args.taskType=="cTaskList" or  args.taskType=="ecTaskList":
        combineBed_tasks(e, bgNameSys, winNameSys, all_chrom, recompute=reRun)

    if args.taskType=="pCombine" or args.taskType=="Combine":
        combine_bed_output(args, p, bgNameSys, winNameSys, all_chrom, recompute=reRun)
    if args.taskType=="eCombine" or args.taskType=="Combine":
        combine_bed_output(args, e, bgNameSys, winNameSys, all_chrom, recompute=reRun)

    if args.taskType=="pPostProcess" or args.taskType=="PostProcess":
        elem_postprocessing(args, reRun,  p, bgwin, bamfsInit)
    if args.taskType=="ePostProcess" or args.taskType=="PostProcess":
        elem_postprocessing(args, reRun,  e, bgwin, bamfsInit)
# python process_Dnase.py --nTasks=600  --file_index=$SLURM_ARRAY_TASK_ID

# sbatch -J prWinBed -q batch --ntasks=1 --cpus-per-task=1 --mem=200G --time 70:00:00 -o slurm_dnase_pr.out --array=0-600 run_script.sh
# sbatch -J prInters -q batch --ntasks=1 --cpus-per-task=1 --mem=200G --time 70:00:00 -o slurm_dnase_pr.out --array=0-600 run_script.sh

# sbatch -J enhWinBed -q batch --ntasks=1 --cpus-per-task=1 --mem=200G --time 70:00:00 -o slurm_dnase_enh.out --array=0-600 run_script.sh
# sbatch -J enhInters -q batch --ntasks=1 --cpus-per-task=1 --mem=200G --time 70:00:00 -o slurm_dnase_enh.out --array=0-600 run_script.sh
















