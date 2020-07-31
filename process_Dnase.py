import colored_traceback.always
import glob
import sys, os, shutil, glob, re
import argparse, itertools
import pandas as pd
import math
import preptools as pt

LIMIT = None

# generate element tasklist
def intersect_tasklist(elem,bamfiles, tasktype='all'):
    numTasks = 0
    site_window = elem['siteWindow']
    num_windows = elem['window']-site_window+1

    assert tasktype in ['all','bg','win']
    bg_intersect_bam = []
    if tasktype == 'all' or tasktype=='bg':
        bg_intersect_bam =  [("intersect-bg",-1,bam_file) for bam_file in bamfiles]
        numTasks += len(bamfiles)

    win_intersect_bam = []
    if tasktype == 'all' or tasktype=='win':
        for bam_file in bamfiles:
            win_intersect_bam += [("intersect-win",i,bam_file) for i in range(num_windows)]
        numTasks += num_windows*len(bamfiles)

    tasklist = bg_intersect_bam+win_intersect_bam

    numTasks = LIMIT if LIMIT else numTasks
    tasklist = tasklist[:numTasks]
    return tasklist,numTasks


def generate_winBed(args, reRun, baseDir, elem):
    bed_path = f"{baseDir}/{elem['bed-path']}"
    headers = elem['headers']
    winBedPath = lambda idx: f"{baseDir}/{elem['winBed-path'](idx)}"
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

def intersect_elem(args, reRun, baseDir, elem, bamfiles, intersectOptions, tasktype='all'):
    elem_tasks,elem_numTasks = intersect_tasklist(elem, bamfiles, tasktype='all')

    def process_elem(task):
        task_type, winIdx, bam_file = task
        assert task_type in ["intersect-win","intersect-bg"]
        if task_type == "intersect-win":
            inputBed = f"{baseDir}/{elem['winBed-path'](winIdx)}"
            outBed = f"{baseDir}/{elem['intersect-out'](bam_file, winIdx)}"
        else:
            inputBed = f"{baseDir}/{elem['bg-path']}"
            outBed = f"{baseDir}/{elem['bg-intersect-out'](bam_file)}"
        if reRun or not os.path.exists(outBed):
            pt.process_peaks(a=inputBed,b=bam_file, out_path = outBed,options = intersectOptions)
        # else:
            # print("not executing..")
        return 0

    args.nTasks = min(args.nTasks,elem_numTasks)

    output =  pt.distribute_task(task_list = elem_tasks, 
        nTasks = int(args.nTasks), file_index=int(args.file_index), 
        func= process_elem, num_tasks=elem_numTasks,
        dry_run=False)
    return output


def elem_postprocessing(args, reRun, baseDir, elem, bg_window, bamfiles):
    bed_path = f"{baseDir}/{elem['bed-path']}"
    headers = elem['headers']
    window = elem['window']
    site_window = elem['siteWindow']

    bg_inter = elem['bg-intersect-out']
    out_p = elem['dnase-out']

    intersect_winfiles = elem['intersect-out']

    elem_bed = pd.read_csv(bed_path,delimiter="\t",names=headers)
    out_dim = (elem_bed.shape[0],window-site_window+1)

    def post_processing(bam_file):
        outputf = f"{baseDir}/{out_p(bam_file)}"
        if reRun or not os.path.exists(outputf):
            pt.post_process(f"{baseDir}/{bg_inter(bam_file)}", bg_window, headers, lambda X:f"{baseDir}/{intersect_winfiles(bam_file,X)}", out_dim, outputf)
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
    The data is sorted according to ['chrom','name']
    """
    import preprocessing as prep

    args = pt.process_inputArgs(input_parse=sys.argv[1:])
    # args = pt.process_inputArgs(input_parse=['--file_index',0,'--nTasks',52,'--taskType','pWin'])

    assert args.taskType in ["pWin","pIntersect","pPostProcess","eWin","eIntersect","ePostProcess"]

    baseDir = prep.baseDataDir
    p = prep.promoter
    e = prep.enhancer
    bgwin = prep.bgWindow

    bamfs = prep.bamfiles
    interesctOpt = prep.intersectOptions

    reRun = prep.reRun

    if args.taskType=="pWin":
        generate_winBed(args, reRun, baseDir, p)
    elif args.taskType=="pIntersect":
        intersect_elem(args, reRun, baseDir, p, bamfs, interesctOpt, tasktype='all')
    elif args.taskType=="pPostProcess":
        elem_postprocessing(args, reRun, baseDir,  p, bgwin, bamfs)
    elif args.taskType=="eWin":
        generate_winBed(args, reRun, baseDir, e)
    elif args.taskType=="eIntersect":
        intersect_elem(args, reRun, baseDir, e, bamfs, interesctOpt, tasktype='all')
    else:#ePostProcess
        elem_postprocessing(args, reRun, baseDir,  e, bgwin, bamfs)
# python process_Dnase.py --nTasks=600  --file_index=$SLURM_ARRAY_TASK_ID

# sbatch -J prWinBed -q batch --ntasks=1 --cpus-per-task=1 --mem=200G --time 70:00:00 -o slurm_dnase_pr.out --array=0-600 run_script.sh
# sbatch -J prInters -q batch --ntasks=1 --cpus-per-task=1 --mem=200G --time 70:00:00 -o slurm_dnase_pr.out --array=0-600 run_script.sh

# sbatch -J enhWinBed -q batch --ntasks=1 --cpus-per-task=1 --mem=200G --time 70:00:00 -o slurm_dnase_enh.out --array=0-600 run_script.sh
# sbatch -J enhInters -q batch --ntasks=1 --cpus-per-task=1 --mem=200G --time 70:00:00 -o slurm_dnase_enh.out --array=0-600 run_script.sh
















