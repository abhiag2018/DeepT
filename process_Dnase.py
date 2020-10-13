"""
outputs openness score computed from DNase-Seq data for elements.
"""
import colored_traceback.always
import glob
import sys, os, shutil, glob, re
import argparse, itertools
import pandas as pd
import numpy as np
import math
import preptools as pt
import functools
import deeptools.countReadsPerBin as crpb

LIMIT = None

def generateIndex(args, bamfilesInit):
    """
    samtools index <bam file>
    """
    # bamfilesInit = prep.bamfilesInit
    genIdx = lambda f:os.system(f"samtools index {f}")

    args.nTasks = min(args.nTasks,len(bamfilesInit))

    pt.distribute_task(task_list = bamfilesInit, 
        nTasks = int(args.nTasks), file_index=int(args.file_index), 
        func=genIdx, num_tasks=len(bamfilesInit),
        dry_run=False)
    return 0

def convertBAM_to_Array(args, dnaseTmpDir, bamfilesInit, all_chrom, chromLenDict, task_limit = None):
    """
    convert bam file to read count array and save files for each chromosome
    """
    # chromLenDict = prep.chromLen_GRCh37v13
    # dnaseTmpDir = prep.dnaseTmpDir

    tasklist = list(itertools.product(prep.bamfilesInit,['chr'+ch for ch in prep.all_chrom]))[:task_limit]

    saveArray = lambda bam_file, chrom :pt.getBamCounts(bam_file, chrom, chromLenDict[chrom], outputf = f"{dnaseTmpDir(bam_file)}/{chrom}.npz")

    args.nTasks = min(args.nTasks,len(tasklist))

    pt.distribute_task(task_list = tasklist, 
        nTasks = int(args.nTasks), file_index=int(args.file_index), 
        func=lambda t:saveArray(*t), num_tasks=len(tasklist),
        dry_run=False)
    return 0
    
def genDNaseProfile(dnaseTmpDir, bamfilesInit, all_chrom, bed_path, bg_path, headers, bgWin, outputf, limit_frame=None):
    """
    generate Dnase Profile for regulatory element
    """
    # all_chrom = prep.all_chrom
    # headers = prep.promoter['headers']
    # bed_path = prep.promoter['bed-path']
    # bg_path = prep.promoter['bg-path']

    def generate_dnase(tmpBamDir, all_chrom, bed_path, bg_path, headers, bgWin, output_fn, limit_frame):
        chrData = {}
        for c in all_chrom:
            chrData['chr'+c] = tmpBamDir+'/chr'+c+'.npz'
            chrData['chr'+c] = np.load(chrData['chr'+c],allow_pickle=True)['count']

        prDF = pd.read_csv(bed_path, delimiter='\t', names=headers)[:limit_frame]
        prDF_bg = pd.read_csv(bg_path, delimiter='\t', names=headers)[:limit_frame]
        
        dnase_profile = np.array(prDF.apply(lambda df:chrData[df[0]][df[1]:df[2]].flatten(), axis=1).tolist())
        bg_profile = np.array(prDF_bg.apply(lambda df:np.sum(chrData[df[0]][df[1]:df[2]]),axis=1).tolist())

        _dnase_profile = np.sum(dnase_profile,axis=1)
        assert np.sum(_dnase_profile[bg_profile==0])==0
        assert np.logical_not(np.logical_xor((bg_profile<1), (bg_profile==0))).all()
        bg_profile[bg_profile==0]=1
        out = dnase_profile / bg_profile[:,None] * bgWin

        np.savez(f"{tmpBamDir}/{output_fn}.npz",expr = out)
        # np.load(f"{tmpBamDir}/{output_fn}.npz",allow_pickle=True)['expr']
        return 0

    tasklist = [prep.dnaseTmpDir(f) for f in prep.bamfilesInit]

    args.nTasks = min(args.nTasks,len(tasklist))

    objfunc = lambda t:generate_dnase(t, all_chrom, bed_path, bg_path, headers, bgWin, outputf, limit_frame)

    pt.distribute_task(task_list = tasklist, 
        nTasks = int(args.nTasks), file_index=int(args.file_index), 
        func=objfunc, num_tasks=len(tasklist),
        dry_run=False)
    return 0
    

if __name__=="__main__":
    import preprocessing as prep

    args = pt.process_inputArgs(input_parse=sys.argv[1:])
    # args = pt.process_inputArgs(input_parse=['--file_index','0','--nTasks','52','--taskType','pWin'])

    taskTypes = ['genIndex', 'bamToArray' , 'pgenProfile' , 'egenProfile' ]
    assert args.taskType in taskTypes

    p = prep.promoter
    e = prep.enhancer
    bgwin = prep.bgWindow

    chromLenDict = prep.chromLen_GRCh37v13
    dnaseTmpDir = prep.dnaseTmpDir
    bamfsInit = prep.bamfilesInit

    all_chrom = prep.all_chrom

    reRun = prep.reRun

    if args.taskType=="genIndex":
        generateIndex(args, bamfsInit)

    if args.taskType=="bamToArray":
        convertBAM_to_Array(args, dnaseTmpDir, bamfsInit, all_chrom, chromLenDict)

    if args.taskType=="pgenProfile":
        genDNaseProfile(dnaseTmpDir, bamfsInit, all_chrom, p['bed-path'], p['bg-path'], p['headers'], bgwin, 'promoter', limit_frame=None)

    if args.taskType=="egenProfile":
        genDNaseProfile(dnaseTmpDir, bamfsInit, all_chrom, e['bed-path'], e['bg-path'], e['headers'], bgwin, 'enhancer', limit_frame=None)

# python process_Dnase.py --nTasks=600  --file_index=$SLURM_ARRAY_TASK_ID

# sbatch -J prWinBed -q batch --ntasks=1 --cpus-per-task=1 --mem=200G --time 70:00:00 -o slurm_dnase_pr.out --array=0-600 run_script.sh
# sbatch -J prInters -q batch --ntasks=1 --cpus-per-task=1 --mem=200G --time 70:00:00 -o slurm_dnase_pr.out --array=0-600 run_script.sh

# sbatch -J enhWinBed -q batch --ntasks=1 --cpus-per-task=1 --mem=200G --time 70:00:00 -o slurm_dnase_enh.out --array=0-600 run_script.sh
# sbatch -J enhInters -q batch --ntasks=1 --cpus-per-task=1 --mem=200G --time 70:00:00 -o slurm_dnase_enh.out --array=0-600 run_script.sh
















