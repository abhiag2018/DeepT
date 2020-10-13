import colored_traceback.always
import traceback
import re, os, sys, shutil, argparse
from numba import jit, prange   
import copy
import bisect
import math
import itertools
import csv
import pickle
import h5py

import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import glob
import subprocess
from sklearn.preprocessing import OneHotEncoder
import scipy.signal as sg

import gtfparse
import deeptools.countReadsPerBin as crpb


MAXATTEMPTS = 3

@jit(nopython=True, parallel=True) # Set "nopython" mode for best performance, equivalent to @njit
def numbacountparallel_chunk(bs):
    CR = ord('\n')
    c = 0
    for i in prange(len(bs)):
        if bs[i] == CR:
            c += 1
    return c

def numbacountparallel(filename):
    f = open(filename, "rb")
    total = 0
    while True:
        chunk = f.read(1024*1024*10)
        lines = numbacountparallel_chunk(np.frombuffer(chunk, dtype=np.uint8))
        total += lines
        if not chunk:
            break
    return total

def concat_files(file_list,output_path):
    try:
        output = open(output_path, "wb")
        for f in file_list:
            shutil.copyfileobj(open(f, "rb"), output)
        output.close()
    except e:
        print(e)
        return 1
    return 0

def splitCSV(input_path, out_paths, readArgs = {}, writeArgs = {}):
    """
    split an input csv file into N(=len(out_paths)) parts
    out_paths = [f"output_path_{i}.csv" for i in range(N)] #output path name function
    """
    N = len(out_paths)

    df = pd.read_csv(input_path,**readArgs) # reading file
    max_lines = len(df)
    numlines = int(np.ceil(max_lines/N))
    low = np.arange(0,max_lines,numlines)
    high = np.concatenate((low[1:],[max_lines]))

    for i in range(N):
        df_new = df[low[i]:high[i]] # subsetting DataFrame based on index
        df_new.to_csv(out_paths[i],**writeArgs) # output file 
    return 0

def makedirs(dirpath,exist_ok = False):
    """
    Creates all directories in dirpath. If dirpath already exists, old directories are deleted if exist_ok = False, 
        otherwise dirpath is over-written.
    input : 
        dirpath: absolute directory path
        exist_ok: create directory option
    """
    try:
        os.makedirs(dirpath,exist_ok=exist_ok)
    except FileExistsError:
        shutil.rmtree(dirpath)
        os.makedirs(dirpath,exist_ok=exist_ok)
    return 0

def fasta_seq_length(fastaf):
    with open(fastaf,'r') as ff:
        for line in ff: 
            if re.match(">([A-Za-z0-9-:]+):",line):
                pass
            else:
                window = len(line)-1
                break
    return window

def fasta_to_onehot(fastaf,maxlines=None,outp=None):
    maxread = None
    if maxlines:
        numlines = maxlines
        maxread = maxlines*2
    else:
        numlines = numbacountparallel(fastaf)//2
    window = fasta_seq_length(fastaf)

    pr_name = [0]*numlines
    pr_loc = [0]*numlines
    pr_seq = [0]*numlines
    pr_onehot = np.zeros((numlines,window,4))

    enc = OneHotEncoder(handle_unknown='ignore')
    init = np.array(list('ATCG')).reshape(4,1)
    tmp = enc.fit(init)

    num_prom=0
    with open(fastaf,'r') as ff:
        for line in itertools.islice(ff,0,maxread):
            new_pr = re.match(">(.*)::(chr[0-9XY]+:[0-9]+-[0-9]+)",line)
            if new_pr:
                pr_name[num_prom] = new_pr.group(1)
                pr_loc[num_prom] = new_pr.group(2)
                num_prom +=1
            else:
                seq = line[:-1].upper()
                if len(set(list(seq)) - set(list("ATCGN")))!=0:
                    print(seq)
                _seq = np.array(list(seq)).reshape(len(seq),1)
                pr_onehot[num_prom-1] = enc.transform(_seq).toarray()
                pr_seq[num_prom-1] = seq
    if outp:
        # assert pr_onehot.shape[1]==1200 or pr_onehot.shape[1]==2200
        # assert pr_onehot.shape[2]==4

        np.savez(outp,sequence=pr_onehot.reshape(pr_onehot.shape[0],pr_onehot.shape[1]*pr_onehot.shape[2]),name=pr_name,loc=pr_loc)
        # np.load(outp, allow_pickle=True)['sequence']

    return pr_name,pr_loc,pr_seq,pr_onehot

def prep_fasta(chrom=list(range(1,23))+list('XY'),path_lambda=lambda x:f"{x}.fa"):
    for chrid in chrom:
        fasta_path = path_lambda(chrid)
        new_fasta_path = path_lambda('chr'+str(chrid))
        nf = open(new_fasta_path,'w')
        with open(fasta_path,'r') as ff:
            for line in ff:
                if re.match(">",line):
                    tmp=nf.write(re.sub(r">([\dXY])", r">chr\1", line))
                else:
                    tmp=nf.write(line);
        nf.close()

def concat_fasta(fasta_path_fn=lambda x:f"{x}.fa",re_process=True):
    combined_fa_path = dna_seq('chr1-22XY')
    chrom = []
    for i in list(range(1,23))+list('XY'):
        if re_process or (not os.path.exists(f"{dna_seq(i)}")):
            chrom.append(i)
    prep_fasta(chrom=chrom,path_lambda=dna_seq)
    fasta_f = [f"{dna_seq('chr'+str(i))}" for i in list(range(1,23))+list('XY')]
    concat_files(fasta_f,combined_fa_path)
    return combined_fa_path

def generate_elem_fa(hg19_fa,elem_bed,out_fa="out.fa",dry_run=False):
    cmd = f"bedtools getfasta -fi {hg19_fa} -bed {elem_bed} -name -fo {out_fa}"
    print("executing :bedtools", end="..",flush=True)
    if not dry_run:
        result = os.system(cmd)
        # result = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
    print(".",flush=True)
    # print(f"fasta saved to :{out_fa} ")
    return numbacountparallel(out_fa)

def process_promoter_bed(promoter_allfield,out_path,all_headers,window=1000):
    # all_headers = ['chrom', 'txStart', 'txEnd', 'name', 'score', 'strand',
    #        'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts',
    #        'exonEnds']


    all_headers = copy.deepcopy(all_headers)
    loc = all_headers.index("--")
    all_headers.remove("--")

    promoter = pd.read_csv(promoter_allfield,delimiter="\t",header=0).loc[:,all_headers]
    all_chr = ['chr'+str(x) for x in list(range(1,23))+list('XY')]
    promoter = promoter[promoter.apply(lambda df:df['chrom'] in all_chr,axis=1)]
    promoter.insert(loc, 'itemRgb', 0, allow_duplicates=False)
    # func_name = lambda df:df['name']+"::"+df['chrom']+":"+str(df['txStart'])+"-"+str(df['txEnd'])
    # promoter['name'] = promoter.apply(func_name,axis=1)

    negstr = promoter['strand']=='-'
    promoter['txStart'] = promoter['txStart'] - window//2
    promoter.loc[negstr,'txStart'] = promoter.loc[negstr,'txEnd'] - window//2
    promoter['txStart'] = promoter['txStart'].clip(lower=0)
    promoter['txEnd'] = promoter['txStart'] + window

    promoter.to_csv(out_path,sep='\t',header=False,index=False)
    return out_path


def process_enhancer_bed(enh_allfield,out_path,headers,window=2000):
    # headers = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
    #        'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes',
    #        'blockStarts']
    enh = pd.read_csv(enh_allfield,delimiter="\t",names=headers)
    all_chr = ['chr'+str(x) for x in list(range(1,23))+list('XY')]
    enh = enh[enh.apply(lambda df:df['chrom'] in all_chr,axis=1)]

    enh['chromStart'] = enh['thickStart']-window//2
    enh['chromStart'] = enh['chromStart'].clip(lower=0)
    enh['chromEnd'] = enh['chromStart'] + window

    enh.to_csv(out_path,sep='\t',header=False,index=False)
    return out_path

def getChrDict(dna_out_file):
    """
    convert elements list into a dictionary keyed by chromosome; containing the (mid site, name ) tuple lists
    """
    def extractElem(elems):
        """for each element (promoters or enhancers) extract the chromosome info, start bp index, end bp index, and mid bp index"""
        matchObj = [re.match('chr([0-9XY]+):([0-9]+)-([0-9]+)',pr) for pr in elems]
        Chr = [obj.group(1) for obj in matchObj]
        Start = [int(obj.group(2)) for obj in matchObj]
        End = [int(obj.group(3)) for obj in matchObj]
        Mid = [(x+y)//2 for x,y in zip(Start,End)]
        return Chr, Start, End, Mid

    elemData = np.load(dna_out_file,allow_pickle=True)
    elemChr, _, _, elemMid = extractElem(elemData['loc'])
    elemDF = pd.DataFrame({'chr':elemChr,'mid':elemMid, 'name':elemData['name']})
    # dictionary of (mid index, name) for array of elements ; dictionary keys = chromosomes; 
    # the arrays for each chromosome are sorted by mid
    return elemDF.groupby(['chr'])[['mid','name']].apply(lambda g: sorted(g.values.tolist(), key = lambda t: t[0])).to_dict()


def concat_PCHiC_PE(hicTSV,promoter_dna,enhancer_dna,selectCell='MK',threshold = 5, outputF=None, sampleFrac=None):
    """
    input : PCHiC tsv file, promoter .fa file, enhancer .fa file
    output : csv file with columns, baitPr, baitEnh, oePr, and oeEnh. 
        Corresponding to promoters and enhancers mid site intersecting with the bait and oe regions
    """

    prChrDict =  getChrDict(promoter_dna)
    enhChrDict = getChrDict(enhancer_dna)

    pchicDF = pd.read_csv(hicTSV,delimiter="\t",dtype={'baitChr':str,'oeChr':str})
    pchicDF = pchicDF[pchicDF[selectCell]>=threshold]
    if sampleFrac:
        pchicDF = pchicDF.sample(frac=sampleFrac, axis=0)

    def intersectSortedElem(St,En,elemList):
        """
        returns all elements in elemList lying between St and En 
        i.e. element x in output array IFF St<=x<En
        """
        _elemList = np.array(elemList)[:,0].astype(np.int32)
        stIdx = bisect.bisect_left(_elemList,St)
        enIdx = bisect.bisect_left(_elemList,En)
        return elemList[stIdx:enIdx]

    def applyFunc(df,loci_type,elem_chr):
        Start = loci_type+'Start'
        End = loci_type+'End'
        chrom = loci_type+'Chr'
        # pdb.set_trace()
        if df[chrom] in elem_chr.keys():
            return intersectSortedElem(df[Start],df[End],elem_chr[df[chrom]])
        return []

    pchicDF['baitPr'] = pchicDF.apply(lambda df:applyFunc(df,'bait',prChrDict),axis=1) 
    pchicDF['baitEnh'] = pchicDF.apply(lambda df:applyFunc(df,'bait',enhChrDict),axis=1) 
    pchicDF['oePr'] = pchicDF.apply(lambda df:applyFunc(df,'oe',prChrDict),axis=1) 
    pchicDF['oeEnh'] = pchicDF.apply(lambda df:applyFunc(df,'oe',enhChrDict),axis=1) 

    if outputF:
        print(f"saving converted file to {outputF}..",end=" ",flush=True)
        pchicDF.to_pickle(outputF)
        # pd.read_pickle(outputF)
        print(".",flush=True)
    return pchicDF


def GroupGeneSymbol(gtf, elemTransript, elemMid):
    """
    input:
        gtf: gene annotation file
        elemTransript: transcripts list for a single PCHiC interaction 
        elemMid: some data corresponding to the trancripts (in the same order in the list)
    output : 
        grouped (gene, transcripts, elemMid)
    """
    # gtf = gtfparse.read_gtf(gtfF)
    # find transcripts in gtf database. If not found : ignore!
    gtfTS = gtf[gtf['transcript_id'].isin(elemTransript)].copy()

    # check if duplicates gene_names for transcript_id exist in gtf file 
    dupl = gtfTS.groupby('transcript_id')['gene_name'].apply(lambda x:len(np.unique(x))>1)
    assert sum(dupl)==0 # they shouldn't

    valDict = {}
    for t,v in zip(elemTransript,elemMid):
        valDict[t] = v
    gtfTS.insert(0,'mid',gtfTS['transcript_id'].apply(lambda x:valDict[x]))

    gs = gtfTS.groupby('transcript_id').first().reset_index() # since there are no duplicates this statement is useless
    uniqL = lambda series:list(np.unique(series))
    gs = gs.groupby('gene_name').agg({'gene_name':uniqL, 'transcript_id':list, 'mid':list})
    return list(zip(gs['gene_name'],gs['mid'],gs['transcript_id']))

def decompress_dir_recursive(glob_expr):
    for file in glob.glob(glob_expr):
        print("converting", os.path.basename(file),"..",end=" ")
        cmd = f"gzip -d {file}"
        result = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
        print(".")
        print(result,end="\n\n")

def distribute_task(task_list = [],nTasks = 10, file_index = 0, func=None,num_tasks=None,dry_run=False):
    if not num_tasks:
        num_tasks = len(task_list)

    nJobs_per_task = math.ceil(num_tasks/nTasks)
    start_idx = file_index*nJobs_per_task
    end_idx = min((1+file_index)*nJobs_per_task,num_tasks)

    print(f"processing ({start_idx}-{end_idx})/{num_tasks}.")
    if start_idx<end_idx:
        for dcmf in itertools.islice(task_list,start_idx,end_idx):
            if not dry_run:
                func(dcmf)
    else:
        print('exit',file_index,flush=True)
    return 0

def replace_suffix(f,ext=".bed"):
    basename = os.path.basename(f)
    filen = re.compile(r"^(.*)\..*$").match(basename)
    if filen:
        basename = filen.group(1)
    if os.path.dirname(f)=="":
        out_path = f"{basename}{ext}"
    else:
        out_path = f"{os.path.dirname(f)}/{basename}{ext}"
    return out_path

def process_peaks(a="hg19_promoter_TSS.bed",b="Dnase-Seq/*.bam",out_path=None,out_suffix="_P.bed",options="",max_try = MAXATTEMPTS, recompute=True):
    for _ in range(max_try):
        try:
            if not recompute and out_path and os.path.exists(out_path):
                result=2
                break

            # run algorithm
            if not out_path:
                out_path = replace_suffix(b,ext=out_suffix)
            cmd = f"bedtools intersect {options} -a {a} -b {b} > {out_path}"

            print("processing",os.path.basename(a), " & ", os.path.basename(b), end =".. ",flush=True)
            result = os.system(cmd)
            # result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
            print(".",flush=True)

            # check if output is correct
            linecount_a = numbacountparallel(a)
            if os.path.exists(out_path):
                linecount_out = numbacountparallel(out_path)
            else:
                continue
            if linecount_out<linecount_a:
                continue


            break
        except Exception:
            print(cmd,"::",result)
            print(traceback.print_exc())
            result = 1
    return result


def split_bam(bam_file,options="-reference"):
    cmd = f"bamtools split -in {bam_file} {options}"
    result = os.system(cmd)
    return result

def generate_elem_window_bed(elemBed_path, Head, winIdx, site_window, ouputfile):
    exit_status = 0
    st_col = 1
    end_col = 2
    try:
        elem_bed = pd.read_csv(elemBed_path,delimiter="\t",names=Head)

        init_start = elem_bed.iloc[:, st_col]
        init_end = elem_bed.iloc[:, end_col] 

        elem_bed.iloc[:, st_col] = init_start+winIdx
        elem_bed.iloc[:, end_col] = elem_bed.iloc[:, st_col] + site_window
            
        elem_bed.to_csv(ouputfile, sep='\t', header=False, index=False)
    except Exception:
        print(traceback.print_exc())
        exit_status = 1
    return exit_status


def combine_intersectBed(bedfiles, Head, outputf, remove=False):
    all_headers = Head + ['bamReads']
    try:
        X = pd.Series(pd.read_csv(bedfiles[0],delimiter='\t', names=all_headers)['bamReads'])
        # X[X == '+']=0
        X = X.apply(int)

        for bedf in itertools.islice(bedfiles,1,None):
            bedDF = pd.read_csv(bedf,delimiter='\t', names=all_headers)
            tmpX =  pd.Series(bedDF['bamReads'])
            # tmpX[tmpX == '+']=0
            tmpX = tmpX.apply(int)
            X += tmpX
        bedDF['bamReads'] = X
    except Exception as e:
        print(bedf)
        raise Exception(bedf)
    bedDF.to_csv(outputf,sep='\t',header=False,index=False)
    if remove:
        for f in bedfiles:
            os.remove(f)
    return outputf

def combine(input_files, hic_out):
    """
    combine pandas DataFrames  into another pandas DataFrame and save pickled output.
    input :
        input_files: list of *pickled* pandas DataFrame of the same format
        hic_out : pickled output path
    """
    print(f"saving converted file to {hic_out}..",end=" ",flush=True)
    pd.concat([pd.read_pickle(inp) for inp in input_files]).to_pickle(hic_out)
    print(".",flush=True)
    return 0

def getBamCounts(bam_file, chrom, chromLen, outputf = None):
    """
    read a bam file and return the readCounts for position - chrom:0-chromLen
    inputs: 
        chrom : chromosome to read
        chromLen : scalar int
    """
    cr = crpb.CountReadsPerBin([bam_file], binLength=1, stepSize=1)
    arr = cr.count_reads_in_region(chrom, 0, chromLen)[0]
    if outputf:
        np.savez(outputf, count=arr)
        # np.load(outputf,allow_pickle=True)['count']
    return arr

def post_process(bg_intersect_out, bg_win, Head, intersect_out, out_dim, out_path):
    """
    combine intersectBed output into numpy array
    """
    # intersect_out = lambda X:eval("f'{}'".format(intersectBed_fStr[0]))
    intersectBed_files = [intersect_out(i) for i in range(out_dim[1])]

    all_headers = Head + ['bamReads']
    select_elem = lambda df,elem,elem_val: df.loc[df[elem] == elem_val]['bamReads'].item()


    Dnase_reads_bg = pd.read_csv(bg_intersect_out, delimiter='\t', names=all_headers)
    Dnase_reads_bg = Dnase_reads_bg.sort_values(by=['chrom','name'],axis=0)
    Dnase_data = np.zeros(out_dim)

    for i,intersectf in enumerate(intersectBed_files):
        Dnase_reads_win = pd.read_csv(intersectf, delimiter='\t', names=all_headers)
        Dnase_reads_win = Dnase_reads_win.sort_values(by=['chrom','name'],axis=0)

        Dnase_reads_win['bgBamReads'] = Dnase_reads_bg['bamReads']
        # Dnase_reads_win['bgBamReads'] = Dnase_reads_win.apply(lambda df: select_elem(Dnase_reads_bg,'name',df['name']),axis=1)

        Dnase_reads_win['opennessScore'] = Dnase_reads_win['bamReads']/Dnase_reads_win['bgBamReads']*bg_win

        Dnase_data[:,i] = Dnase_reads_win['opennessScore'].to_numpy(dtype=float)
    np.savez(out_path,expr=Dnase_data)
    return out_path

def process_inputArgs(input_parse=sys.argv, argType = {'file_index':None, 'nTasks':None, 'taskType':None } ):

    if 'file_index' in argType.keys() and argType['file_index'] is None:
        argType['file_index'] = {'type':int,  'help':'index of task to execute'}

    if 'nTasks' in argType.keys() and argType['nTasks'] is None:
        argType['nTasks'] = {'type':int, 'help':'total number of tasks'}

    if 'taskType' in argType.keys() and argType['taskType'] is None:
        argType['taskType'] = {'type':str, 'help':'task type'}

    if 'cellType' in argType.keys() and argType['cellType'] is None:
        argType['cellType'] = {'type':str, 'help':"cell type in ['tB','tCD4','nCD4','FoeT','Mon','tCD8']"}


    parser = argparse.ArgumentParser()
    for t in argType:
        parser.add_argument(f"--{t}", **argType[t])

    args = parser.parse_args(input_parse)

    if 'file_index' in argType.keys() and args.file_index and args.file_index>args.nTasks:
        sys.exit("file_index > nTasks")
    return args

def SelectElements_in_TrainingData(cell_trainData, pr_data, enh_data, all_prList, all_enhList, prout, enhout, savez_key, transform = lambda x:x, skip_assert = False):
    """
    check if elements in training data are present in the element list. If yes, then select those elements, and save the inpur features corresponding to those elements
    input : 
        cell_trainData : training label csv file
        pr_data : input features for promoter
        enh_data : input features for promoter
        all_prList : all promoters list
        all_enhList : all enhancers list
        prout : selected promoters' features output file 
        enhout : selected enhancers' features output file 
        savez_key : output file key
    """

    trainData = pd.read_csv(cell_trainData)
    train_promoters = trainData['promoter_name'].apply(lambda x:x.split(":")[0])
    train_enhancers = trainData['enhancer_name']

    tfunc = lambda x: pd.Series([x.split(':')[0]]+list(map(int, x.split(':')[1].split('-'))))
    all_enh_DF = all_enhList.apply(tfunc).rename(columns={0:'chr',1:'st',2:'en'})


    if not skip_assert:
        # make sure you for each promoter in training data you find unique match in promoter list
        Count_Matched_Promoters = lambda x:sum(all_prList==x)
        assert len(train_promoters) == sum(train_promoters.apply(Count_Matched_Promoters))
        assert len(train_promoters) == sum(train_promoters.apply(Count_Matched_Promoters)>=1)

        # make sure you for each enhancer in training data you find unique match in enhancer list
        trainenhMid = train_enhancers.apply(lambda x:sum(map(int,x.split(':')[1].split('-')))//2)
        trainenhChr = train_enhancers.apply(lambda x:x.split(':')[0])
        train_enh_DF = pd.DataFrame({'mid':trainenhMid,'chr':trainenhChr})
        count_match =  train_enh_DF.apply(lambda x: sum((all_enh_DF['st']<x['mid']) & (x['mid']<=all_enh_DF['en']) & (x['chr']==all_enh_DF['chr'])),axis=1)
        assert len(train_enhancers) == sum(count_match>=1)
        assert len(train_enhancers) == sum(count_match)


    # # select the DNase data for promoter corresponding to the Training Data
    # pr_data = np.load(f"{bam_dir}/promoter.npz",allow_pickle=True)['expr']
    # enh_data = np.load(f"{bam_dir}/enhancer.npz",allow_pickle=True)['expr']

    def select_from_extPr(data, all_prList, row):
        bin_sel =  all_prList==row['promoter_name'].split(":")[0]
        stidx = row['promoter_start']-row['pr_base_start']
        elem_len = row['promoter_end']-row['promoter_start']
        return data[bin_sel][0][stidx:stidx+elem_len]

    def select_from_extEnh(data, all_enh_DF, row):
        row_chr = row['enhancer_name'].split(':')[0]
        row_mid = sum(map(int,row['enhancer_name'].split(':')[1].split('-')))//2

        bin_sel =  (all_enh_DF['st']<=row_mid) & (row_mid<=all_enh_DF['en']) & ( row_chr==all_enh_DF['chr'])
        stidx = row['enhancer_start']-row['enh_base_start']
        elem_len = row['enhancer_end']-row['enhancer_start']

        return data[bin_sel][0][stidx:stidx+elem_len]

    dnase_pr_train = trainData.apply(lambda row: select_from_extPr(pr_data, all_prList, row),axis=1)
    dnase_pr_train = np.array(dnase_pr_train.tolist())

    dnase_enh_train = trainData.apply(lambda row: select_from_extEnh(enh_data, all_enh_DF, row),axis=1)
    dnase_enh_train = np.array(dnase_enh_train.tolist())
    # dnase_enh_train =  train_enh_DF.apply(lambda x: enh_data[(all_enh_DF['st']<=x['mid']) & (x['mid']<=all_enh_DF['en']) & (x['chr']==all_enh_DF['chr'])][0],axis=1)

    np.savez(prout, **{savez_key:transform(dnase_pr_train)})
    np.savez(enhout, **{savez_key:transform(dnase_enh_train)})

    # np.load(f"{dirDNase}/enhancerTrain.npz", allow_pickle=True)['expr']
    return 0


def npz_to_hdf5(npz_fpath, transform=lambda x:x):
    """save npz file as hdf5"""
    path_list = npz_fpath.split("/")
    out_path = '/'.join(path_list[:-1])+'/'+path_list[-1].split('.')[0]+'.hdf5'
    npz_file = np.load(npz_fpath, allow_pickle=True)

    with h5py.File(out_path, "w") as data_file:
        for k in npz_file.files:
            data_file.create_dataset(k, data=transform(npz_file[k]))
    print(f"{out_path} saved.")


flatten = lambda x:list(itertools.chain.from_iterable(x))

if __name__=="__main__":
    base_dataDir = "/projects/li-lab/agarwa/CUBE/DeepTact/dataset"
    # decompress_dir_recursive(f"{base_dataDir}/Dnase-Seq/**/**/*.gz")

    # option_pr = {'window':int(1e6),'pr-bed':"hg19_promoter_allFields.bed",'bed-out':"hg19_promoter_BG.bed"}
    # generate_promoter_dna(base_dataDir,hg19_fa="hg19.fa", tss_bed = None, out_fa = "promoterBG_dna.fa", dryRun=True, options = option_pr)














