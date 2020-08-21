import traceback
import bisect
import math
import numpy as np
import pandas as pd
import csv
import itertools
import re, os, sys, shutil, argparse
import glob
import subprocess
import matplotlib.pyplot as plt
import pickle
import matplotlib.pyplot as plt
from sklearn.preprocessing import OneHotEncoder
from numba import jit, prange   
import copy
import gtfparse


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
        assert pr_onehot.shape[1]==1000 or pr_onehot.shape[1]==2000
        assert pr_onehot.shape[2]==4

        np.savez(outp,sequence=pr_onehot.reshape(pr_onehot.shape[0],pr_onehot.shape[1]*pr_onehot.shape[2]),name=pr_name,loc=pr_loc)
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

def concat_PCHiC_PE(hicTSV,promoter_dna,enhancer_dna,selectCell='MK',threshold = 5, outputF=None, sampleFrac=None):
    """
    input : PCHiC tsv file, promoter .fa file, enhancer .fa file
    output : csv file with columns, baitPr, baitEnh, oePr, and oeEnh. 
        Corresponding to promoters and enhancers mid site intersecting with the bait and oe regions
    """

    def extractElem(elems):
        """for each element (promoters or enhancers) extract the chromosome info, start bp index, end bp index, and mid bp index"""
        Chr = [re.match('chr([0-9XY])+:([0-9]+)-([0-9]+)',pr).group(1) for pr in elems]
        Start = [int(re.match('chr([0-9XY])+:([0-9]+)-([0-9]+)',pr).group(2)) for pr in elems]
        End = [int(re.match('chr([0-9XY])+:([0-9]+)-([0-9]+)',pr).group(3)) for pr in elems]
        Mid = [(x+y)//2 for x,y in zip(Start,End)]
        return Chr, Start, End, Mid

    def getChrDict(dna_out_file):
        """
        convert elements list into a dictionary keyed by chromosome; containing the (mid site, name ) tuple lists
        """
        elemData = np.load(dna_out_file,allow_pickle=True)
        elemChr, _, _, elemMid = extractElem(elemData['loc'])
        elemDF = pd.DataFrame({'chr':elemChr,'mid':elemMid, 'name':elemData['name']})
        # dictionary of (mid index, name) for array of elements ; dictionary keys = chromosomes; 
        # the arrays for each chromosome are sorted by mid
        return elemDF.groupby(['chr'])[['mid','name']].apply(lambda g: sorted(g.values.tolist(), key = lambda t: t[0])).to_dict()

    prChrDict =  getChrDict(promoter_dna)
    enhChrDict = getChrDict(enhancer_dna)

    pchicDF = pd.read_csv(hicTSV,delimiter="\t",dtype={'baitChr':str,'oeChr':str})
    pchicDF = pchicDF[pchicDF[selectCell]>=threshold]
    if sampleFrac:
        pchicDF = pchicDF.sample(frac=sampleFrac, axis=0)

    def intersectElem(St,En,elemList):
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
            return intersectElem(df[Start],df[End],elem_chr[df[chrom]])
        return []

    pchicDF['baitPr'] = pchicDF.apply(lambda df:applyFunc(df,'bait',prChrDict),axis=1) 
    pchicDF['baitEnh'] = pchicDF.apply(lambda df:applyFunc(df,'bait',enhChrDict),axis=1) 
    pchicDF['oePr'] = pchicDF.apply(lambda df:applyFunc(df,'oe',prChrDict),axis=1) 
    pchicDF['oeEnh'] = pchicDF.apply(lambda df:applyFunc(df,'oe',enhChrDict),axis=1) 

    if outputF:
        pchicDF.to_pickle(outputF)
        # pd.read_pickle(outputF)
    return pchicDF


def GroupGeneSymbol(gtf,elemTransript, elemMid):
    # gtf = gtfparse.read_gtf(gtfF)


    # find transcripts in gtf database. If not found : ignore!
    gtfTS = gtf[gtf['transcript_id'].isin(elemTransript)].copy()

    # check if duplicates gene_names for transcript_id exist
    dupl = gtfTS.groupby('transcript_id')['gene_name'].apply(lambda x:len(np.unique(x))>1)
    assert sum(dupl)==0

    valDict = {}
    for t,v in zip(elemTransript,elemMid):
        valDict[t] = v
    gtfTS.insert(0,'mid',gtfTS['transcript_id'].apply(lambda x:valDict[x]))

    gs = gtfTS.groupby('transcript_id').first().reset_index()
    uniqL = lambda series:list(np.unique(series))
    gs = gs.groupby('gene_name').agg({'gene_name':uniqL, 'transcript_id':list, 'mid':list})
    return list(zip(gs['gene_name'],gs['mid'],gs['transcript_id']))


# def training_PE(base_dataDir,EP_dataPath,selectCell=None,threshold = 5,intTypeList=None,output_name='PCHiC'):

#     df = concat_PCHiC_PE(base_dataDir,EP_dataPath,selectCell=selectCell)
#     baitPr = df.apply(lambda df:len(eval(df['baitPr'])),axis=1)
#     baitEnh = df.apply(lambda df:len(eval(df['baitEnh'])),axis=1)
#     oeEnh = df.apply(lambda df:len(eval(df['oeEnh'])),axis=1)
#     oePr = df.apply(lambda df:len(eval(df['oePr'])),axis=1)
#     for intType in intTypeList:
#         assert intType in ['PE','PP']
#         if intType == 'PE':
#             EP_df =  df[((baitPr==0) & (baitEnh==1) & (oePr==1) & (oeEnh==0)) | ((baitPr==1) & (baitEnh==0) & (oePr==0) & (oeEnh==1)) ]
#         else:
#             EP_df =  df[((baitPr==1) & (baitEnh==0) & (oePr==1) & (oeEnh==0))]
#         EP_df.to_csv(f"{base_dataDir}/{output_name}_{intType}_{selectCell}.csv",index=False)
#     return

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

def process_inputArgs(input_parse=sys.argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--file_index', help='index of task to execute',type=int)
    parser.add_argument('--nTasks', help='total number of tasks',type=int)
    parser.add_argument('--taskType', help='total number of tasks',type=str)
    args = parser.parse_args(input_parse)

    if args.file_index and args.file_index>args.nTasks:
        sys.exit("file_index > nTasks")
    return args

if __name__=="__main__":
    base_dataDir = "/projects/li-lab/agarwa/CUBE/DeepTact/dataset"
    # decompress_dir_recursive(f"{base_dataDir}/Dnase-Seq/**/**/*.gz")

    # option_pr = {'window':int(1e6),'pr-bed':"hg19_promoter_allFields.bed",'bed-out':"hg19_promoter_BG.bed"}
    # generate_promoter_dna(base_dataDir,hg19_fa="hg19.fa", tss_bed = None, out_fa = "promoterBG_dna.fa", dryRun=True, options = option_pr)














