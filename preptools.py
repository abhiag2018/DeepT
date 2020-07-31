import traceback
import bisect
import math
import numpy as np
import pandas as pd
import csv
import itertools
import re, os, sys, argparse
import glob
import subprocess
import matplotlib.pyplot as plt
import pickle
import matplotlib.pyplot as plt
from sklearn.preprocessing import OneHotEncoder
from numba import jit, prange   
import copy
# import pybedtools

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
        with open(output_path, 'w') as outfile:
            for fname in file_list:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
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

def fasta_to_onehot(fastaf,maxlines=None):
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
    print("executing :bedtools.. ", end="")
    if not dry_run:
        result = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
    else:
        result = 1
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


# def generate_promoter_dna(base_dir, hg19_fa=None, tss_bed = None, out_fa = "promoterTSS_dna.fa", dryRun=False):
#     # base_dir = "/projects/li-lab/agarwa/CUBE/DeepTact/dataset"

#     # input to bedtools : construct combined genome seqeunce file for all chromosomes;
#     if not hg19_fa:
#         dna_seq = lambda x:f"{base_dir}/genome_seq/Homo_sapiens.GRCh37.75.dna.chromosome.{x}.fa"
#         hg19_fa_abs = concat_fasta(fasta_path_fn=dna_seq,re_process=False)
#         hg19_fa = os.path.relpath(hg19_fa_abs, base_dir)

#     # input to bedtools : construct bed file for promoters
#     if not tss_bed:
#         raise NameError("run `import preprocessing` for promoter_preprocessing")

#     return generate_elem_fa(f"{base_dir}/{hg19_fa}",f"{base_dir}/{tss_bed}",out_fa=f"{base_dir}/{out_fa}",dry_run=dryRun)



# def generate_enhancer_dna(base_dir, hg19_fa=None, enhWin_bed = None, out_fa = "enhancer_dna.fa", dryRun=False, options = {}):
#     # base_dir = "/projects/li-lab/agarwa/CUBE/DeepTact/dataset"
#     dna_seq = lambda x:f"genome_seq/Homo_sapiens.GRCh37.75.dna.chromosome.{x}.fa"
#     combined_fasta = dna_seq('chr1-22XY')

#     # input to bedtools : construct combined genome seqeunce file for all chromosomes;
#     if not hg19_fa:
#         dna_seq = lambda x:f"{base_dir}/genome_seq/Homo_sapiens.GRCh37.75.dna.chromosome.{x}.fa"
#         hg19_fa_abs = concat_fasta(fasta_path_fn=dna_seq,re_process=False)
#         hg19_fa = os.path.relpath(hg19_fa_abs, base_dir)


#     # input to bedtools : construct bed file for promoters
#     if not enhWin_bed:
#         raise NameError("run `import parameters` for enhancer_preprocessing")

#     return generate_elem_fa(f"{base_dir}/{hg19_fa}",f"{base_dir}/{enhWin_bed}",out_fa=f"{base_dir}/{out_fa}",dry_run=dryRun)



def concat_PCHiC_PE(base_dataDir,EP_dataPath,selectCell='MK',threshold = 5):
    EP_dataPath_abs = f"{base_dataDir}/{EP_dataPath}"
    promoter_seq = f"{base_dataDir}/promoterTSS_dna.fa"
    enhancer_seq = f"{base_dataDir}/enhancer_dna.fa"
    _,enhLoc,_,_ = fasta_to_onehot(f"{base_dataDir}/promoterTSS_dna.fa")
    _,prLoc,_,_ = fasta_to_onehot(enhancer_seq)

    def extractElem(elems):
        Chr = [re.match('chr([0-9XY])+:([0-9]+)-([0-9]+)',pr).group(1) for pr in elems]
        Start = [int(re.match('chr([0-9XY])+:([0-9]+)-([0-9]+)',pr).group(2)) for pr in elems]
        End = [int(re.match('chr([0-9XY])+:([0-9]+)-([0-9]+)',pr).group(3)) for pr in elems]
        Mid = [(x+y)//2 for x,y in zip(Start,End)]
        return Chr, Start, End, Mid
    prChr, prStart, prEnd, prTSS =  extractElem(prLoc)
    enhChr, enhStart, enhEnd, enhMid =  extractElem(enhLoc)

    EP_pos_train = pd.read_csv(EP_dataPath_abs,delimiter="\t")
    EP_pos_train = EP_pos_train[EP_pos_train[selectCell]>=threshold]
    # EP_pos_train = EP_pos_train[EP_pos_train.apply(lambda df:selectCell in set(df['cellType(s)'].split(',')),axis=1)]
    EP_pos_train['baitPr'] = [np.empty(0,dtype=float)]*len(EP_pos_train)
    EP_pos_train['baitEnh'] = [np.empty(0,dtype=float)]*len(EP_pos_train)
    EP_pos_train['oePr'] = [np.empty(0,dtype=float)]*len(EP_pos_train)
    EP_pos_train['oeEnh'] = [np.empty(0,dtype=float)]*len(EP_pos_train)

    def addElem(interDF,elemChr,elemMid,elemLoc,pr_enh = 'Pr',loci_type='bait'):
        assert pr_enh in ['Pr','Enh']
        col_name = loci_type+pr_enh
        df = pd.DataFrame({'chr':elemChr,'mid':elemMid,'loc':elemLoc})

        elem_chr = df.groupby(['chr'])[['mid','loc']].apply(lambda g: sorted(g.values.tolist(), key = lambda t: t[0])).to_dict()
        
        Start = loci_type+'Start'
        End = loci_type+'End'
        chrom = loci_type+'Chr'
        def appendElems(St,En,elemList):
            _elemList = np.array(elemList)[:,0].astype(np.int32)
            stIdx = bisect.bisect_left(_elemList,St)
            enIdx = bisect.bisect_left(_elemList,En)
            return elemList[stIdx:enIdx]
        def applyFunc(df):
            if df[chrom+'idx'] in elem_chr.keys():
                return appendElems(df[Start],df[End],elem_chr[df[chrom+'idx']])
            return []
        # import pdb
        # pdb.set_trace()
        interDF[chrom+'idx'] = interDF[chrom].str.strip().str[-1]
        interDF[col_name] = interDF.apply(applyFunc,axis=1)

    addElem(EP_pos_train,prChr,prTSS,prLoc,'Pr',loci_type='bait')
    addElem(EP_pos_train,enhChr,enhMid,enhLoc,'Enh',loci_type='bait')
    addElem(EP_pos_train,prChr,prTSS,prLoc,'Pr',loci_type='oe')
    addElem(EP_pos_train,enhChr,enhMid,enhLoc,'Enh',loci_type='oe')

    return EP_pos_train


def training_PE(base_dataDir,EP_dataPath,selectCell=None,threshold = 5,intTypeList=None,output_name='PCHiC'):

    df = concat_PCHiC_PE(base_dataDir,EP_dataPath,selectCell=selectCell)
    baitPr = df.apply(lambda df:len(df['baitPr']),axis=1)
    baitEnh = df.apply(lambda df:len(df['baitEnh']),axis=1)
    oeEnh = df.apply(lambda df:len(df['oeEnh']),axis=1)
    oePr = df.apply(lambda df:len(df['oePr']),axis=1)
    for intType in intTypeList:
        assert intType in ['PE','PP']
        if intType == 'PE':
            EP_df =  df[((baitPr==0) & (baitEnh==1) & (oePr==1) & (oeEnh==0)) | ((baitPr==1) & (baitEnh==0) & (oePr==0) & (oeEnh==1)) ]
        else:
            EP_df =  df[((baitPr==1) & (baitEnh==0) & (oePr==1) & (oeEnh==0))]
        EP_df.to_csv(f"{base_dataDir}/{output_name}_{intType}_{selectCell}.csv",index=False)
    return

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

def process_peaks(a="hg19_promoter_TSS.bed",b="Dnase-Seq/*.bam",out_path=None,out_suffix="_P.bed",options="",max_try = MAXATTEMPTS):
    for _ in range(max_try):
        try:
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














