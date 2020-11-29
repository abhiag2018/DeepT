import colored_traceback.always
import sys, os, argparse
import random, itertools

import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import gtfparse

import preptools as pt

SUBSAMPLE = None

def HiCMatch(args, HiCParts, tmp_out, promoter_dna, enhancer_dna, cell, th, sampleFrac = SUBSAMPLE):
    """
    wrapper for matching PCHiC data in tsv file to promoter and enhancer lists
    """
    # HiCParts = prep.HiCParts
    # codeTmpDir = prep.codeTmpDir
    numTasks = len(HiCParts)
    tasklist = list(zip(HiCParts,tmp_out))

    def match_pcHiC(task):
        hicTSV = task[0]
        outpkl = task[1]
        pt.concat_PCHiC_PE(hicTSV,promoter_dna, enhancer_dna, selectCell=cell, threshold = th, outputF=outpkl, sampleFrac=sampleFrac)

    args.nTasks = min(args.nTasks,numTasks)
    print(f"matching PCHiC Data to promoter & enhancers (file_index = {args.file_index})..", end=" ", flush=True)
    pt.distribute_task(task_list = tasklist, 
        nTasks = int(args.nTasks), file_index=int(args.file_index), 
        func=match_pcHiC, num_tasks=numTasks,
        dry_run=False)
    print(".", flush=True)
    return 0


def transcriptsToGene(gtfFile, hicMatched, hic_out):
    """
    Convert Promoter Transcript IDs in matches from PCHiC file to gene-symbols. And save output as a pickled pandas DataFrame. The function uses the matching to group Promoters and select a unique promoter per gene-symbol.
    Note: there may still be multiple gene symbols per interaction region (bait/oe). But the mapping from transcript IDs to gene symbols is a function. Note the assert statement in pt.GroupGeneSymbol
    input : 
        gtfFile: the input gtf file from parameters.py
        hicMatched: output of the function self.HiCMatch
        hic_out: pickled output path
    """
    pchicDF = pd.read_pickle(hicMatched)

    pchicDF['baitPr'] =  pchicDF['baitPr'].apply(np.array)
    pchicDF['baitEnh'] =  pchicDF['baitEnh'].apply(np.array)
    pchicDF['oePr'] =  pchicDF['oePr'].apply(np.array)
    pchicDF['oeEnh'] =  pchicDF['oeEnh'].apply(np.array)

    gtf = gtfparse.read_gtf(gtfFile)
    # group Transcripts by Gene Symbols
    groupTS = lambda elemBaitPr:pt.GroupGeneSymbol(gtf,elemBaitPr[:,1],elemBaitPr[:,0]) if elemBaitPr.size>0 else []
    
    print("converting transcripts IDs in  PCHiC to GeneSymbols..", end=" ", flush=True)
    pchicDF['baitPr'] = pchicDF['baitPr'].apply(groupTS)
    print("..", end=" ", flush=True)
    pchicDF['oePr'] = pchicDF['oePr'].apply(groupTS)
    print('.',flush=True)

    print(f"saving converted file to {hic_out}..",end=" ",flush=True)
    pchicDF.to_pickle(hic_out)
    # pd.read_pickle(hic_out)
    print(".",flush=True)
    return 0

def hicUniqueMatch(hicGroupMatched, hic_out_PE=None, hic_out_EP=None):
    """
    select only those interactions where both the elements (bait promoter & oe enhancer; bait enhancer & oe promoter) 
    have a unique match to gene symbols and enhancer
    input :  
        hicGroupMatched: the input file from the previous processing step (processing order in __main__ )
        hic_out_PE : output for promoter = bait; enhancer = oe
        hic_out_EP  : output for promoter = oe; enhancer = bait
    """
    pchicDF = pd.read_pickle(hicGroupMatched)

    baitprF = pchicDF['baitPr'].apply(lambda x:np.array(x).size>0)
    baitenhF = pchicDF['baitEnh'].apply(lambda x:np.array(x).size>0)
    oeprF = pchicDF['oePr'].apply(lambda x:np.array(x).size>0)
    oeenhF = pchicDF['oeEnh'].apply(lambda x:np.array(x).size>0)

    pehic = pchicDF[baitprF&oeenhF].copy()
    ephic = pchicDF[baitenhF&oeprF].copy()

    def uniqMatchPr(st,en,elem):
        MID = (st+en)//2
        
        _mat = []
        for m in elem:
            minidx = min(range(len(m[1])),key=lambda i:int(m[1][i])-MID)
            _mat.append([m[0][0],m[1][minidx],m[2][minidx]])
        return _mat

    # def uniqMatchEnh(st,en,elem):
    #     MID = (st+en)//2
    #     minidx = min(range(len(elem)),key=lambda i:int(elem[i][0])-MID)
    #     return [elem[minidx][0],elem[minidx][1]]

    pehic['baitPr'] = pehic.apply(lambda df:uniqMatchPr(df['baitStart'],df['baitEnd'],df['baitPr']),axis=1)
    ephic['oePr'] = ephic.apply(lambda df:uniqMatchPr(df['oeStart'],df['oeEnd'],df['oePr']),axis=1)

    # pehic['oeEnh'] = pehic.apply(lambda df:uniqMatchEnh(df['oeStart'],df['oeEnd'],df['oeEnh']),axis=1)
    # ephic['baitEnh'] = ephic.apply(lambda df:uniqMatchEnh(df['baitStart'],df['baitEnd'],df['baitEnh']),axis=1)

    pehic = pehic[pehic['baitPr'].apply(len)==1]
    pehic = pehic[pehic['oeEnh'].apply(len)==1]

    ephic = ephic[ephic['baitEnh'].apply(len)==1]
    ephic = ephic[ephic['oePr'].apply(len)==1]

    def apply_cols(df,cols,func):
        for col in cols:
            df[col] = df[col].apply(func)
        return 0
    apply_cols(pehic,['baitPr','oeEnh'],lambda x:np.array(x).flatten())
    apply_cols(ephic,['oePr','baitEnh'],lambda x:np.array(x).flatten())

    if hic_out_PE:
        print(f"saving converted file to {hic_out_PE}..",end=" ",flush=True)
        pehic.to_pickle(hic_out_PE)
        print(".",flush=True)
    if hic_out_EP:
        print(f"saving converted file to {hic_out_EP}..",end=" ",flush=True)
        ephic.to_pickle(hic_out_EP)
        print(".",flush=True)

    return pehic,ephic


def HiCTrainingPosLabel(hicUniquePE, hicUniqueEP, prwin, enhwin, traindata_out = None):
    """
    select the data from the columns of the training data and put them in appropriate format
    for input to DeepTACT. The output contains only positive labels.
        traindata_out: csv output file path
    if traindata_out is None the output is not saved to a file
    """
    def niceHiCDF(pickle_path, prwin, enhwin, prefixPr, prefixEnh):
        pchicDF = pd.read_pickle(pickle_path)
        # pchicDF = pd.read_pickle(f"/projects/li-lab/agarwa/CUBE/DeepTact/code/tmp_data/PCHiC_{cell}_transcript.pkl")
        pchicDF['oeChr'] = 'chr'+pchicDF['oeChr']
        pchicDF['baitChr'] = 'chr'+pchicDF['baitChr']
        pchicDF['baitPr'] = pchicDF['baitPr'].apply(np.array)
        pchicDF['oePr'] = pchicDF['oePr'].apply(np.array)

        train_out = {}

        train_out['enhancer_chrom'] = pchicDF[prefixEnh+'Chr']
        train_out['enhancer_start'] = pchicDF[prefixEnh+'Enh'].apply(lambda x:int(x[0])-enhwin//2)
        train_out['enhancer_end'] = pchicDF[prefixEnh+'Enh'].apply(lambda x:int(x[0])+enhwin//2)
        train_out['enhancer_name'] = pchicDF[prefixEnh+'Enh'].apply(lambda x:x[1])


        train_out['promoter_chrom'] = pchicDF[prefixPr+'Chr']
        train_out['promoter_start'] = pchicDF[prefixPr+'Pr'].apply(lambda x:int(x[1])-prwin//2)
        train_out['promoter_end'] = pchicDF[prefixPr+'Pr'].apply(lambda x:int(x[1])+prwin//2)
        train_out['promoter_name'] = pchicDF[prefixPr+'Pr'].apply(lambda x:x[2]+":"+x[0])

        train_out['label'] = 1
        return pd.DataFrame(train_out)


    dfPE = niceHiCDF(hicUniquePE, prwin, enhwin, prefixPr="bait", prefixEnh="oe")
    dfEP = niceHiCDF(hicUniqueEP, prwin, enhwin, prefixPr="oe", prefixEnh="bait")

    trainDF =  pd.concat([dfPE,dfEP])
    if traindata_out:
        print(f"saving converted file to {traindata_out}..",end=" ",flush=True)
        trainDF.to_csv(traindata_out, index=False)
        print(".",flush=True)
    return trainDF


def intersectHiCElem(pchicDF, elemMid, chrom, intersectWith = 'bait'):
    """
    helper function to HiCTrainingNegLabel(..)
    """

    hicDF = pchicDF[pchicDF[intersectWith+'Chr'] == chrom]

    assert intersectWith in ['bait','oe']
    if intersectWith=='bait':
        intersectOther = 'oe'
    else:
        intersectOther = 'bait'

    st = intersectWith+'Start'
    en = intersectWith+'End'

    stOth = intersectOther+'Start'
    enOth = intersectOther+'End'

    
    intersectedRows =  hicDF[(hicDF[st]<=elemMid) & (elemMid<hicDF[en])]
    return intersectedRows[stOth], intersectedRows[enOth]

# function for checking whether there is any element within the range of start and end values
checkInter =lambda st,en,elem:sum((st<=elem) & (elem<en))>0

def HiCTrainingNegLabel(hicTSV, cell, th, pos_training, promoter_dna, enhancer_dna, prwin, enhwin,hic_out = None, numSamples=None):
    """
    add negative training labels according to histogram of positive training labels
    ignore those randomly chosen negative elements that have threshold higher than th. 
    input :
        th: selection threshold for Hi-C interaction. Above which it is not considered for a negative label. (should be 0)
    """
    # pos_training = prep.HiC_Training('MK')
    # promoter_dna = prep.promoter['dna-out']
    # enhancer_dna = prep.enhancer['dna-out']
    pchicDF = pd.read_csv(hicTSV,delimiter="\t",dtype={'baitChr':str,'oeChr':str})
    pchicDF = pchicDF[pchicDF[cell]>=th]


    LINF = -2e6 #2 Mbp is the max distance
    RINF = 2e6
    side_prob = 0.1

    trainData = pd.read_csv(pos_training)
    TrainLbl = {'enhancer_chrom': str, 'enhancer_start':int, 'enhancer_end':int, 'enhancer_name':str, 'promoter_chrom':str, 'promoter_start':int, 'promoter_end':int, 'promoter_name':str, 'label':int}
    trainData = trainData.astype(TrainLbl)
    if numSamples is None:
        numSamples = trainData.shape[0]
    else:
        numSamples = trainData.shape[0]*numSamples

    _freq, _bins, _ = plt.hist((trainData.enhancer_start + trainData.enhancer_end)//2  - (trainData.promoter_start + trainData.promoter_end)//2)
    bins = np.concatenate(([LINF], _bins, [RINF]))
    freq = np.concatenate(([side_prob/2], _freq/sum(_freq), [side_prob/2]))

    prDict = pt.getChrDict(promoter_dna)
    enhDict = pt.getChrDict(enhancer_dna)

    all_chrom = list([str(x) for x in range(1,23)]) + list('XY')

    promoters = []
    enhancers=[]
    chroms = []
    for _ in range(numSamples):
        chrom = np.random.choice(all_chrom)

        prMid = pd.Series(sorted(prDict[chrom],key=lambda x:x[0])).apply(lambda x:[int(x[0]),x[1]])
        enhMid = pd.Series(sorted(enhDict[chrom],key=lambda x:x[0])).apply(lambda x:[int(x[0]),x[1]])

        promoter = np.random.choice(prMid)

        p = np.zeros(len(enhMid))
        part = np.searchsorted(enhMid.apply(lambda x:x[0]).to_numpy() - promoter[0], bins)
        
        for binID in range(1,len(part)):
            if part[binID]-part[binID-1] > 0:
                val = freq[binID-1]/(part[binID]-part[binID-1])
                p[part[binID-1]:part[binID]] = val
        p = p/sum(p)

        _st1, _en1 = intersectHiCElem(pchicDF, promoter[0], chrom, intersectWith = 'bait')
        _st2, _en2 = intersectHiCElem(pchicDF, promoter[0], chrom, intersectWith = 'oe')
        enhancer = np.random.choice(enhMid,p=p)
        while checkInter(_st1,_en1,enhancer[0]) or checkInter(_st2,_en2,enhancer[0]):
            print(".",end="",flush=True)
            enhancer = np.random.choice(enhMid,p=p)

        chroms.append('chr'+chrom)
        promoters.append(promoter)
        enhancers.append(enhancer)

    print(f"generated {numSamples} negative samples.",flush=True)

    prMids = pd.Series(np.array(promoters)[:,0]).apply(int).to_numpy()
    prNames = np.array(promoters)[:,1]
    enhMids = pd.Series(np.array(enhancers)[:,0]).apply(int).to_numpy()
    enhNames = np.array(enhancers)[:,1]

    negTrain={'enhancer_chrom': chroms, 'enhancer_start': enhMids - enhwin//2, 'enhancer_end': enhMids + enhwin//2, 'enhancer_name': enhNames, 
        'promoter_chrom': chroms, 'promoter_start': prMids - prwin//2, 'promoter_end': prMids + prwin//2, 'promoter_name': prNames, 
        'label': [0]*numSamples}

    negTrain = pd.DataFrame(data=negTrain)
    negTrain = negTrain.astype(TrainLbl)
    Data = pd.concat([trainData,negTrain])
    if hic_out:
        print(f"saving converted file to {hic_out}..",end=" ",flush=True)
        Data.to_csv(hic_out,index=False)
        print(".",flush=True)
    return Data

def train_augment(cell_traindata, cell_train_aug, aug_len_enh, aug_len_pr, aug_step_enh, aug_step_pr, enh_length, pr_length, mult_fac = 20):
    """
    augment training data csv
    """
    # cell_traindata = prep.hictrain(cell)
    # cell_train_aug = prep.hictrain_augment(cell)
    # aug_len_enh = prep.enhancer['augment_len']
    # aug_len_pr = prep.promoter['augment_len']
    # aug_step_enh = prep.enhancer['aug_step']
    # aug_step_pr = prep.promoter['aug_step']

    # enh_length = prep.enhancer['window']
    # pr_length = prep.promoter['window']

    def genDF(index,row):
        enh_shifts = random.choices(np.arange(0,aug_len_enh,aug_step_enh),k=mult_fac)
        pr_shifts = random.choices(np.arange(0,aug_len_pr,aug_step_pr),k=mult_fac)

        return pd.DataFrame(
            {'enhancer_chrom': row.enhancer_chrom,
            'enhancer_start': [row.enhancer_start+x for x in enh_shifts],
            'enhancer_end': [row.enhancer_start+enh_length+x for x in enh_shifts],
            'enhancer_name': row.enhancer_name,
            'promoter_chrom': row.promoter_chrom,
            'promoter_start': [row.promoter_start+x for x in pr_shifts],
            'promoter_end': [row.promoter_start+pr_length+x for x in pr_shifts],
            'promoter_name': row.promoter_name,
            'label': row.label,
            'enh_base_start': row.enhancer_start,
            'pr_base_start': row.promoter_start
            })

    random.seed(42)
    train_data = pd.read_csv(cell_traindata)

    train_augDF = pd.concat(map(lambda x:genDF(*x), train_data.iterrows()), ignore_index=True)
    train_augDF.to_csv(cell_train_aug, index=False)
    return 0

def SelectElements_DNase(args, cell_trainData, bamfiles_cell, dnaseTmpDir, all_prList, all_enhList, outname = "Train", skip_assert = False):
    """
    select DNase elements
    """
    # cell_trainData = prep.hictrain_augment(cell)
    # bamfiles_cell = list(itertools.chain.from_iterable(prep.DnaseCells[cell]))
    # dnaseTmpDir = prep.dnaseTmpDir
    # all_prList = pd.read_csv(prep.promoter['bed-path'], delimiter='\t', names=prep.promoter['headers']).name
    # all_enhList = pd.read_csv(prep.enhancer['bed-path'], delimiter='\t', names=prep.enhancer['headers']).name

    tasklist = [dnaseTmpDir(bamf) for bamf in bamfiles_cell]
    objfunc = lambda bam_dir:pt.SelectElements_in_TrainingData(cell_trainData, np.load(f"{bam_dir}/promoter.npz",allow_pickle=True)['expr'],
        np.load(f"{bam_dir}/enhancer.npz",allow_pickle=True)['expr'],
        all_prList, all_enhList, 
        f"{bam_dir}/promoter{outname}.npz",
        f"{bam_dir}/enhancer{outname}.npz", 'expr', skip_assert = skip_assert)

    args.nTasks = min(args.nTasks,len(tasklist))

    print(f"selecting reg elements' DNase  from PCHi-C training data for: {args.file_index}/{args.nTasks}..", end=" ", flush=True)
    pt.distribute_task(task_list = tasklist, 
        nTasks = int(args.nTasks), file_index=int(args.file_index), 
        func=objfunc, num_tasks=len(tasklist),
        dry_run=False)
    print(".", flush=True)

    return 0

def combine_DNase(bamfs_cell, outpath): 
    # bamfs_cell = [f"{prep.tmpBaseDir}/{bamf.split('.')[0]}/promoter_Train_Augment.npz" for bamf in pt.flatten(prep.DnaseCells[cell])]
    # outpath = f"{prep.tmpBaseDir_hic}/{cell}/P-E/promoter_Seq.npz"

    dnase_cell = [np.load(bamf,allow_pickle=True)['expr'] for bamf in bamfs_cell]
    dnase_cell = np.stack(dnase_cell)
    dnase_cell = np.transpose(dnase_cell,axes=(1,0,2))

    np.savez(outpath, expr=dnase_cell)
    return 0


def combine_DNase_Reps(pr_bamfs_cell, enh_bamfs_cell, outpath_pr, outpath_enh):
    """
    combine chromatin openness score from different bam file for one cell type
    """
    # pr_bamfs_cell = [f"{prep.tmpBaseDir}/{bamf.split('.')[0]}/promoter_Train_Augment.npz" for bamf in pt.flatten(prep.DnaseCells[cell])]
    # enh_bamfs_cell = [f"{prep.tmpBaseDir}/{bamf.split('.')[0]}/enhancer_Train_Augment.npz" for bamf in pt.flatten(prep.DnaseCells[cell])]

    # outpath_pr = f"{prep.tmpBaseDir_hic}/{cell}/P-E/promoter_Seq.npz"
    # outpath_enh = f"{prep.tmpBaseDir_hic}/{cell}/P-E/enhancer_Seq.npz"

    tasklist = [(pr_bamfs_cell,outpath_pr),(enh_bamfs_cell,outpath_enh)]

    objfunc = lambda task:combine_DNase(*task)

    args.nTasks = min(args.nTasks,len(tasklist))

    print(f"selecting reg elements' DNase  from PCHi-C training data for: {args.file_index}/{args.nTasks}..", end=" ", flush=True)
    pt.distribute_task(task_list = tasklist, 
        nTasks = int(args.nTasks), file_index=int(args.file_index), 
        func=objfunc, num_tasks=len(tasklist),
        dry_run=False)
    print(".", flush=True)

    return 0




def SelectElements_DNA(cell_trainData, prDNA, enhDNA, all_prList, all_enhList, outdir, out_append = '', skip_assert = True):
    """
    select DNA elements
    """
    # cell_trainData = prep.hictrain_augment(cell)
    # prDNA = prep.promoter['dna-out']
    # enhDNA = prep.enhancer['dna-out']
    # all_prList = pd.read_csv(prep.promoter['bed-path'], delimiter='\t', names=prep.promoter['headers']).name
    # all_enhList = pd.read_csv(prep.enhancer['bed-path'], delimiter='\t', names=prep.enhancer['headers']).name

    reshapeDNA = lambda x:x.reshape(-1,x.shape[-1]//4,4)

    print(f"selecting  reg elements' DNA-sequence from PCHi-C training data for: {args.file_index}/{args.nTasks}..", end=" ", flush=True)

    pt.SelectElements_in_TrainingData(cell_trainData, reshapeDNA(np.load(prDNA, allow_pickle=True)['sequence']),
        reshapeDNA(np.load(enhDNA, allow_pickle=True)['sequence']),
        all_prList, all_enhList, 
        f"{outdir}/promoterDNA_{out_append}.npz",
        f"{outdir}/enhancerDNA_{out_append}.npz", 'sequence', 
        transform = lambda x:x.reshape(-1,x.shape[-1]*x.shape[-2]),
        skip_assert = skip_assert)

    print(".", flush=True)

    return 0


def SelectElements_DNA_parallel(args, cell_trainData, prDNA, enhDNA, all_prList, all_enhList, outdir, tmp_dir = ".",  out_append = '', skip_assert = False, parallel_num = 100, task="combine"):
    """
    select DNA elements parallel
    """

    assert task in ['split','run','combine']
    split_csvs = [f"{tmp_dir}/{out_append}_csv_part_{i}.csv" for i in range(parallel_num)]
    if task == 'split':
        pt.splitCSV(cell_trainData, split_csvs, readArgs= {}, writeArgs= {'index':False})
    elif task == 'run':
        applyfunc = lambda index, augmented_trainData_cell : SelectElements_DNA(augmented_trainData_cell, prDNA, enhDNA, all_prList, all_enhList, 
            outdir, out_append = f"{out_append}_{index}", skip_assert = skip_assert)

        tasklist = list(zip(range(args.nTasks),split_csvs))

        args.nTasks = min(args.nTasks,len(tasklist))

        print(f"{out_append}: applying {applyfunc.__name__} to csv (file_index = {args.file_index}/{len(tasklist)})..", end=" ", flush=True)
        pt.distribute_task(task_list = tasklist, 
            nTasks = int(args.nTasks), file_index=int(args.file_index), 
            func=lambda x:applyfunc(*x), num_tasks=len(tasklist),
            dry_run=False)
        print(".", flush=True)

    else:
        prDNA = np.concatenate([np.load(f"{outdir}/promoterDNA_{cell}_{index}.npz", allow_pickle=True)['sequence'] for index in range(parallel_num)])
        enhDNA = np.concatenate([np.load(f"{outdir}/enhancerDNA_{cell}_{index}.npz", allow_pickle=True)['sequence'] for index in range(parallel_num)])

        label_data = np.array(pd.read_csv(cell_trainData).label)

        np.savez(f"{outdir}/promoterDNA_{cell}.npz", **{'sequence':prDNA, 'label':label_data})
        np.savez(f"{outdir}/enhancerDNA_{cell}.npz", **{'sequence':enhDNA, 'label':label_data})

        # np.load(f"{prep.tmpBaseDir}/promoterDNA_{cell}.npz", allow_pickle=True)

    return 0


def seperate_data(CELL, TYPE='P-E', NUM_SEQ=4):
    if TYPE == 'P-P':
        filename1 = 'promoter1'
        filename2 = 'promoter2'
        RESIZED_LEN = 1000 #promoter
    elif TYPE == 'P-E':
        filename1 = 'enhancer'
        filename2 = 'promoter'
        RESIZED_LEN = 2000 #enhancer
    else:
        print('ERROR')
        sys.exit()

    print(os.path.basename(CELL), TYPE,flush=True)
    print("reading data..",flush=True)
    shape1 = (-1, 1, RESIZED_LEN, NUM_SEQ)
    shape2 = (-1, 1, 1000, NUM_SEQ)
    region1 = np.load(CELL+'/'+TYPE+'/'+filename1+'_Seq.npz')
    region2 = np.load(CELL+'/'+TYPE+'/'+filename2+'_Seq.npz')
    Tlabel = region1['label']
    Tregion1_seq = region1['sequence'].reshape(shape1).transpose(0, 1, 3, 2)
    Tregion2_seq = region2['sequence'].reshape(shape2).transpose(0, 1, 3, 2)


    ## load data: DNase
    region1 = np.load(CELL+'/'+TYPE+'/'+filename1+'_DNase.npz')
    region2 = np.load(CELL+'/'+TYPE+'/'+filename2+'_DNase.npz')
    Tregion1_expr = region1['expr']
    NUM_REP = Tregion1_expr.shape[1]
    shape1 = (-1, 1, NUM_REP, RESIZED_LEN)
    shape2 = (-1, 1, NUM_REP, 1000)

    Tregion1_expr = Tregion1_expr.reshape(shape1)
    Tregion2_expr = region2['expr'].reshape(shape2)

    NUM = Tlabel.shape[0]
    os.makedirs(f'{CELL}/{TYPE}/data',exist_ok=True)
    print("saving data..",flush=True)
    for index in range(NUM):
        np.savez(CELL+'/'+TYPE+'/data/'+f"{index}.npz", enh_seq = Tregion1_seq[index]
                 , pr_seq = Tregion2_seq[index]
                 , enh_dnase = Tregion1_expr[index]
                 , pr_dnase = Tregion2_expr[index], label = Tlabel[index])
    return 0

if __name__=="__main__":
    import preprocessing as prep

    argType = {'file_index':None, 'nTasks':None, 'taskType':None , 'cellType': None, 'num_rep':{'type':int,  'help':'bio repeat samples of DNase for cell type','default':1}} 
    args = pt.process_inputArgs(input_parse=sys.argv[1:], argType=argType)
    # args = pt.process_inputArgs(input_parse=['--file_index','0','--nTasks','52','--taskType','pWin'])

    taskTypes = ['hicMatch', 'hicLabels', 'selectDNase','combineDNaseReps', 'selectDNA', 'sepData']
    # taskTypes = [t+_t for t in ['','p','e'] for _t in _taskTypes]
    assert args.taskType in taskTypes

    gtfFile = prep.gtf
    HiCParts = prep.HiCParts
    codeTmpDir = prep.codeTmpDir
    hicTSV = prep.hicTSV
    enhancer_dna = prep.enhancer['dna-out']
    enhancer_bed = prep.enhancer['bed-path']
    enhancer_bed_bg = prep.enhancer['bg-path']
    promoter_dna = prep.promoter['dna-out']
    promoter_bed = prep.promoter['bed-path']
    promoter_bed_bg = prep.promoter['bg-path']
    th_pos = 5
    th_neg = 0

    numSamples = None

    cell = args.cellType
    # for cell in ['tB','tCD4','nCD4','FoeT','Mon','tCD8']:

    tmp_out = [f"{codeTmpDir}/pchicMatch_{cell}_{i}.pkl" for i in range(len(HiCParts))]
    hic_matched = prep.HiC_Match(cell)
    hic_grpMatch = prep.HiC_GroupMatch(cell)
    hic_TrainPE = prep.HiC_UniqueMatchPE(cell)
    hic_TrainEP = prep.HiC_UniqueMatchEP(cell)
    hicTrainPos = prep.HiC_TrainingPos(cell)
    hicTrain = prep.HiC_Training(cell)

    if args.taskType=="hicMatch":
        HiCMatch(args, HiCParts, tmp_out, promoter_dna, enhancer_dna, cell, th_pos)

    elif args.taskType=="hicLabels":
        # pt.combine(tmp_out, hic_matched)
        # print("STEP 1 COMPLETE",flush=True)
        

        # transcriptsToGene(gtfFile, hic_matched, hic_grpMatch)
        # print("STEP 2 COMPLETE",flush=True)
        

        # hicUniqueMatch(hic_grpMatch, hic_TrainPE, hic_TrainEP)
        # print("STEP 3 COMPLETE",flush=True)
        

        # HiCTrainingPosLabel(hic_TrainPE, hic_TrainEP, prep.promoter['window'], prep.enhancer['window'], traindata_out = hicTrainPos)
        # print("STEP 4 COMPLETE",flush=True)
        
    
        # HiCTrainingNegLabel(hicTSV, cell, th_neg, hicTrainPos, promoter_dna, enhancer_dna, prep.promoter['window'], prep.enhancer['window'], hic_out = hicTrain, numSamples=numSamples)
        # os.makedirs(f"{prep.tmpBaseDir_hic}/{cell}/P-E",exist_ok=True)
        # os.rename(prep.HiC_Training(cell), f"{prep.tmpBaseDir_hic}/{cell}/P-E/pairs.csv")
        # print("STEP 5 COMPLETE",flush=True)
        

        # cmd = f"""python ../DeepTACT-master/DataPrepare.py {prep.tmpBaseDir_hic}/{cell} P-E"""
        # os.system(cmd)
        # print("STEP 6 COMPLETE",flush=True)
        

        train_augment(prep.hictrain(cell), prep.hictrain_augment(cell), 
            prep.enhancer['ext-window']-prep.enhancer['window'], prep.promoter['ext-window']-prep.promoter['window'], 
            prep.enhancer['aug_step'],  prep.promoter['aug_step'], 
            prep.enhancer['window'], prep.promoter['window'])
        print("STEP 7 COMPLETE",flush=True)
        

    elif args.taskType=="selectDNase":
        all_prList = pd.read_csv(prep.promoter['bed-path'], delimiter='\t', names=prep.promoter['headers']).name
        all_enhList = pd.read_csv(prep.enhancer['bed-path'], delimiter='\t', names=prep.enhancer['headers']).name

        bamfiles_cell = pt.flatten(prep.DnaseCells[cell])
        SelectElements_DNase(args, prep.hictrain_augment(cell), bamfiles_cell, prep.dnaseTmpDir, all_prList, all_enhList, outname = "_Train_Augment", skip_assert = False)

    elif args.taskType=="combineDNaseReps":
        pr_bamfs_cell = [f"{prep.tmpBaseDir}/{bamf.split('.')[0]}/promoter_Train_Augment.npz" for bamf in pt.flatten(prep.DnaseCells[cell])]
        enh_bamfs_cell = [f"{prep.tmpBaseDir}/{bamf.split('.')[0]}/enhancer_Train_Augment.npz" for bamf in pt.flatten(prep.DnaseCells[cell])]

        outpath_pr = f"{prep.tmpBaseDir_hic}/{cell}/P-E/promoter_DNase.npz"
        outpath_enh = f"{prep.tmpBaseDir_hic}/{cell}/P-E/enhancer_DNase.npz"

        combine_DNase_Reps(pr_bamfs_cell, enh_bamfs_cell, outpath_pr, outpath_enh)

    elif args.taskType=="selectDNA":
        all_prList = pd.read_csv(prep.promoter['bed-path'], delimiter='\t', names=prep.promoter['headers']).name
        all_enhList = pd.read_csv(prep.enhancer['bed-path'], delimiter='\t', names=prep.enhancer['headers']).name

        # SelectElements_DNA(prep.hictrain_augment(cell), prep.promoter['dna-out'], prep.enhancer['dna-out'], all_prList, all_enhList, prep.tmpBaseDir, out_append = cell,  skip_assert = False)

        parallel_num = 100
        SelectElements_DNA_parallel(args, prep.hictrain_augment(cell), prep.promoter['dna-out'], prep.enhancer['dna-out'], all_prList, all_enhList, 
            prep.tmpBaseDir, tmp_dir = prep.codeTmpDir,  out_append = cell, skip_assert = False, parallel_num = parallel_num, task="combine")

    elif args.taskType=="sepData":
        CELL = f"{prep.tmpBaseDir_hic}/{cell}"
        # NUM_REP = args.num_rep
        seperate_data(CELL, TYPE='P-E', NUM_SEQ=4)







