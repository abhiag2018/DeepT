import sys, itertools
import numpy as np
import argparse
import gtfparse
import preptools as pt
import pandas as pd
import matplotlib.pyplot as plt
import random

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

def SelectElements_in_TrainingData(cell_trainData, bam_dir, all_prList, all_enhList):

    trainData = pd.read_csv(cell_trainData)
    train_promoters = trainData['promoter_name'].apply(lambda x:x.split(":")[0])
    train_enhancers = trainData['enhancer_name']

    tfunc = lambda x: pd.Series([x.split(':')[0]]+list(map(int, x.split(':')[1].split('-'))))
    all_enh_DF = all_enhList.apply(tfunc).rename(columns={0:'chr',1:'st',2:'en'})


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


    # select the DNase data for promoter corresponding to the Training Data
    dnase_promoter = np.load(f"{bam_dir}/promoter.npz",allow_pickle=True)['expr']
    dnase_enhancer = np.load(f"{bam_dir}/enhancer.npz",allow_pickle=True)['expr']

    dnase_pr_train = train_promoters.apply(lambda x: dnase_promoter[all_prList==x][0])

    dnase_enh_train =  train_enh_DF.apply(lambda x: dnase_enhancer[(all_enh_DF['st']<=x['mid']) & (x['mid']<=all_enh_DF['en']) & (x['chr']==all_enh_DF['chr'])][0],axis=1)

    np.savez(f"{bam_dir}/promoterTrain.npz",expr=np.array(dnase_pr_train.tolist()))
    np.savez(f"{bam_dir}/enhancerTrain.npz",expr=np.array(dnase_enh_train.tolist()))

    # np.load(f"{dirDNase}/enhancerTrain.npz", allow_pickle=True)['expr']
    return 0

def SelectElements(args, cell_trainData, bamfiles_cell, dnaseTmpDir, all_prList, all_enhList):
    # cell_trainData = prep.HiC_Training(cell)
    # bamfiles_cell = list(itertools.chain.from_iterable(prep.DnaseCells[cell]))
    # dnaseTmpDir = prep.dnaseTmpDir
    # all_prList = pd.read_csv(prep.promoter['bed-path'], delimiter='\t', names=prep.promoter['headers']).name
    # all_enhList = pd.read_csv(prep.enhancer['bed-path'], delimiter='\t', names=prep.enhancer['headers']).name

    tasklist = [dnaseTmpDir(bamf) for bamf in bamfiles_cell]
    objfunc = lambda t:SelectElements_in_TrainingData(cell_trainData, t, all_prList, all_enhList)

    args.nTasks = min(args.nTasks,len(tasklist))

    print(f"selecting promoter & enhancers from PCHi-C training data for: {args.file_index}/{args.nTasks}..", end=" ", flush=True)
    pt.distribute_task(task_list = tasklist, 
        nTasks = int(args.nTasks), file_index=int(args.file_index), 
        func=objfunc, num_tasks=len(tasklist),
        dry_run=False)
    print(".", flush=True)

    return 0


if __name__=="__main__":
    import preprocessing as prep

    argType = {'file_index':None, 'nTasks':None, 'taskType':None , 'cellType': None} 
    args = pt.process_inputArgs(input_parse=sys.argv[1:], argType=argType)
    # args = pt.process_inputArgs(input_parse=['--file_index','0','--nTasks','52','--taskType','pWin'])

    taskTypes = ["hicMatch", "hicLabels"]
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

    if args.taskType=="hicLabels":
        # pt.combine(tmp_out, hic_matched)

        # transcriptsToGene(gtfFile, hic_matched, hic_grpMatch)

        # hicUniqueMatch(hic_grpMatch, hic_TrainPE, hic_TrainEP)

        # HiCTrainingPosLabel(hic_TrainPE, hic_TrainEP, prep.promoter['window'], prep.enhancer['window'], traindata_out = hicTrainPos)
    
        # HiCTrainingNegLabel(hicTSV, cell, th_neg, hicTrainPos, promoter_dna, enhancer_dna, prep.promoter['window'], prep.enhancer['window'], hic_out = hicTrain, numSamples=numSamples)

        bamfiles_cell = list(itertools.chain.from_iterable(prep.DnaseCells[cell]))
        all_prList = pd.read_csv(prep.promoter['bed-path'], delimiter='\t', names=prep.promoter['headers']).name
        all_enhList = pd.read_csv(prep.enhancer['bed-path'], delimiter='\t', names=prep.enhancer['headers']).name
        SelectElements(args, prep.HiC_Training(cell), bamfiles_cell, prep.dnaseTmpDir, all_prList, all_enhList)










