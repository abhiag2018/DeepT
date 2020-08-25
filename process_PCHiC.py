import sys
import numpy as np
import argparse
import gtfparse
import preptools as pt
import pandas as pd
import random


def HiCMatch(hicTSV, hic_out, promoter_dna, enhancer_dna, cell, th):
    """
    wrapper for matching PCHiC data in tsv file to promoter and enhancer lists
    """
    print('matching PCHiC Data to promoter & enhancers..',end=" ",flush=True)
    pchicDF = pt.concat_PCHiC_PE(hicTSV,promoter_dna,enhancer_dna, selectCell=cell, threshold = th,
                                 outputF=hic_out,
                                 sampleFrac=.01)
    print('.',flush=True)
    return 0

def transcriptsToGene(gtfFile, hicMatched, hic_out):
    """
    convert Promoter Transcript IDs in matches from PCHiC file to gene-symbols
    """
    pchicDF = pd.read_pickle(hicMatched)

    pchicDF['baitPr'] =  pchicDF['baitPr'].apply(np.array)
    pchicDF['baitEnh'] =  pchicDF['baitEnh'].apply(np.array)
    pchicDF['oePr'] =  pchicDF['oePr'].apply(np.array)
    pchicDF['oeEnh'] =  pchicDF['oeEnh'].apply(np.array)

    gtf = gtfparse.read_gtf(gtfFile)
    # group Transcripts by Gene Symbols
    groupTS = lambda elemBaitPr:pt.GroupGeneSymbol(gtf,elemBaitPr[:,1],elemBaitPr[:,0]) if elemBaitPr.size>0 else []
    
    print('converting transcripts IDs in  PCHiC to GeneSymbols..',end=" ",flush=True)
    pchicDF['baitPr'] = pchicDF['baitPr'].apply(groupTS)
    print('..',end=" ",flush=True)
    pchicDF['oePr'] = pchicDF['oePr'].apply(groupTS)
    print('.',flush=True)

    pchicDF.to_pickle(hic_out)
    # pd.read_pickle(hic_out)
    return 0

def hicUniqueMatch(hicGroupMatched, hic_out_PE=None, hic_out_EP=None):
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
        pehic.to_pickle(hic_out_PE)
    if hic_out_EP:
        ephic.to_pickle(hic_out_EP)

    return pehic,ephic


def hicTrainingCSV(hicUniquePE, hicUniqueEP, prwin, enhwin, traindata_out = None):
    """
    select the data from the columns of the training data and put them in appropriate format
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
        train_out['promoter_name'] = pchicDF[prefixPr+'Pr'].apply(lambda x:x[0])

        train_out['label'] = 1
        return pd.DataFrame(train_out)


    dfPE = niceHiCDF(hicUniquePE, prwin, enhwin, prefixPr="bait", prefixEnh="oe")
    dfEP = niceHiCDF(hicUniqueEP, prwin, enhwin, prefixPr="oe", prefixEnh="bait")

    trainDF =  pd.concat([dfPE,dfEP])
    if traindata_out:
        trainDF.to_csv(traindata_out, index=False)
    return trainDF


if __name__=="__main__":
    import preprocessing as prep

    args = pt.process_inputArgs(input_parse=sys.argv[1:])
    # args = pt.process_inputArgs(input_parse=['--file_index','0','--nTasks','52','--taskType','pWin'])

    taskTypes = ["Convert"]
    # taskTypes = [t+_t for t in ['','p','e'] for _t in _taskTypes]
    assert args.taskType in taskTypes

    gtfFile = prep.gtf
    hicTSV = prep.hicTSV
    promoter_dna = prep.promoter['dna-out']
    enhancer_dna = prep.enhancer['dna-out']
    cell = 'MK'
    th = 5

    hic1 = prep.HiC_Match(cell)
    hic2 = prep.HiC_GroupMatch(cell)
    hic_out_PE = prep.HiC_UniqueMatchPE(cell)
    hic_out_EP = prep.HiC_UniqueMatchEP(cell)
    hicTraining = prep.HiC_Training(cell)
    prwin = prep.promoter['window']
    enhwin =  prep.enhancer['window']
    # for cell in ['tB','tCD4','nCD4','FoeT','Mon','tCD8']:

    # HiCMatch(hicTSV, hic1, promoter_dna, enhancer_dna, cell, th)

    # transcriptsToGene(gtfFile, hic1, hic2)

    # hicUniqueMatch(hic2, hic_out_PE, hic_out_EP)

    hicTrainingCSV(hic_out_PE, hic_out_EP, prwin, enhwin, traindata_out = hicTraining)


