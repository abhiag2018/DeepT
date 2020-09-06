import re, os
import numpy as np
import pandas as pd
import preptools as pt
import sys


def SelectElements_in_TrainingData(cell, cell_trainData, all_prList, all_enhList, prDNA, enhDNA, outDir):
    # get training data
    trainData = pd.read_csv(cell_trainData(cell))
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


    # select the DNA Sequence data for elements corresponding to the Training Data
    dnaSeq_pr_train = train_promoters.apply(lambda x: prDNA[all_prList==x][0])
    dnaSeq_enh_train =  train_enh_DF.apply(lambda x: enhDNA[(all_enh_DF['st']<=x['mid']) & (x['mid']<=all_enh_DF['en']) & (x['chr']==all_enh_DF['chr'])][0],axis=1)


    np.savez(f"{outDir}/promoterDNA_{cell}.npz",sequence=np.array(dnaSeq_pr_train.tolist()))
    np.savez(f"{outDir}/enhancerDNA_{cell}.npz",sequence=np.array(dnaSeq_enh_train.tolist()))

    # np.load(f"{outDir}/promoterDNA_{cell}.npz", allow_pickle=True)['sequence']
    return 0



if __name__=="__main__":
    import preprocessing as prep

    argType = {'file_index':None, 'nTasks':None, 'taskType':None , 'cellType':{'type':str, 'help':"cell type in ['tB','tCD4','nCD4','FoeT','Mon','tCD8']", 'default':None}} 
    args = pt.process_inputArgs(input_parse=sys.argv[1:], argType=argType)

    taskTypes = [ 'pDNA' , 'eDNA' , 'selectDNA' ]
    assert args.taskType in taskTypes

    hg19 = prep.hg19
    p = prep.promoter
    e = prep.enhancer
    reRun = prep.reRun

    if args.taskType=="pDNA":
        if reRun or not os.path.exists(p['fa-out']):
            pt.generate_elem_fa(hg19, p['bed-path'], out_fa = p['fa-out'])
        prName,prLoc,prSeq,prOhc = pt.fasta_to_onehot(p['fa-out'],outp=p['dna-out'])

    if args.taskType=="eDNA":
        if reRun or not os.path.exists(e['fa-out']):
            pt.generate_elem_fa(hg19, e['bed-path'], out_fa = e['fa-out'])
        enhName,enhLoc,enhSeq,enhOhc = pt.fasta_to_onehot(e['fa-out'],outp=e['dna-out'])

    if args.taskType=="selectDNA":
        cell = args.cellType
        prDNA = np.load(p['dna-out'], allow_pickle=True)['sequence']
        enhDNA = np.load(e['dna-out'], allow_pickle=True)['sequence']
        all_prList = pd.read_csv(p['bed-path'], delimiter='\t', names=p['headers']).name
        all_enhList = pd.read_csv(e['bed-path'], delimiter='\t', names=e['headers']).name

        SelectElements_in_TrainingData(cell, prep.HiC_Training, all_prList, all_enhList, prDNA, enhDNA, prep.tmpBaseDir)


