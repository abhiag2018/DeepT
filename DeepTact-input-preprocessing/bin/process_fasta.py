import re, os
import numpy as np
import pandas as pd
import preptools as pt
import sys


if __name__=="__main__":
    import preprocessing as prep

    argType = {'file_index':None, 'nTasks':None, 'taskType':None , 'cellType':{'type':str, 'help':"cell type in ['tB','tCD4','nCD4','FoeT','Mon','tCD8']", 'default':None}} 
    args = pt.process_inputArgs(input_parse=sys.argv[1:], argType=argType)

    taskTypes = [ 'pDNA' , 'eDNA' ]
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