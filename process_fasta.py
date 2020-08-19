import re, os
import preptools as pt
import sys

if __name__=="__main__":
    import preprocessing as prep
    args = pt.process_inputArgs(input_parse=sys.argv[1:])

    _taskTypes = ["DNA"]
    taskTypes = [t+_t for t in ['','p','e'] for _t in _taskTypes]
    assert args.taskType in taskTypes

    hg19 = prep.hg19
    p = prep.promoter
    e = prep.enhancer
    reRun = prep.reRun

    if args.taskType=="pDNA" or args.taskType=="DNA":
        if reRun or not os.path.exists(p['fa-out']):
            pt.generate_elem_fa(hg19, p['bed-path'], out_fa = p['fa-out'])
        prName,prLoc,prSeq,prOhc = pt.fasta_to_onehot(p['fa-out'],outp=p['dna-out'])

    if args.taskType=="eDNA" or args.taskType=="DNA":
        if reRun or not os.path.exists(e['fa-out']):
            pt.generate_elem_fa(hg19, e['bed-path'], out_fa = e['fa-out'])
        enhName,enhLoc,enhSeq,enhOhc = pt.fasta_to_onehot(e['fa-out'],outp=e['dna-out'])



