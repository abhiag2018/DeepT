import sys
import numpy as np
import argparse
import mygene
import preptools as pt

def transcriptsToGene(hicTSV, hic_out, promoter_dna, enhancer_dna, cell, th):
    """
    convert Transcript IDs in matches from PCHiC file to gene-symbols
    """
    print('constructing PCHiC DataFrame..',end=" ",flush=True)
    pchicDF = pt.concat_PCHiC_PE(hicTSV,promoter_dna,enhancer_dna, selectCell=cell, threshold = th,
                                 outputF=None,
                                 sampleFrac=.001)
    print('.',flush=True)

    pchicDF['baitPr'] =  pchicDF['baitPr'].apply(np.array)
    pchicDF['baitEnh'] =  pchicDF['baitEnh'].apply(np.array)

    mg = mygene.MyGeneInfo()
    # group Transcripts by Gene Symbols
    groupTS = lambda elemBaitPr:pt.GroupGeneSymbol(mg,elemBaitPr[:,1],elemBaitPr[:,0]) if elemBaitPr.size>0 else []
    
    print('converting transcripts IDs in  PCHiC to GeneSymbols..',end=" ",flush=True)
    pchicDF['baitPr'] = pchicDF['baitPr'].apply(groupTS)
    print('.',flush=True)

    pchicDF['PrCount'] = pchicDF['baitPr'].apply(len)
    pchicDF['EnhCount'] = pchicDF['baitEnh'].apply(len)
    pchicDF.to_pickle(hic_out)
    # pd.read_pickle(hic_out)
    return pchicDF

if __name__=="__main__":
    import preprocessing as prep

    args = pt.process_inputArgs(input_parse=sys.argv[1:])
    # args = pt.process_inputArgs(input_parse=['--file_index','0','--nTasks','52','--taskType','pWin'])

    taskTypes = ["Convert"]
    # taskTypes = [t+_t for t in ['','p','e'] for _t in _taskTypes]
    assert args.taskType in taskTypes

    hicTSV = prep.hicTSV
    promoter_dna = prep.promoter['dna-out']
    enhancer_dna = prep.enhancer['dna-out']
    cell = 'MK'
    th = 5

    # for cell in ['tB','tCD4','nCD4','FoeT','Mon','tCD8']:
    hic_out = prep.HiCIntersections(cell)
    transcriptsToGene(hicTSV, hic_out, promoter_dna, enhancer_dna, cell, th)
