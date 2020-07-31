import colored_traceback.always
import os, glob
import pandas as pd
import preptools as pt

baseDataDir = "/projects/li-lab/agarwa/CUBE/DeepTact/dataset"
bgWindow = int(1e6)

hg19 = "hg19.fa"

## promoter parameters
promoter = {'headers': ['chrom', 'txStart', 'txEnd', 'name', 'score', 'strand',
               'cdsStart', 'cdsEnd', '--', 'exonCount', 'exonStarts',
               'exonEnds'],
    'window': 1000,
    'siteWindow': 200,
    'allfield-bed': "hg19_promoter_allFields.bed"
}

    
## enhancer parameters
enhancer = {'headers': ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
           'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes',
           'blockStarts'],
    'window': 2000,
    'siteWindow': 200,
    'allfield-bed': "enhancers_fantom5/human_permissive_enhancers_phase_1_and_2.bed"
    }



## DNase input parameters
bamDir = f"{baseDataDir}/Dnase-Seq/cellTypes"
bamfilesInit = glob.glob(f"{bamDir}/**/*.bam",recursive=True)
intersectOptions="-c -e -f 0.5 -F 0.5"

clearRun = True
reRun = False




