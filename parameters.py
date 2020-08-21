import colored_traceback.always
import os, re, glob
import pandas as pd
import preptools as pt

LIMIT = 1

baseDataDir = "/projects/li-lab/agarwa/CUBE/DeepTact/dataset"
codeTmpDir = "/projects/li-lab/agarwa/CUBE/DeepTact/code/tmp_data"
bgWindow = int(1e6)

hg19 = f"{baseDataDir}/hg19.fa"

## promoter parameters
promoter = {'headers': ['chrom', 'txStart', 'txEnd', 'name', 'score', 'strand',
               'cdsStart', 'cdsEnd', '--', 'exonCount', 'exonStarts',
               'exonEnds'],
    'window': 1000,
    'siteWindow': 200,
    'allfield-bed': f"{baseDataDir}/hg19_promoter_allFields.bed"
}

    
## enhancer parameters
enhancer = {'headers': ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
           'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes',
           'blockStarts'],
    'window': 2000,
    'siteWindow': 200,
    'allfield-bed': f"{baseDataDir}/enhancers_fantom5/human_permissive_enhancers_phase_1_and_2.bed"
    }



## DNase input parameters
bamDir = f"{baseDataDir}/Dnase-Seq/cellTypes"
bamfilesInit = [fn for fn in glob.glob(f"{bamDir}/**/*.bam",recursive=True) if re.match("[^.]*\.bam$",fn)][:LIMIT]
intersectOptions="-c -e -f 0.5 -F 0.5"

clearRun = True
reRun = False



## PCHiC Training Data

hicTSV = f"{baseDataDir}/Javierre_ref_18/DATA_S1/PCHiC_peak_matrix_cutoff5.tsv"

# gtf = f"{baseDataDir}/gencode.v19.annotation.gtf"
gtf = f"{baseDataDir}/gencode.v3c.annotation.GRCh37.gtf" #ref release. 09/2009




















