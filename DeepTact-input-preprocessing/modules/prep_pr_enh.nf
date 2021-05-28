nextflow.preview.dsl=2

params.dataDir="/projects/li-lab/agarwa/CUBE/DeepTact/dataset"

// regulatory element preprocessing options
params.bgWindow = 1000000 // background length for regulatory elements = 1MB (both promoters and enhancers)
params.augment_length = 1000
params.augment_step = 50

params.promoter_bedfile = "${params.dataDir}/hg19_promoter_allFields.bed"
params.promoter_window = 1000
params.promoter_headers = "['chrom', 'txStart', 'txEnd', 'name', 'score', 'strand', 'cdsStart', 'cdsEnd', '--', 'exonCount', 'exonStarts', 'exonEnds']"

params.enhancer_bedfile = "${params.dataDir}/enhancers_fantom5/human_permissive_enhancers_phase_1_and_2.bed"
params.enhancer_window = 2000
params.enhancer_headers = "['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']"

ch_promoter_bed = Channel.fromPath( params.promoter_bedfile )
ch_enhancer_bed = Channel.fromPath( params.enhancer_bedfile )

// promoter pre-processing
process PREPROCESS_PROMOTER {
    storeDir "${params.store_dir}"

    input:
    path(input)

    output:
    // stdout result
    tuple file("promoter.bed"), file("promoter_bg.bed")


    script:  
    """
    python -c "from preptools import process_promoter_bed; \
    process_promoter_bed('$input', \
        'promoter.bed', \
        $params.promoter_headers, \
        window = $params.promoter_window+$params.augment_length)"

    python -c "from preptools import process_promoter_bed; \
    process_promoter_bed('$input', \
        'promoter_bg.bed', \
        $params.promoter_headers, \
        window = $params.bgWindow)"
    """
    // func( allfield_bed, bg_path, headers, window=bg_window)
}

// enhancer pre-processing
process PREPROCESS_ENHANCER {
    storeDir "${params.store_dir}"

    input:
    path(input)

    output:
    // stdout result
    tuple file('enhancer.bed'), file('enhancer_bg.bed')


    script:  
    """
    python -c "from preptools import process_enhancer_bed; \
    process_enhancer_bed('$input', \
        'enhancer.bed', \
        $params.enhancer_headers, \
        window = $params.enhancer_window+$params.augment_length)"

    python -c "from preptools import process_enhancer_bed; \
    process_enhancer_bed('$input', \
        'enhancer_bg.bed', \
        $params.enhancer_headers, \
        window = $params.bgWindow)"
    """
}

workflow {
    ch_promoter_bed_prep = PREPROCESS_PROMOTER( ch_promoter_bed )
    ch_enhancer_bed_prep = PREPROCESS_ENHANCER( ch_enhancer_bed )
}
