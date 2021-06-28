nextflow.preview.dsl=2
include { prep_gen_seq } from "${params.codeDir}/NN-feature-prep/modules/genome-seq-prep"

ch_cellTypes = Channel.fromList(params.cellTypes)
ch_hic_input = Channel.fromPath( params.dev ? params.hic_dev : params.hic_input)
ch_gtf_input = Channel.fromPath(params.gtf_transcript_to_gene)

hic_split_num = params.dev ? 2 : params.hic_split_process

// PCHi-C processing : step 1
// split input tsv file for PCHi-C data into params.hic_split_process parts
process SPLIT_HIC {
    input:
    path(hic_input)

    output:
    path 'hic*.tsv'

    script:
    """
    python -c "import preptools as pt; \
    pt.splitCSV('$hic_input', \
        [], \
        readArgs = {'delimiter':'\t','dtype':{'baitChr':str,'oeChr':str}}, \
        writeArgs = {'header':True,'sep':'\t'}, prefix='hic', \
        split_num=$hic_split_num, \
        suffix='.tsv')"
    """

}


// PCHi-C processing : step 2
// match the elements in the tsv file to elements in the regulatory elements lists
process MATCH_HIC_TO_REG_ELEMENTS {
    input:
    tuple path(enhancer_npz), path(promoter_npz), path(hic_input), val(cellType)

    output:
    tuple val(cellType), path('*.pkl')

    script:  
    """
    gzip -df $enhancer_npz $promoter_npz
    python -c "import preptools as pt; \
    import os; \
    pt.concat_PCHiC_PE('$hic_input', \
        '$promoter_npz'[:-3], \
        '$enhancer_npz'[:-3], \
        selectCell='$cellType', \
        threshold = $params.pos_threshold, \
        outputF='pchic_match.${cellType}.'+'$hic_input'.split('.')[0]+'.pkl', \
        sampleFrac=None)
    os.system('rm $promoter_npz'[:-3])
    os.system('rm $enhancer_npz'[:-3])"
    """ 
}

// PCHi-C processing : step 3
// combine matched hic .tsv files
process COMBINE_MATCHED_HIC {
    input:
    tuple val(cellType), path(matched_hic_parts)

    output:
    tuple val(cellType), path('*.pkl')

    script:  
    """
    python -c "import preptools as pt; \
    pt.combine('$matched_hic_parts'.split(), \
        'hic_matched.${cellType}.pkl')"
    """ 
}

// PCHi-C processing : step 4
// combine transcript IDs under the same gene
process COMBINE_PROMOTERS {
    input:
    tuple val(cellType), path(hic_matched), path(gtf_input)

    output:
    tuple val(cellType), path("hic.unique.${cellType}.csv")

    script:  
    """
    python -c "import process_PCHiC as pchic; \
    pchic.transcriptsToGene('$gtf_input', \
        '$hic_matched', \
        'hic.${cellType}.csv')"

    python -c "import process_PCHiC as pchic; \
    pchic.hicUniqueMatch('hic.${cellType}.csv', \
        hic_out_PE= 'hic.${cellType}.PE.csv', \
        hic_out_EP= 'hic.${cellType}.EP.csv')"

    python -c "import process_PCHiC as pchic; \
    pchic.HiCTrainingPosLabel('hic.${cellType}.PE.csv', \
        'hic.${cellType}.EP.csv', \
        $params.promoter_window, \
        $params.promoter_window, \
        traindata_out = 'hic.unique.${cellType}.csv')"
    """ 
}

pos_neg_interac_ratio = params.dev ? 1 : params.pos_neg_interac_ratio

// PCHi-C processing : step 5
// generate negative labels
process GEN_NEG_LABEL {
    input:
    tuple path(hic_input), path(enhancer_dnaseq_gz), path(promoter_dnaseq_gz), val(cellType), path(hic_gtfed)


    output:
    tuple val(cellType), path("hic.pos.neg.${cellType}.csv")

    script:  
    enhancer_dnaseq = enhancer_dnaseq_gz.baseName
    promoter_dnaseq = promoter_dnaseq_gz.baseName
    """
    gzip -df $enhancer_dnaseq_gz $promoter_dnaseq_gz
    python -c "import process_PCHiC as pchic; \
    pchic.HiCTrainingNegLabel('$hic_input', \
        '$cellType', $params.neg_threshold,  \
        '$hic_gtfed', \
        '$promoter_dnaseq', \
        '$enhancer_dnaseq',  \
        $params.promoter_window, \
        $params.enhancer_window, \
        hic_out = 'hic.pos.neg.${cellType}.csv', \
        numSamples=$pos_neg_interac_ratio )"
    rm $enhancer_dnaseq $promoter_dnaseq
    """ 
}


hic_augment_factor = params.dev ? 2 : params.hic_augment_factor 

// PCHi-C processing : step 6
// augment data
process GEN_AUGMENTED_LABEL {
    storeDir "${params.store_dir}"

    input:
    tuple val(cellType), path(hic_pos_neg)

    output:
    tuple val(cellType), path("hic.aug.${cellType}.csv")

    script:  
    """
    echo $cellType, $hic_pos_neg

    python -c "import process_PCHiC as pchic; \
    pchic.train_augment('$hic_pos_neg', \
        'hic.aug.${cellType}.csv', \
        $params.augment_length, \
        $params.augment_length, \
        $params.augment_step, \
        $params.augment_step, \
        $params.enhancer_window, \
        $params.promoter_window, \
        mult_fac = $hic_augment_factor )"
    """ 
}

workflow prep_pchic{
    ch_gen_seq = prep_gen_seq()
    ch_hic_cell_regElements = SPLIT_HIC(ch_hic_input).flatten().combine(ch_cellTypes)

    ch_hic_match = MATCH_HIC_TO_REG_ELEMENTS(ch_gen_seq.combine(ch_hic_cell_regElements))
    ch_hic_matched = COMBINE_MATCHED_HIC(ch_hic_match.groupTuple(by: 0)).transpose()
    ch_hic_gtfed = COMBINE_PROMOTERS(ch_hic_matched.combine(ch_gtf_input))
    ch_hic_augment = GEN_NEG_LABEL(ch_hic_input.combine(ch_gen_seq).combine(ch_hic_gtfed)) | GEN_AUGMENTED_LABEL
    emit:
    ch_hic_augment
}

workflow{
    prep_pchic().view()
}

