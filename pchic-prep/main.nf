nextflow.preview.dsl=2
// include { prep_gen_seq } from "${params.codeDir}/NN-feature-prep/modules/genome-seq-prep"

if (params.species == "hg" ){
    ch_cellTypes = Channel.fromList(params.cellTypes_hg)
    ch_hic_input = Channel.fromPath( params.dev ? params.hic_dev_hg : params.hic_input_hg)
    ch_gtf_input = Channel.fromPath(params.gtf_transcript_to_gene_hg)
    // ch_gen_seq = Channel.fromPath("$params.store_dir/enhancer_DNAseq.hg19.npz.gz").combine(Channel.fromPath("$params.store_dir/promoter_DNAseq.hg19.npz.gz"))
    all_chrom = params.all_chrom_hg
}
else if (params.species == "mm" ){
    ch_cellTypes = Channel.fromList(params.cellTypes_mm)
    ch_hic_input = Channel.fromPath( params.dev ? params.hic_dev_mm : params.hic_input_mm)
    ch_gtf_input = Channel.fromPath(params.gtf_transcript_to_gene_mm)
    // ch_gen_seq = Channel.fromPath("$params.store_dir/enhancer_DNAseq.mm9.npz").combine(Channel.fromPath("$params.store_dir/promoter_DNAseq.mm9.npz"))
    all_chrom = params.all_chrom_mm
}
else{
    println "species (: $params.species) should be hg or mm; "
    exit 1
}

Channel.fromPath("dnaseq.csv")
  .splitCsv(header:true, sep:',')
  .map { row -> [ file("$params.store_dir/$row.enhancer", checkIfExists: true), file("$params.store_dir/$row.promoter", checkIfExists: true) ]  }
  .set {ch_gen_seq}


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
// 6 min
process MATCH_HIC_TO_REG_ELEMENTS {
    input:
    tuple path(enhancer_dnaseq_gz), path(promoter_dnaseq_gz), path(hic_input), val(cellType)

    output:
    tuple val(cellType), path('*.pkl')

    script:  
    ext = enhancer_dnaseq_gz.extension
    if (ext=="gz") {
        enhancer_npz = enhancer_dnaseq_gz.baseName
        promoter_npz = promoter_dnaseq_gz.baseName
    }
    else{
        enhancer_npz = enhancer_dnaseq_gz
        promoter_npz = promoter_dnaseq_gz
    }
    """
    if [ "$ext" == ".gz" ]; then
        gzip -df $enhancer_dnaseq_gz
        gzip -df $promoter_dnaseq_gz
    fi
    python -c "import preptools as pt; \
    import os; \
    pt.concat_PCHiC_PE('$hic_input', \
        '$promoter_npz', \
        '$enhancer_npz', \
        selectCell='$cellType', \
        threshold = $params.pos_threshold, \
        outputF='pchic_match.${cellType}.'+'$hic_input'.split('.')[0]+'.pkl', \
        sampleFrac=None)
    os.system('rm $promoter_npz')
    os.system('rm $enhancer_npz')"
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
    ext = enhancer_dnaseq_gz.extension
    if (ext=="gz") {
        enhancer_dnaseq = enhancer_dnaseq_gz.baseName
        promoter_dnaseq = promoter_dnaseq_gz.baseName
    }
    else{
        enhancer_dnaseq = enhancer_dnaseq_gz
        promoter_dnaseq = promoter_dnaseq_gz
    }
    """
    if [ "$ext" == ".gz" ]; then
        gzip -df $enhancer_dnaseq_gz
        gzip -df $promoter_dnaseq_gz
    fi
    python -c "import process_PCHiC as pchic; \
    pchic.HiCTrainingNegLabel('$hic_input', \
        '$cellType', $params.neg_threshold,  \
        '$hic_gtfed', \
        '$promoter_dnaseq', \
        '$enhancer_dnaseq',  \
        $params.promoter_window, \
        $params.enhancer_window, \
        hic_out = 'hic.pos.neg.${cellType}.csv', \
        numSamples=$pos_neg_interac_ratio, \
        all_chrom = [x[3:] for x in $all_chrom] )"
    rm $enhancer_dnaseq $promoter_dnaseq
    """ 
}

// PCHi-C processing : step 6
// split data into test, train and validation 
process SPLIT_TEST {
    input:
    tuple val(cellType), path(hic_pos_neg)

    output:
    tuple val(cellType), path("hic.pos.neg.${cellType}.train.csv"), path("hic.pos.neg.${cellType}.test.csv"), path("hic.pos.neg.${cellType}.val.csv")

    script:
    """
    python -c "import preptools as pt; \
    pt.splitCSVparts('$hic_pos_neg', \
        [], \
        readArgs = {}, \
        writeArgs = {'header':True,'index':False}, \
        prefix='hic.pos.neg.${cellType}.', \
        split_num=[$params.trainFrac, $params.valFrac, $params.testFrac], \
        suffix='.csv')"
        
    mv hic.pos.neg.${cellType}.0.csv hic.pos.neg.${cellType}.train.csv 
    mv hic.pos.neg.${cellType}.1.csv hic.pos.neg.${cellType}.val.csv 
    mv hic.pos.neg.${cellType}.2.csv hic.pos.neg.${cellType}.test.csv 
    """
}

// PCHi-C processing : step 7
// augment data
process GEN_AUGMENTED_LABEL {
    storeDir "${params.save_dir}"

    input:
    tuple val(cellType), path(hic_pos_neg), val(aug_len), val(hic_aug_factor)

    output:
    path("${hic_pos_neg.baseName}.aug.csv")

    script:  
    """
    #!/usr/bin/env python
    import process_PCHiC as pchic
    pchic.train_augment('$hic_pos_neg',
        '${hic_pos_neg.baseName}.aug.csv', \
        $aug_len, \
        $aug_len, \
        $params.augment_step, \
        $params.augment_step, \
        $params.enhancer_window, \
        $params.promoter_window, \
        mult_fac = $hic_aug_factor )
    """ 
}

workflow prep_pchic{
    // ch_gen_seq = prep_gen_seq()
    ch_hic_cell_regElements = SPLIT_HIC(ch_hic_input).flatten().combine(ch_cellTypes)

    ch_hic_match = MATCH_HIC_TO_REG_ELEMENTS(ch_gen_seq.combine(ch_hic_cell_regElements))
    ch_hic_matched = COMBINE_MATCHED_HIC(ch_hic_match.groupTuple(by: 0)).transpose()
    ch_hic_gtfed = COMBINE_PROMOTERS(ch_hic_matched.combine(ch_gtf_input))
    ch_hic_split = GEN_NEG_LABEL(ch_hic_input.combine(ch_gen_seq).combine(ch_hic_gtfed)) | SPLIT_TEST
    ch_hic_split_filt = ch_hic_split
        .map{it -> [it[0],[it[1],it[2],it[3]]]}
        .transpose()

    // ch_hic_split_valtest = ch_hic_split_filt.filter{ it[1].getBaseName() =~ /^((?!train$).)*$/ }  // all splits that are not train
    ch_hic_split_test = ch_hic_split_filt.filter{ it[1].getBaseName() =~ /test$/ }
    ch_hic_split_val = ch_hic_split_filt.filter{ it[1].getBaseName() =~ /val$/ }
    ch_hic_split_train = ch_hic_split_filt.filter{ it[1].getBaseName() =~ /train$/ }

    ch_hic_split_test = ch_hic_split_test.combine(Channel.value([1,1]))
    ch_hic_split_val = ch_hic_split_val.combine(Channel.value([1,1]))
    ch_hic_split_train = ch_hic_split_train.combine(Channel.value([params.augment_length,params.hic_augment_factor]))
    // ch_hic_split_valtest_ = ch_hic_split_valtest.map(it -> it[1]).collectFile(storeDir:"${params.store_dir}")

    ch_hic_augment = GEN_AUGMENTED_LABEL(ch_hic_split_test.mix(ch_hic_split_val).mix(ch_hic_split_train))

    emit:
    ch_hic_augment
}

workflow{
    prep_pchic().view()

    // ch_hic_input = Channel.value("test").merge(Channel.fromPath("test1.csv"))

    // ch_hic_split = SPLIT_TEST(ch_hic_input)
    //     .map{it -> [it[0],[it[1],it[2],it[3]]]}.transpose()

    // ch_hic_split_valtest = ch_hic_split.filter{ it[1].getBaseName() =~ /^((?!train$).)*$/ } // all splits that are not train
    // ch_hic_split_train = ch_hic_split.filter{ it[1].getBaseName() =~ /train$/ }

    // ch_hic_split_valtest_ = ch_hic_split_valtest.map(it -> it[1]).collectFile(storeDir:"${params.store_dir}")

    // ch_hic_split_valtest_.mix(ch_hic_split_train).view()

}

