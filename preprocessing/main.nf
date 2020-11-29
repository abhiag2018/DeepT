#! /projects/li-lab/agarwa/conda_envs/cube/bin/nextflow

ch_input_fasta = Channel.fromPath(params.species_genome_fasta)
ch_chrom_val = Channel.fromList(params.chromList)
ch_chrLen = Channel.value(params.chromLen_GRCh37v13)

ch_promoter_bed = Channel.fromPath( params.promoter_bedfile )
ch_enhancer_bed = Channel.fromPath( params.enhancer_bedfile )

Channel.fromPath(params.coScore_data)
  .splitCsv(header:true, sep:',')
  .map { row -> [ row.celltype, row.repetition, file(row.bam, checkIfExists: true) ]  }
  .set {ch_input_bam}


ch_cellTypes = Channel.fromList(params.cellTypes)
ch_hic_input = Channel.fromPath(params.hic_input)
    .into { ch_hic_input1; ch_hic_input2}
ch_gtf_input = Channel.fromPath(params.gtf_transcript_to_gene)

// promoter pre-processing
process PREPROCESS_PROMOTER {
    publishDir "${params.outdir}/promoters", mode: params.publish_dir_mode
    input:
    path(input) from ch_promoter_bed

    output:
    // stdout result
    tuple file("promoter.bed"), file("promoter_bg.bed")  into ch_promoter_bed_prep


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
    publishDir "${params.outdir}/enhancers", mode: params.publish_dir_mode
    input:
    path(input) from ch_enhancer_bed

    output:
    // stdout result
    tuple file('enhancer.bed'), file('enhancer_bg.bed') into ch_enhancer_bed_prep


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

// DNase/ATAC preprocessing step 1
// generate .bam.bai index file from DNase/ATAC-seq .bam files
process GEN_BAM_INDEX {
    // publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    tuple val(cell), val(rep), path(bamfile) from ch_input_bam

    output:
    tuple val(cell), val(rep), file("${bamfile}.bai"), file("${bamfile}") into ch_bam_index

    script:  
    """
      samtools index $bamfile ${bamfile}.bai
    """
}

ch_bam_index
    .combine(ch_chrom_val)
    .set {ch_bam_chrom}


// DNase/ATAC preprocessing step 2
// chromatin openness score construction from .bam file:
// at each genomic location count the number of reads intersecting with that location
// .. for each .bam and for each chromosome
process CHROM_OPENN_SCORE {
    input:
    tuple val(cell), val(rep), path(bamfile_index), path(bamfile), val(chrom) from ch_bam_chrom
    val(chrLen) from ch_chrLen

    output:
    // stdout result
    tuple val(cell), val(rep), path("*.npz") into ch_bam_chrom_readcounts


    script:  
    len_val=chrLen["chr${chrom}"]
    """
    python -c "from preptools import getBamCounts; \
    getBamCounts('${bamfile}', \
        'chr'+'${chrom}', \
        $len_val, \
        outputf='${cell}.rep${rep}.chr${chrom}.npz')"
    """
}

ch_promoter_bed_prep
    .into{ ch_promoter_bed_prep1; ch_promoter_bed_prep2; ch_promoter_bed_prep3; ch_promoter_bed_prep4 }

ch_enhancer_bed_prep
    .into{ ch_enhancer_bed_prep1; ch_enhancer_bed_prep2; ch_enhancer_bed_prep3; ch_enhancer_bed_prep4 }

ch_bam_chrom_readcounts
    .groupTuple(by: [0,1])
    .into{ ch_bam_readcounts1;ch_bam_readcounts2 }

ch_bam_readcounts1
    .combine(ch_promoter_bed_prep1)
    .set {  ch_bam_readcounts_promoter }

ch_bam_readcounts2
    .combine(ch_enhancer_bed_prep1)
    .set {  ch_bam_readcounts_enhancer }

ch_input_fasta
    .into {ch_input_fasta1; ch_input_fasta2 }

ch_promoter_bed_prep2
    .combine(ch_input_fasta1)
    .set { ch_fasta_promoter }

ch_enhancer_bed_prep2
    .combine(ch_input_fasta2)
    .set { ch_fasta_enhancer }

// DNase/ATAC preprocessing step 3.1
// generate chromatin openness core profile (COscore) for promoter
process CHROM_OPENN_SCORE_PROFILE_PROMOTER {
    publishDir "${params.outdir}/promoters", mode: params.publish_dir_mode
    input:
    tuple val(cell), val(rep), path(npzlist), path(promoter), path(promoter_bg) from ch_bam_readcounts_promoter

    output:
    tuple val(cell), val(rep), path("promoter_COscore.${cell}.rep${rep}.npz") into ch_promoter_COscore

    script:  
    // npzlist_ = npzlist.toList()
    """
    python -c "import preptools as pt; \
    pt.generate_dnase('$npzlist'.split(), \
        '$promoter', \
        '$promoter_bg', \
        $params.promoter_headers, \
        $params.bgWindow, \
        'promoter_COscore.${cell}.rep${rep}')"
    """
}

// DNase/ATAC preprocessing step 3.2
// generate chromatin openness core profile (COscore) for enhancer
process CHROM_OPENN_SCORE_PROFILE_ENHANCER {
    publishDir "${params.outdir}/enhancers", mode: params.publish_dir_mode
    input:
    tuple val(cell), val(rep), path(npzlist), path(enhancer), path(enhancer_bg) from ch_bam_readcounts_enhancer

    output:
    tuple val(cell), val(rep), path("enhancer_COscore.${cell}.rep${rep}.npz") into ch_enhancer_COscore

    script:  
    """
    python -c "import preptools as pt; \
    pt.generate_dnase('$npzlist'.split(), \
        '$enhancer', \
        '$enhancer_bg', \
        $params.enhancer_headers, \
        $params.bgWindow, \
        'enhancer_COscore.${cell}.rep${rep}')"
    """
}


// promoter genomic sequence computation 
process GENOMIC_SEQUENCE_PROMOTER {
    publishDir "${params.outdir}/promoters", mode: params.publish_dir_mode
    input:
    tuple path(promoter), path(promoter_bg), path(fasta) from ch_fasta_promoter

    output:
    file("promoter_DNAseq.npz") into ch_promoter_DNAseq

    script:  
    """
    python -c "import preptools as pt; \
    pt.generate_elem_fa('$fasta', \
        '$promoter', \
        out_fa = 'promoter.fa' )"

    python -c "import preptools as pt; \
    pt.fasta_to_onehot('promoter.fa', \
        outp='promoter_DNAseq.npz')"
    """
}



// enhancer genomic sequence computation 
process GENOMIC_SEQUENCE_ENHANCER {
    publishDir "${params.outdir}/enhancers", mode: params.publish_dir_mode
    input:
    tuple path(enhancer), path(enhancer_bg), path(fasta) from ch_fasta_enhancer

    output:
    file("enhancer_DNAseq.npz") into ch_enhancer_DNAseq

    script:  
    """
    python -c "import preptools as pt; \
    pt.generate_elem_fa('$fasta', \
        '$enhancer', \
        out_fa = 'enhancer.fa' )"

    python -c "import preptools as pt; \
    pt.fasta_to_onehot('enhancer.fa', \
        outp='enhancer_DNAseq.npz')"
    """
}

// PCHi-C processing : step 1
// split input tsv file for PCHi-C data into params.hic_split_process parts
process SPLIT_HIC {
    input:
    path(hic_input) from ch_hic_input1

    output:
    path 'hic?.tsv' into ch_hic_parts mode flatten

    script:
    """
    python -c "import preptools as pt; \
    pt.splitCSV('$hic_input', \
        [], \
        readArgs = {'delimiter':'\t','dtype':{'baitChr':str,'oeChr':str}}, \
        writeArgs = {'header':True,'sep':'\t'}, prefix='hic', \
        split_num=$params.hic_split_process, \
        suffix='.tsv')"
    """

}

ch_enhancer_DNAseq
    .concat( ch_promoter_DNAseq )
    .collect()
    .into { ch_regElement_DNAseq1; ch_regElement_DNAseq2; ch_regElement_DNAseq3 }

ch_regElement_DNAseq1
    .combine(ch_hic_parts)
    .combine(ch_cellTypes)
    .set{ ch_hic_cell_regElements }

// PCHi-C processing : step 2
// match the elements in the tsv file to elements in the regulatory elements lists
process MATCH_HIC_TO_REG_ELEMENTS {
    input:
    tuple path(enhancer_npz), path(promoter_npz), path(hic_input), val(cellType) from ch_hic_cell_regElements

    output:
    tuple val(cellType), path('*.pkl') into ch_hic_match

    script:  
    """
    python -c "import preptools as pt; \
    pt.concat_PCHiC_PE('$hic_input', \
        '$promoter_npz', \
        '$enhancer_npz', \
        selectCell='$cellType', \
        threshold = $params.pos_threshold, \
        outputF='pchic_match.${cellType}.'+'$hic_input'.split('.')[0]+'.pkl', \
        sampleFrac=None)"
    """ 
}

ch_hic_match
    .groupTuple(by: 0)
    .set{ ch_hic_match_grouped }

// PCHi-C processing : step 3
// combine matched hic .tsv files
process COMBINE_MATCHED_HIC {
    input:
    tuple val(cellType), path(matched_hic_parts) from ch_hic_match_grouped

    output:
    tuple val(cellType), path('*.pkl') into ch_hic_matched mode flatten

    script:  
    """
    python -c "import preptools as pt; \
    pt.combine('$matched_hic_parts'.split(), \
        'hic_matched.${cellType}.pkl')"
    """ 
}

ch_hic_gtf = ch_hic_matched.combine( ch_gtf_input )

// PCHi-C processing : step 4
// combine transcript IDs under the same gene
process COMBINE_PROMOTERS {
    input:
    tuple val(cellType), path(hic_matched), path(gtf_input) from ch_hic_gtf

    output:
    tuple val(cellType), path("hic.unique.${cellType}.csv") into ch_hic_gtfed

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

ch_hic_input2
    .concat(ch_regElement_DNAseq2)
    .collect()
    .combine(ch_hic_gtfed)
    .set{ ch_hic_intermed }

// PCHi-C processing : step 5
// generate negative labels
process GEN_NEG_LABEL {
    publishDir "${params.outdir}/pchic", mode: params.publish_dir_mode

    input:
    tuple path(hic_input), path(enhancer_dnaseq), path(promoter_dnaseq), val(cellType), path(hic_gtfed) from ch_hic_intermed

    output:
    tuple val(cellType), path("hic.pos.neg.${cellType}.csv") into ch_hic_pre_augment

    script:  
    """
    python -c "import process_PCHiC as pchic; \
    pchic.HiCTrainingNegLabel('$hic_input', \
        '$cellType', $params.neg_threshold,  \
        '$hic_gtfed', \
        '$promoter_dnaseq', \
        '$enhancer_dnaseq',  \
        $params.promoter_window, \
        $params.enhancer_window, \
        hic_out = 'hic.pos.neg.${cellType}.csv', \
        numSamples=$params.pos_neg_interac_ratio )"
    """ 
}


// PCHi-C processing : step 6
// augment data
process GEN_AUGMENTED_LABEL {
    publishDir "${params.outdir}/pchic", mode: params.publish_dir_mode

    input:
    tuple val(cellType), path(hic_pos_neg) from ch_hic_pre_augment

    output:
    tuple val(cellType), path("hic.aug.${cellType}.csv") into ch_hic_augment

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
        mult_fac = $params.hic_augment_factor )"
    """ 
}

ch_hic_augment
    .into { ch_hic_augment1; ch_hic_augment2; ch_hic_augment3 }

ch_enhancer_COscore
    .combine( ch_enhancer_bed_prep3 )
    .set{ ch_enhancer_COscore_bedprep }

ch_promoter_COscore
    .combine( ch_promoter_bed_prep3 )
    .set{ ch_promoter_COscore_bedprep }

ch_enhancer_COscore_bedprep
    .join( ch_promoter_COscore_bedprep, by:[0,1] )
    .combine(ch_hic_augment1, by:0)
    .set{ ch_hic_COscore }

// Combining Data : step 1
// combine PCHi-C interactions to get CO score for each element in the list
process COMBINE_PCHIC_CO_SCORE {
    input:
    tuple val(cellType), val(rep), \
        path(enhancer_COscore), path(enhancer_bed), path(enhancer_bg), \
        path(promoter_COscore), path(promoter_bed), path(promoter_bg), \
        path(hic_aug) from ch_hic_COscore

    output:
    tuple val(cellType), val(rep), path("enhancer_hic.${cellType}.rep${rep}.npz"), path("promoter_hic.${cellType}.rep${rep}.npz") into ch_hic_features_COscore

    script:  
    """
    python -c "import preptools as pt; \
    import pandas as pd; \
    import numpy as np; \
    pt.SelectElements_in_TrainingData('$hic_aug', 
        np.load('$promoter_COscore',allow_pickle=True)['expr'],
        np.load('$enhancer_COscore',allow_pickle=True)['expr'],
        pd.read_csv('$promoter_bed', delimiter='\t', names=$params.promoter_headers).name, 
        pd.read_csv('$enhancer_bed', delimiter='\t', names=$params.enhancer_headers).name, 
        'promoter_hic.${cellType}.rep${rep}.npz',
        'enhancer_hic.${cellType}.rep${rep}.npz', 
        'expr', 
        skip_assert = True)"

    """ 
}


ch_hic_features_COscore
    .groupTuple(by: 0)
    .set { ch_hic_features_COscore_reps }

// Combining Data : step 2
// combine PCHi-C interactions to get CO score for each element in the list
process COMBINE_CO_SCORE_REPS {
    publishDir "${params.outdir}/pchic", mode: params.publish_dir_mode

    input:
    tuple val(cellType), val(reps), val(enhancer_COscore_reps), val(promoter_COscore_reps) from ch_hic_features_COscore_reps

    output:
    tuple val(cellType), path("enhancer_hic_COscore.${cellType}.npz"), path("promoter_hic_COscore.${cellType}.npz") into ch_hic_COscore_features_out

    script:  
    """
    python -c "import process_PCHiC as pchic; \
    pchic.combine_DNase('$enhancer_COscore_reps'[1:-1].split(', '), 'enhancer_hic_COscore.${cellType}.npz');
    pchic.combine_DNase('$promoter_COscore_reps'[1:-1].split(', '), 'promoter_hic_COscore.${cellType}.npz')"
    """ 
}

// tuple val(cellType), path(hic_aug) from ch_hic_augment2
// tuple path(enhancer_DNAseq), path(promoter_DNAseq) from ch_regElement_DNAseq3
// tuple path(promoter_bed), path(promoter_bg)  from ch_promoter_bed_prep4
// tuple path(enhancer_bed), path(enhancer_bg)  from ch_enhancer_bed_prep4

ch_hic_augment2
    .combine(ch_regElement_DNAseq3)
    .combine(ch_enhancer_bed_prep4)
    .combine(ch_promoter_bed_prep4)
    .set{ ch_hic_DNA_seq }

// Combining Data : step 3
// combine PCHi-C interactions to get DNA sequence for each element in the list
process COMBINE_PCHIC_DNA_SEQ {
    publishDir "${params.outdir}/pchic", mode: params.publish_dir_mode

    input:
    tuple val(cellType), path(hic_aug), path(enhancer_DNAseq), path(promoter_DNAseq), path(enhancer_bed), path(enhancer_bg), path(promoter_bed), path(promoter_bg)  from ch_hic_DNA_seq

    output:
    tuple val(cellType), path("enhancer_hic_DNA_seq.${cellType}.npz"), path("promoter_hic_DNA_seq.${cellType}.npz") into ch_hic_DNA_seq_features_out

    script:  
    """
    python -c "import preptools as pt; \
    import numpy as np; \
    import pandas as pd; \
    reshapeDNA = lambda x:x.reshape(-1,x.shape[-1]//4,4) ;\
    pt.SelectElements_in_TrainingData('$hic_aug', 
            reshapeDNA(np.load('$promoter_DNAseq', allow_pickle=True)['sequence']),
            reshapeDNA(np.load('$enhancer_DNAseq', allow_pickle=True)['sequence']),
            pd.read_csv('$promoter_bed', delimiter='\t', names=$params.promoter_headers).name, 
            pd.read_csv('$enhancer_bed', delimiter='\t', names=$params.enhancer_headers).name, 
            'promoter_hic_DNA_seq.${cellType}.npz',
            'enhancer_hic_DNA_seq.${cellType}.npz', 
            'sequence', 
            transform = lambda x:x.reshape(-1,x.shape[-1]*x.shape[-2]),
            skip_assert = True)"
    """ 
}

ch_hic_COscore_features_out
    .join( ch_hic_DNA_seq_features_out, by:0 )
    .join( ch_hic_augment3, by:0 )
    .set{ ch_hic_features_out }
    // .view()

// Combining Data : step 4
// separate the data in ch_hic_DNA_seq_features_out, ch_hic_COscore_features_out into data points
process SEPARATE_DATA {
    publishDir "${params.outdir}/pchic", mode: params.publish_dir_mode

    input:
    tuple val(cellType), path(enhancer_hic_CO), path(promoter_hic_CO), path(enhancer_hic_DNAseq), path(promoter_hic_DNAseq), path(hic_aug) from ch_hic_features_out

    output:
    tuple val(cellType), path("features.${cellType}.tar.gz") into ch_hic_features_tar_out

    script:  
    """
    python -c "import preptools as pt; \
    pt.seperate_data('$enhancer_hic_DNAseq', \
        '$promoter_hic_DNAseq', \
        '$enhancer_hic_CO', \
        '$promoter_hic_CO', \
        $params.enhancer_window, \
        $params.promoter_window, \
        '$hic_aug', \
        'data', NUM_SEQ=4)"
    tar -zcf features.${cellType}.tar.gz data
    """ 
}



// result.view { it.trim() }

//out_bed_ch
//.collectFile(name: 'blast_output_combined.txt', storeDir: params.dataDir)

















