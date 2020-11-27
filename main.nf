#! /projects/li-lab/agarwa/conda_envs/cube/bin/nextflow

ch_input = Channel.fromPath(params.input)
ch_chrom_val = Channel.fromList(params.chromList)
ch_chrLen = Channel.value(params.chromLen_GRCh37v13)

ch_promoter_bed = Channel.fromPath( params.promoter_bedfile )
ch_enhancer_bed = Channel.fromPath( params.enhancer_bedfile )

ch_input
  .splitCsv(header:true, sep:',')
  .map { row -> [ row.celltype, row.repetition, file(row.bam, checkIfExists: true) ]  }
  .set {ch_input_bam}

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
    python -c "from preptools import process_promoter_bed; process_promoter_bed('$input', 'promoter.bed', $params.promoter_headers, window = $params.promoter_window+$params.augment_length)"
    python -c "from preptools import process_promoter_bed; process_promoter_bed('$input', 'promoter_bg.bed', $params.promoter_headers, window = $params.bgWindow)"
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
    python -c "from preptools import process_enhancer_bed; process_enhancer_bed('$input', 'enhancer.bed', $params.enhancer_headers, window = $params.enhancer_window+$params.augment_length)"
    python -c "from preptools import process_enhancer_bed; process_enhancer_bed('$input', 'enhancer_bg.bed', $params.enhancer_headers, window = $params.bgWindow)"
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
    python -c "from preptools import getBamCounts; getBamCounts('${bamfile}', 'chr'+'${chrom}', $len_val, outputf='${cell}.rep${rep}.chr${chrom}.npz')"
    """
}

ch_bam_chrom_readcounts
    .groupTuple(by: [0,1])
    .into{ ch_bam_readcounts1;ch_bam_readcounts2 }

ch_bam_readcounts1
    .combine(ch_promoter_bed_prep)
    .set {  ch_bam_readcounts_promoter }

ch_bam_readcounts2
    .combine(ch_enhancer_bed_prep)
    .set {  ch_bam_readcounts_enhancer }


// DNase/ATAC preprocessing step 3.1
// generate chromatin openness core profile (COscore) for promoter
process CHROM_OPENN_SCORE_PROFILE_PROMOTER {
    publishDir "${params.outdir}/promoters", mode: params.publish_dir_mode
    input:
    tuple val(cell), val(rep), path(npzlist), path(promoter), path(promoter_bg) from ch_bam_readcounts_promoter

    output:
    file("promoter_COscore.npz") into ch_promoter_COscore

    script:  
    // npzlist_ = npzlist.toList()
    """
    python -c "import preptools as pt; pt.generate_dnase('$npzlist'.split(), '$promoter', '$promoter_bg', $params.promoter_headers, $params.bgWindow, 'promoter_COscore')"
    """
}

// DNase/ATAC preprocessing step 3.2
// generate chromatin openness core profile (COscore) for enhancer
process CHROM_OPENN_SCORE_PROFILE_ENHANCER {
    publishDir "${params.outdir}/enhancers", mode: params.publish_dir_mode
    input:
    tuple val(cell), val(rep), path(npzlist), path(enhancer), path(enhancer_bg) from ch_bam_readcounts_enhancer

    output:
    file("enhancer_COscore.npz") into ch_enhancer_COscore

    script:  
    """
    python -c "import preptools as pt; pt.generate_dnase('$npzlist'.split(), '$enhancer', '$enhancer_bg', $params.enhancer_headers, $params.bgWindow, 'enhancer_COscore')"
    """
}



// result.view { it.trim() }

//out_bed_ch
  //.collectFile(name: 'blast_output_combined.txt', storeDir: params.dataDir)