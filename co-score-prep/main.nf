nextflow.preview.dsl=2
// include { promoter_bed; enhancer_bed } from "${params.codeDir}/NN-feature-prep/modules/pr-enh-prep"


if (params.species == "hg" ){
  all_chrom = "$params.all_chrom_hg"
  ch_enhancer_bed_prep = Channel.fromPath("$params.store_dir/enhancer_hg.bed")
  ch_enhancer_bed_prep_bg = Channel.fromPath("$params.store_dir/enhancer_bg_hg.bed")
  ch_promoter_bed_prep = Channel.fromPath("$params.store_dir/promoter_hg.bed")
  ch_promoter_bed_prep_bg = Channel.fromPath("$params.store_dir/promoter_bg_hg.bed")
  promoter_headers = "$params.promoter_headers_hg"
  enhancer_headers = "$params.enhancer_headers"
  insertRGB="True"
  ch_chrLen = Channel.value(params.chromLen_hg)
}
else if (params.species == "mm" ){
  all_chrom = "$params.all_chrom_mm"
  ch_enhancer_bed_prep = Channel.fromPath("$params.store_dir/enhancer_mm.bed")
  ch_enhancer_bed_prep_bg = Channel.fromPath("$params.store_dir/enhancer_bg_mm.bed")
  ch_promoter_bed_prep = Channel.fromPath("$params.store_dir/promoter_mm.bed")
  ch_promoter_bed_prep_bg = Channel.fromPath("$params.store_dir/promoter_bg_mm.bed")
  promoter_headers = "$params.promoter_headers_mm"
  enhancer_headers = "$params.enhancer_headers"
  insertRGB="False"
  ch_chrLen = Channel.value(params.chromLen_mm)
}
else{
    println "species (: $params.species) should be hg or mm; "
    exit 1
}


ch_chrom_val = Channel.fromList(Eval.me(all_chrom))

Channel.fromPath("bam-input.csv", checkIfExists: true)
  .splitCsv(header:true, sep:',')
  .map { row -> [ row.celltype, row.repetition, file(row.bam, checkIfExists: true) ]  }
  .set {ch_input_bam}


// DNase/ATAC preprocessing step 1
// generate .bam.bai index file from DNase/ATAC-seq .bam files
process GEN_BAM_INDEX {
    input:
    tuple val(cell), val(rep), path(bamfile)

    output:
    tuple val(cell), val(rep), file("${bamfile}.bai"), file("${bamfile}")

    script:  
    """
      samtools index $bamfile ${bamfile}.bai
    """
}

// DNase/ATAC preprocessing step 2
// chromatin openness score construction from .bam file:
// at each genomic location count the number of reads intersecting with that location
// .. for each .bam and for each chromosome
process CHROM_OPENN_SCORE {
    input:
    tuple val(cell), val(rep), path(bamfile_index), path(bamfile), val(chrom), val(chrLen)

    output:
    tuple val(cell), val(rep), path("*.npz")


    script:  
    len_val=chrLen["${chrom}"]
    """
    python -c "from preptools import getBamCounts; \
    getBamCounts('${bamfile}', \
        '${chrom}', \
        $len_val, \
        outputf='${cell}.rep${rep}.${chrom}.npz')"
    """
}

// DNase/ATAC preprocessing step 3.1
// generate chromatin openness core profile (COscore) for promoter
process CHROM_OPENN_SCORE_PROFILE_PROMOTER {
    storeDir "${params.save_dir}"

    input:
    tuple val(cell), val(rep), path(npzlist), path(promoter), path(promoter_bg)

    output:
    tuple val(cell), val(rep), path("promoter_COscore.${cell}.rep${rep}.npz.gz")

    script:  
    // npzlist_ = npzlist.toList()
    """
    python -c "import preptools as pt; \
    pt.generate_dnase('$npzlist'.split(), \
        '$promoter', \
        '$promoter_bg', \
        $promoter_headers, \
        $params.bgWindow, \
        'promoter_COscore.${cell}.rep${rep}')"
    gzip promoter_COscore.${cell}.rep${rep}.npz
    """
}

// DNase/ATAC preprocessing step 3.2
// generate chromatin openness core profile (COscore) for enhancer
process CHROM_OPENN_SCORE_PROFILE_ENHANCER {
    storeDir "${params.save_dir}"

    input:
    tuple val(cell), val(rep), path(npzlist), path(enhancer), path(enhancer_bg)

    output:
    tuple val(cell), val(rep), path("enhancer_COscore.${cell}.rep${rep}.npz.gz")

    script:  
    """
    python -c "import preptools as pt; \
    pt.generate_dnase('$npzlist'.split(), \
        '$enhancer', \
        '$enhancer_bg', \
        $enhancer_headers, \
        $params.bgWindow, \
        'enhancer_COscore.${cell}.rep${rep}')"
    gzip enhancer_COscore.${cell}.rep${rep}.npz
    """
}



workflow _prep_co_score{
    ch_bam_index = GEN_BAM_INDEX(ch_input_bam) 
    ch_bam_chrom = ch_bam_index.combine(ch_chrom_val)
    ch_bam_chrom_readcounts = CHROM_OPENN_SCORE(ch_bam_chrom.combine(ch_chrLen))

    ch_bam_readcounts_grouped = ch_bam_chrom_readcounts.groupTuple(by: [0,1])

    emit: ch_bam_readcounts_grouped
}

workflow prep_co_score_enh{
    ch_bam_readcounts_grouped = _prep_co_score()
    // ch_enhancer_bed_prep = enhancer_bed()

    ch_enhancer_COscore = CHROM_OPENN_SCORE_PROFILE_ENHANCER(ch_bam_readcounts_grouped.combine(ch_enhancer_bed_prep.merge(ch_enhancer_bed_prep_bg)))
    emit: ch_enhancer_COscore
}

workflow prep_co_score_pr{
    ch_bam_readcounts_grouped = _prep_co_score()
    // ch_promoter_bed_prep = promoter_bed()
    ch_promoter_COscore = CHROM_OPENN_SCORE_PROFILE_PROMOTER(ch_bam_readcounts_grouped.combine(ch_promoter_bed_prep.merge(ch_promoter_bed_prep_bg)))
    emit: ch_promoter_COscore
}

workflow {
    prep_co_score_enh().view()
    prep_co_score_pr().view()
}
// nextflow  -C DeepTact-input-preprocessing/CO_score.config -C promoter-enhancer-preprocesssing/nextflow.config run DeepTact-input-preprocessing/CO_score_reg_elem.nf -profile slurm -w "/fastscratch/agarwa/nf-tmp/work" -with-timeline -resume 












