nextflow.preview.dsl=2
include { PREPROCESS_PROMOTER; PREPROCESS_ENHANCER } from './modules/pr-enh-prep'

ch_promoter_bed = Channel.fromPath( params.promoter_bedfile )
ch_enhancer_bed = Channel.fromPath( params.enhancer_bedfile )

ch_chrom_val = Channel.fromList(params.chromList)
ch_chrLen = Channel.value(params.chromLen_GRCh37v13)

Channel.fromPath(params.coScore_data)
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
    tuple val(cell), val(rep), path(bamfile_index), path(bamfile), val(chrom)
    val(chrLen)

    output:
    tuple val(cell), val(rep), path("*.npz")


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

// DNase/ATAC preprocessing step 3.1
// generate chromatin openness core profile (COscore) for promoter
process CHROM_OPENN_SCORE_PROFILE_PROMOTER {
    storeDir "${params.store_dir}"

    input:
    tuple val(cell), val(rep), path(npzlist)
    tuple path(promoter), path(promoter_bg)

    output:
    tuple val(cell), val(rep), path("promoter_COscore.${cell}.rep${rep}.npz.gz")

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
    gzip promoter_COscore.${cell}.rep${rep}.npz
    """
}

// DNase/ATAC preprocessing step 3.2
// generate chromatin openness core profile (COscore) for enhancer
process CHROM_OPENN_SCORE_PROFILE_ENHANCER {
    storeDir "${params.store_dir}"

    input:
    tuple val(cell), val(rep), path(npzlist)
    tuple path(enhancer), path(enhancer_bg)

    output:
    tuple val(cell), val(rep), path("enhancer_COscore.${cell}.rep${rep}.npz.gz")

    script:  
    """
    python -c "import preptools as pt; \
    pt.generate_dnase('$npzlist'.split(), \
        '$enhancer', \
        '$enhancer_bg', \
        $params.enhancer_headers, \
        $params.bgWindow, \
        'enhancer_COscore.${cell}.rep${rep}')"
    gzip enhancer_COscore.${cell}.rep${rep}.npz
    """
}



workflow {
    ch_promoter_bed_prep = PREPROCESS_PROMOTER( ch_promoter_bed )
    ch_enhancer_bed_prep = PREPROCESS_ENHANCER( ch_enhancer_bed )
    ch_bam_index = GEN_BAM_INDEX(ch_input_bam) 
    ch_bam_chrom = ch_bam_index.combine(ch_chrom_val)
    ch_bam_chrom_readcounts = CHROM_OPENN_SCORE(ch_bam_chrom, ch_chrLen)

    ch_bam_readcounts_grouped = ch_bam_chrom_readcounts.groupTuple(by: [0,1])

    ch_promoter_COscore = CHROM_OPENN_SCORE_PROFILE_PROMOTER(ch_bam_readcounts_grouped, ch_promoter_bed_prep)
    ch_enhancer_COscore = CHROM_OPENN_SCORE_PROFILE_ENHANCER(ch_bam_readcounts_grouped, ch_enhancer_bed_prep)
    ch_promoter_COscore.view()
    ch_enhancer_COscore.view()
}

// nextflow  -C DeepTact-input-preprocessing/CO_score.config -C promoter-enhancer-preprocesssing/nextflow.config run DeepTact-input-preprocessing/CO_score_reg_elem.nf -profile slurm -w "/fastscratch/agarwa/nf-tmp/work" -with-timeline -resume 












