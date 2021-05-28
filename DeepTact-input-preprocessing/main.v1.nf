nextflow.preview.dsl=2
include { PREPROCESS_PROMOTER; PREPROCESS_ENHANCER } from './modules/prep_pr_enh'

ch_promoter_bed = Channel.fromPath( params.promoter_bedfile )
ch_enhancer_bed = Channel.fromPath( params.enhancer_bedfile )

workflow {
    ch_promoter_bed_prep = PREPROCESS_PROMOTER( ch_promoter_bed )
    ch_enhancer_bed_prep = PREPROCESS_ENHANCER( ch_enhancer_bed )
}
