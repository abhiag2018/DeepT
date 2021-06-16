nextflow.preview.dsl=2
// include { promoter_bed; enhancer_bed } from "$projectDir/modules/pr-enh-prep" // step 1
// include { prep_co_score_enh; prep_co_score_pr } from "$projectDir/modules/co-score-prep" // step 2
// include { prep_gen_seq } from "$projectDir/modules/genome-seq-prep" // step 3
// include { prep_pchic } from "$projectDir/modules/pchic-prep" // step 4

include {SPLIT_HIC_AUG; COMBINE_PCHIC_CO_SCORE; COMBINE_PCHIC_CO_SCORE_ENH; COMBINE_PCHIC_CO_SCORE_PR; \
    COMBINE_CO_SCORE_REPS_ENH; COMBINE_CO_SCORE_REPS_PR; \
    COMBINE_PCHIC_DNA_SEQ; COMBINE_PCHIC_OUT_ENHANCER; COMBINE_PCHIC_OUT_PROMOTER; \
    SAVE_ENH_DNA_SEQ; SAVE_PR_DNA_SEQ; SAVE_ENH_CO_SCORE; SAVE_PR_CO_SCORE; SEPARATE_DATA; CONVERT_TAR_XZ; COMBINE_DATA_TAR} from "$projectDir/modules/combine" 

ch_enhancer_bed_prep = Channel.fromPath("$params.store_dir/enhancer.bed").combine(Channel.fromPath("$params.store_dir/enhancer_bg.bed"))
ch_promoter_bed_prep = Channel.fromPath("$params.store_dir/promoter.bed").combine(Channel.fromPath("$params.store_dir/promoter_bg.bed"))

ch_dnaseq = Channel.fromPath("$params.store_dir/enhancer_DNAseq.hg19.npz").combine(Channel.fromPath("$params.store_dir/promoter_DNAseq.hg19.npz"))

Channel.fromPath("$projectDir/co_score_pr.csv")
  .splitCsv(header:true, sep:',')
  .map { row -> [ row.celltype, row.repetition, file("$params.store_dir/$row.npz", checkIfExists: true) ]  }
  .set {ch_pr_co_score}

Channel.fromPath("$projectDir/co_score_enh.csv")
  .splitCsv(header:true, sep:',')
  .map { row -> [ row.celltype, row.repetition, file("$params.store_dir/$row.npz", checkIfExists: true) ]  }
  .set {ch_enh_co_score}

Channel.fromPath("$projectDir/pchic_data.csv")
    .splitCsv(header:true, sep:',')
    .map { row -> [ row.celltype, file("$params.store_dir/$row.hic", checkIfExists: true) ]  }
    .set { ch_hic_aug }

//
workflow combine_data{

    ch_enh = ch_enh_co_score.combine( ch_enhancer_bed_prep )
    ch_pr = ch_pr_co_score.combine( ch_promoter_bed_prep )

    ch_hic_aug_split = SPLIT_HIC_AUG(ch_hic_aug).transpose().take( params.dev ? 2 : -1)

    ch_combined = ch_enh.join(ch_pr, by:[0,1]).combine(ch_hic_aug_split, by:0)
    ch_combined_ = COMBINE_PCHIC_CO_SCORE(ch_combined).groupTuple(by:[0,1])

    // ch_hic_aug_split_.take(1).view()
    // ch_hic_aug_split.view()
    // ch_enh.join(ch_pr, by:[0,1]).view()
    // ch_enh.join(ch_pr, by:[0,1]).combine(ch_hic_aug_split, by:0).view()
    // ch_combined_.view()

    ch_hic_coscore = COMBINE_PCHIC_CO_SCORE_ENH(ch_combined_)
        .join(COMBINE_PCHIC_CO_SCORE_PR(ch_combined_), by:[0,1])
        .groupTuple(by: 0)
    ch_hic_coscore_reps = COMBINE_CO_SCORE_REPS_ENH(ch_hic_coscore).join(COMBINE_CO_SCORE_REPS_PR(ch_hic_coscore), by:0)

    ch__ = ch_hic_aug_split.combine(ch_dnaseq).combine(ch_enhancer_bed_prep).combine(ch_promoter_bed_prep)
    ch_hic_dnaseq = COMBINE_PCHIC_DNA_SEQ(ch__).groupTuple(by:0)

    ch_enh_dnaseq = COMBINE_PCHIC_OUT_ENHANCER(ch_hic_dnaseq) | SAVE_ENH_DNA_SEQ //| UNZIP1
    ch_pr_dnaseq = COMBINE_PCHIC_OUT_PROMOTER(ch_hic_dnaseq) | SAVE_PR_DNA_SEQ //| UNZIP2

    ch_enh_coscore = SAVE_ENH_CO_SCORE(ch_hic_coscore_reps) //| UNZIP3
    ch_pr_coscore = SAVE_PR_CO_SCORE(ch_hic_coscore_reps) //| UNZIP4

    ch_data = ch_enh_coscore
        .join(ch_pr_coscore, by:0)
        .join(ch_enh_dnaseq, by:0)
        .join(ch_pr_dnaseq, by:0)
        .combine(ch_hic_aug_split, by:0)

    // ch_data_gz.view()
    // ch_data = UNZIP(ch_data_gz)
    ch_data.view()
    ch_data_tar = SEPARATE_DATA(ch_data).transpose() | CONVERT_TAR_XZ
    ch_deept_features = COMBINE_DATA_TAR(ch_data_tar.groupTuple(by:0))

    emit: ch_deept_features
}

workflow{
    combine_data()
} 























