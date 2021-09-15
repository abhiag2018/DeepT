nextflow.preview.dsl=2
// include { promoter_bed; enhancer_bed } from "$projectDir/modules/pr-enh-prep" // step 1
// include { prep_co_score_enh; prep_co_score_pr } from "$projectDir/modules/co-score-prep" // step 2
// include { prep_gen_seq } from "$projectDir/modules/genome-seq-prep" // step 3
// include { prep_pchic } from "$projectDir/modules/pchic-prep" // step 4

include {SPLIT_HIC_AUG; COMBINE_PCHIC_CO_SCORE; COMBINE_PCHIC_CO_SCORE_ENH; COMBINE_PCHIC_CO_SCORE_PR; \
    COMBINE_CO_SCORE_REPS_ENH; COMBINE_CO_SCORE_REPS_PR; \
    COMBINE_PCHIC_DNA_SEQ; COMBINE_PCHIC_OUT_ENHANCER; COMBINE_PCHIC_OUT_PROMOTER; \
    SEPARATE_DATA} from "$projectDir/modules/combine.v1" 
include {splitByChromosomeBed; splitByChromosomeHiC; splitByChromosomeHiCcross; splitByChromosomeCOscore; splitByChromosomeDNAseq} from "$projectDir/modules/splitByChrom"

if (params.species == "hg" ){
    ch_enhancer_bed_prep = Channel.fromPath("$params.store_dir/hg19-enh-pr/enhancer.bed").combine(Channel.fromPath("$params.store_dir/hg19-enh-pr/enhancer_bg.bed"))
    ch_promoter_bed_prep = Channel.fromPath("$params.store_dir/hg19-enh-pr/promoter.bed").combine(Channel.fromPath("$params.store_dir/hg19-enh-pr/promoter_bg.bed"))
    // feature input files
    // features for promoter and enhancer list (in the same order as the list)
    ch_dnaseq = Channel.fromPath("$params.store_dir/enhancer_DNAseq.hg19.npz").combine(Channel.fromPath("$params.store_dir/promoter_DNAseq.hg19.npz"))
}
else if (params.species == "mm" ){
    ch_enhancer_bed_prep = Channel.fromPath("$params.store_dir/mm9-enh-pr/enhancer.bed").combine(Channel.fromPath("$params.store_dir/mm9-enh-pr/enhancer_bg.bed"))
    ch_promoter_bed_prep = Channel.fromPath("$params.store_dir/mm9-enh-pr/promoter.bed").combine(Channel.fromPath("$params.store_dir/mm9-enh-pr/promoter_bg.bed"))
    // feature input files
    // features for promoter and enhancer list (in the same order as the list)
    ch_dnaseq = Channel.fromPath("$params.store_dir/enhancer_DNAseq.mm9.npz").combine(Channel.fromPath("$params.store_dir/promoter_DNAseq.mm9.npz"))
}
else{
    println "species (: $params.species) should be hg or mm; "
    exit 1
}


Channel.fromPath("$projectDir/co_score_pr.csv")
  .splitCsv(header:true, sep:',')
  .map { row -> [ row.celltype, row.repetition, file("$params.store_dir/$row.npz", checkIfExists: true) ]  }
  .set {ch_pr_co_score}

Channel.fromPath("$projectDir/co_score_enh.csv")
  .splitCsv(header:true, sep:',')
  .map { row -> [ row.celltype, row.repetition, file("$params.store_dir/$row.npz", checkIfExists: true) ]  }
  .set {ch_enh_co_score}

if( params.dtype=="val" ){
    dtype = "data_val"
    pchic_data = "pchic_data-val.csv"
}
else if( params.dtype=="test" ){
    dtype = "data_test"
    pchic_data = "pchic_data-test.csv"
}
else if( params.dtype=="train" ){
    dtype = "data_train"
    pchic_data = "pchic_data-train.csv"
}
else if( params.dtype=="all" ){
    dtype = "data_all"
    pchic_data = "pchic_data-all.csv"
}

Channel.fromPath("$projectDir/$pchic_data")
    .splitCsv(header:true, sep:',')
    .map { row -> [ row.celltype, file("$params.store_dir/$row.hic", checkIfExists: true) ]  }
    .set { ch_hic_aug }

workflow combine_data{
    take:
    ch_enh_co_score
    ch_enhancer_bed_prep
    ch_pr_co_score
    ch_promoter_bed_prep
    ch_dnaseq
    ch_hic_aug

    main:
    ch_enh = ch_enh_co_score.combine( ch_enhancer_bed_prep )
    ch_pr = ch_pr_co_score.combine( ch_promoter_bed_prep )

    ch_hic_aug_split = SPLIT_HIC_AUG(ch_hic_aug).transpose()//.take( params.dev ? 2 : -1)

    ch_combined = ch_enh.join(ch_pr, by:[0,1]).combine(ch_hic_aug_split, by:0)
    ch_combined_ = COMBINE_PCHIC_CO_SCORE(ch_combined).groupTuple(by:[0,1])

    ch_hic_coscore = COMBINE_PCHIC_CO_SCORE_ENH(ch_combined_)
        .join(COMBINE_PCHIC_CO_SCORE_PR(ch_combined_), by:[0,1])
        .groupTuple(by: 0)
    ch_enh_coscore = COMBINE_CO_SCORE_REPS_ENH(ch_hic_coscore)
    ch_pr_coscore = COMBINE_CO_SCORE_REPS_PR(ch_hic_coscore)

    ch__ = ch_hic_aug_split.combine(ch_dnaseq).combine(ch_enhancer_bed_prep).combine(ch_promoter_bed_prep)
    ch_hic_dnaseq = COMBINE_PCHIC_DNA_SEQ(ch__).groupTuple(by:0)

    ch_enh_dnaseq = COMBINE_PCHIC_OUT_ENHANCER(ch_hic_dnaseq) 
    ch_pr_dnaseq = COMBINE_PCHIC_OUT_PROMOTER(ch_hic_dnaseq) 

    ch_data = ch_enh_coscore
        .join(ch_pr_coscore, by:0)
        .join(ch_enh_dnaseq, by:0)
        .join(ch_pr_dnaseq, by:0)
        .combine(ch_hic_aug, by:0)

    ch_deept_features = SEPARATE_DATA(ch_data)

    emit: ch_deept_features
}

workflow splitByChrom{
    take:
        ch_chr
        ch_enh_co_score
        ch_enhancer_bed_prep
        ch_pr_co_score
        ch_promoter_bed_prep
        ch_dnaseq
        ch_hic_aug
        data_type

    main:
    ch_hic_aug_ChromSplit = ch_chr.combine(ch_hic_aug) | splitByChromosomeHiC

    ch_bed = Channel.value('enh').combine(ch_enhancer_bed_prep.flatten()).concat(Channel.value('pr').combine(ch_promoter_bed_prep.flatten()))

    ch_bedsplit = splitByChromosomeBed(ch_chr.combine(ch_bed))

    ch_bedsplit_filtered = ch_bedsplit
        .filter{ it[4] == 'enhancer_bg' ||  it[4] == 'promoter_bg' }
        .map( it-> [it[0], it[1], it[3]])


    ch_dnaseq_ = ch_chr.combine(Channel.value(['enh','pr']).flatten().merge(ch_dnaseq.flatten()))
    ch_inp_dnaseq_split = ch_bedsplit_filtered.concat(ch_dnaseq_).groupTuple(by:[0,1])
        .map(it -> [it[0], it[1], it[2][0], it[2][1]])


    ch_dnaseq_split = splitByChromosomeDNAseq(ch_inp_dnaseq_split)

    ch_dnaseq_split_enh = ch_dnaseq_split.filter{ it[1]=="enh" }.map(it -> [it[0], it[2]])
    ch_dnaseq_split_pr = ch_dnaseq_split.filter{ it[1]=="pr" }.map(it -> [it[0], it[2]])
    ch_dnaseq_split_ = ch_dnaseq_split_enh.join(ch_dnaseq_split_pr, by:0)

    ch_co_score = Channel.value('enh').combine(ch_enh_co_score)
        .concat(Channel.value('pr').combine(ch_pr_co_score))

    ch_coscore_split = ch_bedsplit_filtered.combine(ch_chr.combine(ch_co_score), by:[0,1]) | splitByChromosomeCOscore
    ch_coscore_split_enh = ch_coscore_split.filter{ it[1]=="enh" }.map(it -> [it[0], it[2], it[3], it[4]])
    ch_coscore_split_pr = ch_coscore_split.filter{ it[1]=="pr" }.map(it -> [it[0], it[2], it[3], it[4]])

    ch_bedsplit_bg = ch_bedsplit.filter{ it[4] == 'enhancer_bg' ||  it[4] == 'promoter_bg' }
    ch_bedsplit_main = ch_bedsplit.filter{ it[4] == 'enhancer' ||  it[4] == 'promoter' }

    ch_enh_bg_enh = ch_bedsplit_bg.filter{ it[1] == "enh" }
        .map( it -> [it[0], it[2]] )
        
    ch_enh_bg_pr = ch_bedsplit_bg.filter{ it[1] == "pr" }
        .map( it -> [it[0], it[2]] )
        
    ch_bedsplit_enh = ch_bedsplit_main.filter{ it[1] == "enh" }
        .map( it -> [it[0], it[2]] )
        .join(ch_enh_bg_enh, by:0)


    ch_bedsplit_pr = ch_bedsplit_main.filter{ it[1] == "pr" }
        .map( it -> [it[0], it[2]] )
        .join(ch_enh_bg_pr, by:0)

    deeptact_features = combine_data(ch_coscore_split_enh.map( it -> [it[1], it[2], it[3]] ),
        ch_bedsplit_enh.map( it -> [it[1], it[2]] ),
        ch_coscore_split_pr.map( it -> [it[1], it[2], it[3]] ),
        ch_bedsplit_pr.map( it -> [it[1], it[2]] ),
        ch_dnaseq_split_.map( it -> [it[1], it[2]] ),
        ch_hic_aug_ChromSplit.map( it -> [it[1], it[2]] ) )

    features_out = RENAME_OUT(ch_chr.combine(data_type).combine(deeptact_features))

    emit: features_out
}

process RENAME_OUT {
    storeDir "${params.store_dir}"
    errorStrategy 'finish'

    input:
    tuple val(chrom), val(dataType), val(cellType), path(data_gz) 

    output:
    path("${dataType}.${chrom}.${cellType}.tar.gz")

    // myFile = file("$data_gz")
    // myFile.renameTo("${dataType}.${chrom}.${cellType}.tar.gz")
    script: 
    """
    ln -s `readlink -f $data_gz` ${dataType}.${chrom}.${cellType}.tar.gz
    """
}

workflow chr1{ splitByChrom(Channel.value("chr1"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chr2{ splitByChrom(Channel.value("chr2"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chr3{ splitByChrom(Channel.value("chr3"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chr4{ splitByChrom(Channel.value("chr4"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chr5{ splitByChrom(Channel.value("chr5"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chr6{ splitByChrom(Channel.value("chr6"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chr7{ splitByChrom(Channel.value("chr7"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chr8{ splitByChrom(Channel.value("chr8"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chr9{ splitByChrom(Channel.value("chr9"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chr10{ splitByChrom(Channel.value("chr10"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chr11{ splitByChrom(Channel.value("chr11"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chr12{ splitByChrom(Channel.value("chr12"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chr13{ splitByChrom(Channel.value("chr13"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chr14{ splitByChrom(Channel.value("chr14"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chr15{ splitByChrom(Channel.value("chr15"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chr16{ splitByChrom(Channel.value("chr16"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chr17{ splitByChrom(Channel.value("chr17"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chr18{ splitByChrom(Channel.value("chr18"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chr19{ splitByChrom(Channel.value("chr19"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chr20{ splitByChrom(Channel.value("chr20"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chr21{ splitByChrom(Channel.value("chr21"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chr22{ splitByChrom(Channel.value("chr22"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chrX{ splitByChrom(Channel.value("chrX"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }
workflow chrY{ splitByChrom(Channel.value("chrY"), ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug, Channel.value(dtype)).view() }

workflow {
    // chr_list = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
    // ch_chr = params.dev ? Channel.value("chr22") : Channel.value(chr_list).flatten()
    if (params.dev && params.splitByChr){
        chr22()
    }else if (params.splitByChr){            
        ch_chr1 = chr1()
        ch_chr2 = chr2()
        ch_chr3 = chr3()
        ch_chr4 = chr4()
        ch_chr5 = chr5()
        ch_chr6 = chr6()
        ch_chr7 = chr7()
        ch_chr8 = chr8()
        ch_chr9 = chr9()
        ch_chr10 = chr10()
        ch_chr11 = chr11()
        ch_chr12 = chr12()
        ch_chr13 = chr13()
        ch_chr14 = chr14()
        ch_chr15 = chr15()
        ch_chr16 = chr16()
        ch_chr17 = chr17()
        ch_chr18 = chr18()
        ch_chr19 = chr19()
        ch_chr20 = chr20()
        ch_chr21 = chr21()
        ch_chr22 = chr22()
        ch_chrX = chrX()
        ch_chrY = chrY()

        data_cross=combine_data(ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, splitByChromosomeHiCcross(ch_hic_aug))
        chr_cross = RENAME_OUT(Channel.value("chrCROSS").combine(Channel.value(dtype)).combine(data_cross))

    }else{
        combine_data(ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug)
    }
}














