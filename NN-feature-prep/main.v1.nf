nextflow.preview.dsl=2
// include { promoter_bed; enhancer_bed } from "$projectDir/modules/pr-enh-prep" // step 1
// include { prep_co_score_enh; prep_co_score_pr } from "$projectDir/modules/co-score-prep" // step 2
// include { prep_gen_seq } from "$projectDir/modules/genome-seq-prep" // step 3
// include { prep_pchic } from "$projectDir/modules/pchic-prep" // step 4

include {SPLIT_HIC_AUG; COMBINE_PCHIC_CO_SCORE; COMBINE_PCHIC_CO_SCORE_ENH; COMBINE_PCHIC_CO_SCORE_PR; \
    COMBINE_CO_SCORE_REPS_ENH; COMBINE_CO_SCORE_REPS_PR; \
    COMBINE_PCHIC_DNA_SEQ; COMBINE_PCHIC_OUT_ENHANCER; COMBINE_PCHIC_OUT_PROMOTER; \
    SAVE_ENH_DNA_SEQ; SAVE_PR_DNA_SEQ; SAVE_ENH_CO_SCORE; SAVE_PR_CO_SCORE; SEPARATE_DATA; CONVERT_TAR_XZ; COMBINE_DATA_TAR} from "$projectDir/modules/combine.v1" 


ch_enhancer_bed_prep = Channel.fromPath("$params.store_dir/enhancer.bed").combine(Channel.fromPath("$params.store_dir/enhancer_bg.bed"))
ch_promoter_bed_prep = Channel.fromPath("$params.store_dir/promoter.bed").combine(Channel.fromPath("$params.store_dir/promoter_bg.bed"))

// feature input files
// features for promoter and enahncer list (in the same order as the list)
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


// process to filter bed file by chromosome 
// outputs the filtered file and the index file for the filtered rows
process splitByChromosomeBed {
    memory '5 GB'

    input:
    tuple val(chr), val(enh_pr), path(bed_file)

    output:
    tuple val(chr), val(enh_pr), path("${bed_file}.$chr"), path("${bed_file}.index"), val(bed_file.baseName)


    script:
    """
    num=`wc -l $bed_file | awk '{print \$1}'`
    seq 0 \$((\$num-1)) > tmp
    paste $bed_file tmp > ${bed_file}.tmp
    awk '{ if(\$1 == "$chr") { print } }' ${bed_file}.tmp > ${bed_file}.tmp.chr
    awk 'NF-=1' OFS='\t' ${bed_file}.tmp.chr > ${bed_file}.$chr
    awk '{ print \$NF }' ${bed_file}.tmp.chr > ${bed_file}.index
    """
}


process splitByChromosomeCOscore {
    memory '5 GB'

    input:
    tuple val(chr), val(enh_pr), path(index_file), val(cellType), val(rep), path(numpy_coscore)

    output:
    tuple val(chr), val(enh_pr), path("${numpy_coscore.baseName}.${chr}.npz")


    script:
    """
    #!/usr/bin/env python
    import numpy as np
    import pandas as pd
    index=list(pd.read_csv("$index_file",header=None)[0])
    coscore_data=np.load("$numpy_coscore")
    np.savez("${numpy_coscore.baseName}.${chr}.npz",expr=coscore_data['sequence'][index,:])
    """
}

process  splitByChromosomeDNAseq{
    memory '5 GB'

    input:
    tuple val(chr), val(enh_pr), path(index_file), path(numpy_dna)

    output:
    tuple val(chr), val(enh_pr), path("${numpy_dna.baseName}.${chr}.npz")


    script:
    """
    #!/usr/bin/env python
    import numpy as np
    import pandas as pd
    index=list(pd.read_csv("$index_file",header=None)[0])
    dna_data=np.load("$numpy_dna")
    np.savez("${numpy_dna.baseName}.${chr}.npz",sequence=dna_data['sequence'][index,:],name=dna_data['name'][index],loc=dna_data['loc'][index])
    """
}

ch_nCD4 = Channel.value('nCD4')
    .combine(Channel.fromPath( '/projects/li-lab/agarwa/CUBE/DeepTact/code/storeDir/nCD4.hic.*/input_feature_part_*.h5.gz' ))
ch_tB = Channel.value('tB')
    .combine(Channel.fromPath( '/projects/li-lab/agarwa/CUBE/DeepTact/code/storeDir/tB.hic.*/input_feature_part_*.h5.gz' ))

workflow combine_data{

    ch_enh = ch_enh_co_score.combine( ch_enhancer_bed_prep )
    ch_pr = ch_pr_co_score.combine( ch_promoter_bed_prep )

    ch_hic_aug_split = SPLIT_HIC_AUG(ch_hic_aug).transpose().take( params.dev ? 2 : -1)

    ch_combined = ch_enh.join(ch_pr, by:[0,1]).combine(ch_hic_aug_split, by:0)
    ch_combined_ = COMBINE_PCHIC_CO_SCORE(ch_combined).groupTuple(by:[0,1])

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
        .combine(ch_hic_aug, by:0)

    // ch_data = UNZIP(ch_data_gz)
    ch_data_tar = SEPARATE_DATA(ch_data).transpose() | CONVERT_TAR_XZ
    // ch_data_tar.take(3).view()
    ch_deept_features = COMBINE_DATA_TAR(ch_data_tar.groupTuple(by:0))

    emit: ch_deept_features
}


workflow{
    // ch_data_tar = CONVERT_TAR_XZ(ch_nCD4.concat(ch_tB))
    // ch_deept_features = COMBINE_DATA_TAR(ch_data_tar.groupTuple(by:0))
    // ch_deept_features.view()

    chr_list = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
    ch_chr = Channel.value(chr_list).flatten().take( params.dev ? 2 : -1)
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

    ch_bedsplit_filtered.join(ch_chr.combine(Channel.value('enh').combine(ch_enh_co_score).concat(Channel.value('pr').combine(ch_pr_co_score))), by:[0,1])
    left = Channel.from(['X', 1], ['X', 1.1], ['Y', 2], ['Z', 3], ['P', 7])
    right= Channel.from(['Z', 6], ['Y', 5], ['X', 4])
    left.cross(right).count().view()

    source = Channel.from( [1, 'alpha'], [2, 'beta'] )
    target = Channel.from( [1, 'x'], [1, 'y'], [1, 'z'], [2,'p'], [2,'q'], [2,'t'] )

    source.cross(target).view()

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

} 























