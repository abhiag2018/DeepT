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

// process to filter bed file by chromosome 
process splitByChromosomeHiC {
    memory '5 GB'

    input:
    tuple val(chr), val(cellType), path(hic_aug) 

    output:
    tuple val(chr), val(cellType), path("${hic_aug}.$chr")


    script:
    """
    head -n 1 $hic_aug > ${hic_aug}.$chr
    tail -n +2 $hic_aug | awk -F',' '{ if(\$1 == "$chr" && \$5 == "$chr") { print }}' >> ${hic_aug}.$chr
    """
}

// process to filter interactions with diff chromosome on diff anchors
process splitByChromosomeHiCcross {
    memory '5 GB'

    input:
    tuple val(cellType), path(hic_aug) 

    output:
    tuple val(cellType), path("${hic_aug}.X")


    script:
    """
    head -n 1 $hic_aug > ${hic_aug}.X
    tail -n +2 $hic_aug | awk -F',' '{ if(\$1 != \$5) { print }}' >> ${hic_aug}.X
    """
}


process splitByChromosomeCOscore {
    memory '5 GB'

    input:
    tuple val(chr), val(enh_pr), path(index_file), val(cellType), val(rep), path(numpy_coscore)

    output:
    tuple val(chr), val(enh_pr), val(cellType), val(rep), path("${numpy_coscore.baseName}.${chr}.npz")


    script:
    """
    #!/usr/bin/env python
    import numpy as np
    import pandas as pd
    index=list(pd.read_csv("$index_file",header=None)[0])
    coscore_data=np.load("$numpy_coscore")
    np.savez("${numpy_coscore.baseName}.${chr}.npz",expr=coscore_data['expr'][index,:])
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

// Combining Data : step 5.3.2
// separate the data into NPZ files.. into $sepdata_split chunks for easy processing in next step
// output : input_feature_part_{0-511}.h5.gz
// time: 1h40m
process SEPARATE_DATA_NPZ {
    clusterOptions = '--qos batch --cpus-per-task 8'
    errorStrategy 'finish'

    input:
    tuple val(cellType), path(enhancer_hic_CO_gz), path(promoter_hic_CO_gz), path(enhancer_hic_DNAseq_gz), path(promoter_hic_DNAseq_gz), path(hic_aug) 

    output:
    tuple val(cellType), path("data_chr.${cellType}.tar.gz")

    script:  
    enhancer_hic_CO = enhancer_hic_CO_gz.baseName
    promoter_hic_CO = promoter_hic_CO_gz.baseName
    enhancer_hic_DNAseq = enhancer_hic_DNAseq_gz.baseName
    promoter_hic_DNAseq = promoter_hic_DNAseq_gz.baseName
    """
    gzip -df $enhancer_hic_CO_gz
    gzip -df $promoter_hic_CO_gz
    gzip -df $enhancer_hic_DNAseq_gz
    gzip -df $promoter_hic_DNAseq_gz
    python -c "import h5py
    import os
    import pandas as pd
    from multiprocessing import Pool
    import numpy as np

    Tlabel = pd.read_csv('$hic_aug')['label']
    Tenh_coscore = h5py.File('$enhancer_hic_CO_gz.baseName','r')['data']
    Tenh_dnaseq = h5py.File('$enhancer_hic_DNAseq_gz.baseName','r')['data']
    Tpr_coscore = h5py.File('$promoter_hic_CO_gz.baseName','r')['data']
    Tpr_dnaseq = h5py.File('$promoter_hic_DNAseq_gz.baseName','r')['data']

    def npz_save(out, enh_coscore_, enh_dnaseq_, pr_coscore_, pr_dnaseq_, label_): 
        np.savez(out, enh_seq = enh_dnaseq_ , 
            pr_seq = pr_dnaseq_ , 
            enh_dnase = enh_coscore_[np.newaxis,:,:] , 
            pr_dnase = pr_coscore_[np.newaxis,:,:] , 
            label = label_)

    os.makedirs('data')
    with Pool() as pool:
        pool.starmap(npz_save, zip([f'data/{i}.npz' for i in range(len(Tlabel))], Tenh_coscore, Tenh_dnaseq, Tpr_coscore, Tpr_dnaseq, Tlabel ))
    "
    # rm $enhancer_hic_CO_gz.baseName $promoter_hic_CO_gz.baseName $enhancer_hic_DNAseq_gz.baseName $promoter_hic_DNAseq_gz.baseName
    tar -czvf data_chr.${cellType}.tar.gz data
    rm -rf data
    """
}

process RENAME_OUT {
    storeDir "${params.store_dir}"
    errorStrategy 'finish'

    input:
    tuple val(chrom), val(cellType), path(data_gz) 

    output:
    path("data.${chrom}.${cellType}.tar.gz")

    script: 
    """
    mv $data_gz data.${chrom}.${cellType}.tar.gz
    """

}

ch_nCD4 = Channel.value('nCD4')
    .combine(Channel.fromPath( '/projects/li-lab/agarwa/CUBE/DeepTact/code/storeDir/nCD4.hic.*/input_feature_part_*.h5.gz' ))
ch_tB = Channel.value('tB')
    .combine(Channel.fromPath( '/projects/li-lab/agarwa/CUBE/DeepTact/code/storeDir/tB.hic.*/input_feature_part_*.h5.gz' ))

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

    ch_deept_features = SEPARATE_DATA_NPZ(ch_data)

    // // ch_data = UNZIP(ch_data_gz)
    // ch_data_tar = SEPARATE_DATA(ch_data).transpose() | CONVERT_TAR_XZ
    // // ch_data_tar.take(3).view()
    // ch_deept_features = COMBINE_DATA_TAR(ch_data_tar.groupTuple(by:0))

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

    emit: deeptact_features
}

workflow chrom_features{   
    take: ch_chr

    main:
    // NOTE
    // assert ch_chr.count() == 1 // since splitByChrom can only handle one input; otherwise the channels get mixed up
    features_chr = splitByChrom(ch_chr, ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug)
    features_out = RENAME_OUT(ch_chr.combine(features_chr))

    emit: features_out
}


workflow chr1{   chrom_features(Channel.value("chr1")).view() }
workflow chr2{   chrom_features(Channel.value("chr2")).view() }
workflow chr3{   chrom_features(Channel.value("chr3")).view() }
workflow chr4{   chrom_features(Channel.value("chr4")).view() }
workflow chr5{   chrom_features(Channel.value("chr5")).view() }
workflow chr6{   chrom_features(Channel.value("chr6")).view() }
workflow chr7{   chrom_features(Channel.value("chr7")).view() }
workflow chr8{   chrom_features(Channel.value("chr8")).view() }
workflow chr9{   chrom_features(Channel.value("chr9")).view() }
workflow chr10{   chrom_features(Channel.value("chr10")).view() }
workflow chr11{   chrom_features(Channel.value("chr11")).view() }
workflow chr12{   chrom_features(Channel.value("chr12")).view() }
workflow chr13{   chrom_features(Channel.value("chr13")).view() }
workflow chr14{   chrom_features(Channel.value("chr14")).view() }
workflow chr15{   chrom_features(Channel.value("chr15")).view() }
workflow chr16{   chrom_features(Channel.value("chr16")).view() }
workflow chr17{   chrom_features(Channel.value("chr17")).view() }
workflow chr18{   chrom_features(Channel.value("chr18")).view() }
workflow chr19{   chrom_features(Channel.value("chr19")).view() }
workflow chr20{   chrom_features(Channel.value("chr20")).view() }
workflow chr21{   chrom_features(Channel.value("chr21")).view() }
workflow chr22{   chrom_features(Channel.value("chr22")).view() }
workflow chrX{   chrom_features(Channel.value("chrX")).view() }
workflow chrY{   chrom_features(Channel.value("chrY")).view() }


workflow {
    // chr_list = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
    // ch_chr = params.dev ? Channel.value("chr22") : Channel.value(chr_list).flatten()
    if (params.dev){
        chr22()
    }else{            
        chr1()
        chr2()
        chr3()
        chr4()
        chr5()
        chr6()
        chr7()
        chr8()
        chr9()
        chr10()
        chr11()
        chr12()
        chr13()
        chr14()
        chr15()
        chr16()
        chr17()
        chr18()
        chr19()
        chr20()
        chr21()
        chr22()
        chrX()
        chrY()
        combine_data(ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, splitByChromosomeHiCcross(ch_hic_aug)).view()
    }
}



















