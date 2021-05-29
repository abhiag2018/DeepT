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
    storeDir "${params.store_dir}"
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
    storeDir "${params.store_dir}"
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
    path 'hic*.tsv' into ch_hic_parts mode flatten

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
    .into { ch_hic_augment1; ch_hic_augment3 }


// Combining Data : step 0
process SPLIT_HIC_AUG{
    input:
    tuple val(cellType), path(hic_aug) from ch_hic_augment1

    output:
    tuple val(cellType), path("hic.aug.${cellType}.*.csv") into ch_hic_augment1_ mode flatten

    script:
    """
    python -c "import preptools as pt; \
    pt.splitCSV('$hic_aug', \
        [], \
        readArgs = {}, \
        writeArgs = {'header':True}, prefix='hic.aug.${cellType}.', \
        split_num=$params.hic_split_combine, \
        suffix='.csv')"
    """
}

ch_hic_augment1_
    .into { ch_hic_augment2_split; ch_hic_augment1_split }

ch_enhancer_COscore
    .combine( ch_enhancer_bed_prep3 )
    .set{ ch_enhancer_COscore_bedprep }

ch_promoter_COscore
    .combine( ch_promoter_bed_prep3 )
    .set{ ch_promoter_COscore_bedprep }

ch_enhancer_COscore_bedprep
    .join( ch_promoter_COscore_bedprep, by:[0,1] )
    .combine(ch_hic_augment1_split, by:0)
    .set{ ch_hic_COscore }

// Combining Data : step 1.1
// combine PCHi-C interactions to get CO score for each element in the list
process COMBINE_PCHIC_CO_SCORE {
    label 'bigmem'
    afterScript "rm *.npz"

    input:
    tuple val(cellType), val(rep), \
        path(enhancer_COscore), path(enhancer_bed), path(enhancer_bg), \
        path(promoter_COscore), path(promoter_bed), path(promoter_bg), \
        path(hic_aug) from ch_hic_COscore

    output:
    tuple val(cellType), val(rep), path("enhancer_hic.${cellType}.*.rep${rep}.h5.gz"), path("promoter_hic.${cellType}.*.rep${rep}.h5.gz") into ch_hic_features_COscore

    script:  
    """
    python -c "import preptools as pt; \
    import pandas as pd; \
    import numpy as np; \
    import os; \
    num='$hic_aug'.split('.')[-2]; \
    pt.SelectElements_in_TrainingData('$hic_aug', 
        np.load('$promoter_COscore',allow_pickle=True)['expr'],
        np.load('$enhancer_COscore',allow_pickle=True)['expr'],
        pd.read_csv('$promoter_bed', delimiter='\t', names=$params.promoter_headers).name, 
        pd.read_csv('$enhancer_bed', delimiter='\t', names=$params.enhancer_headers).name, 
        f'promoter_hic.${cellType}.{num}.rep${rep}.npz',
        f'enhancer_hic.${cellType}.{num}.rep${rep}.npz', 
        'expr', 
        skip_assert = True)"

    python -c "import numpy as np; \
    import h5py; \
    num='$hic_aug'.split('.')[-2]; \
    A = np.load(f'enhancer_hic.${cellType}.{num}.rep${rep}.npz', allow_pickle=True)['expr']
    with h5py.File(f'enhancer_hic.${cellType}.{num}.rep${rep}.h5', 'w') as h5f:
        h5f.create_dataset('data', data=A)"

    python -c "import numpy as np; \
    import h5py; \
    num='$hic_aug'.split('.')[-2]; \
    A = np.load(f'promoter_hic.${cellType}.{num}.rep${rep}.npz', allow_pickle=True)['expr']
    with h5py.File(f'promoter_hic.${cellType}.{num}.rep${rep}.h5', 'w') as h5f:
        h5f.create_dataset('data', data=A)
    os.system(f'gzip enhancer_hic.${cellType}.{num}.rep${rep}.h5')
    os.system(f'gzip promoter_hic.${cellType}.{num}.rep${rep}.h5')"
    """ 
}

ch_hic_features_COscore
    .groupTuple(by:[0,1])
    .into{ ch_hic_features_COscore_1; ch_hic_features_COscore_2 }

// Combining Data : step 1.2a
// combine PCHi-C interactions to get CO score for each element in the list
process COMBINE_PCHIC_CO_SCORE_ENH {
    label 'bigmem'
    input:
    tuple val(cellType), val(rep), path(enhancer_rep_dat), path(promoter_rep_dat) from ch_hic_features_COscore_1


    output:
    tuple val(cellType), val(rep), path("enhancer_hic.${cellType}.rep${rep}.h5.gz") into ch_hic_features_COscore__1

    script:  
    """
    gzip -d $enhancer_rep_dat
    python -c "import numpy as np; 
    import h5py
    from natsort import natsorted
    L = natsorted('$enhancer_rep_dat'.split())
    h5f_list = list(map(lambda file_name:h5py.File(file_name,'r'), L) )
    h5f = h5py.File('enhancer_hic.${cellType}.rep${rep}.h5','w')
    h5f.create_dataset('data',data=np.zeros((np.sum(list(map(lambda i:h5f_list[i]['data'].shape[0],range(len(L))))),h5f_list[0]['data'].shape[1])))
    row=0;
    for i in range(len(L)):
        s = h5f_list[i]['data'].shape
        h5f['data'][row:row+s[0],:]=h5f_list[i]['data']
        row += s[0]
    list(map(lambda f:f.close(), h5f_list))
    h5f.close()"
    gzip enhancer_hic.${cellType}.*.rep${rep}.h5
    gzip enhancer_hic.${cellType}.rep${rep}.h5
    """ 
}

// Combining Data : step 1.2b
// combine PCHi-C interactions to get CO score for each element in the list
process COMBINE_PCHIC_CO_SCORE_PR {
    label 'bigmem'
    input:
    tuple val(cellType), val(rep), path(enhancer_rep_dat), path(promoter_rep_dat) from ch_hic_features_COscore_2


    output:
    tuple val(cellType), val(rep), path("promoter_hic.${cellType}.rep${rep}.h5.gz") into ch_hic_features_COscore__2

    script:  
    """
    gzip -d $promoter_rep_dat
    python -c "import numpy as np; 
    import h5py
    from natsort import natsorted
    L = natsorted('$promoter_rep_dat'.split())
    h5f_list = list(map(lambda file_name:h5py.File(file_name,'r'), L) )
    h5f = h5py.File('promoter_hic.${cellType}.rep${rep}.h5','w')
    h5f.create_dataset('data',data=np.zeros((np.sum(list(map(lambda i:h5f_list[i]['data'].shape[0],range(len(L))))),h5f_list[0]['data'].shape[1])))
    row=0;
    for i in range(len(L)):
        s = h5f_list[i]['data'].shape
        h5f['data'][row:row+s[0],:]=h5f_list[i]['data']
        row += s[0]
    list(map(lambda f:f.close(), h5f_list))
    h5f.close()"
    gzip promoter_hic.${cellType}.*.rep${rep}.h5
    gzip promoter_hic.${cellType}.rep${rep}.h5
    """ 
}


ch_hic_features_COscore__1
    .join(ch_hic_features_COscore__2, by:[0,1])
    .groupTuple(by: 0)
    .into { ch_hic_features_COscore_reps_1; ch_hic_features_COscore_reps_2 }

// Combining Data : step 1.3
// combine PCHi-C interactions to get CO score for each element in the list
process COMBINE_CO_SCORE_REPS_ENH {
    label 'bigmem'
    publishDir "${params.outdir}/pchic", mode: params.publish_dir_mode

    input:
    tuple val(cellType), val(reps), val(enhancer_COscore_reps), val(promoter_COscore_reps) from ch_hic_features_COscore_reps_1

    output:
    // stdout result
    tuple val(cellType), path("enhancer_hic_COscore.${cellType}.h5") into ch_hic_COscore_features_out_enh

    script:  
    """
    python -c "from natsort import natsorted
    import numpy as np
    import h5py
    import math
    L = natsorted('$enhancer_COscore_reps'[1:-1].split(', '))
    dnase_cell = list(map(lambda file_name:h5py.File(file_name,'r'), L))
    h5f = h5py.File('enhancer_hic_COscore.${cellType}.h5','w')
    h5f.create_dataset('data',data=np.zeros((dnase_cell[0]['data'].shape[0], len(L), dnase_cell[0]['data'].shape[1])))
    load_rows=10000
    for j in range(len(L)):
        for i in range(math.ceil(dnase_cell[0]['data'].shape[0]/load_rows)):
            h5f['data'][load_rows*i:load_rows*(i+1),j,:] = dnase_cell[j]['data'][load_rows*i:load_rows*(i+1),:]
    h5f.close()
    for i in range(len(dnase_cell)):
        dnase_cell[i].close()
    "
    gzip enhancer_hic_COscore.${cellType}.h5
    """ 
}

// Combining Data : step 1.3
// combine PCHi-C interactions to get CO score for each element in the list
process COMBINE_CO_SCORE_REPS_PR {
    label 'bigmem'
    publishDir "${params.outdir}/pchic", mode: params.publish_dir_mode

    input:
    tuple val(cellType), val(reps), val(enhancer_COscore_reps), val(promoter_COscore_reps) from ch_hic_features_COscore_reps_2

    output:
    // stdout result
    tuple val(cellType), path("promoter_hic_COscore.${cellType}.h5") into ch_hic_COscore_features_out_pr

    script:  
    """ 
    python -c "from natsort import natsorted
    import numpy as np
    import h5py
    import math
    L = natsorted('$promoter_COscore_reps'[1:-1].split(', '))
    dnase_cell = list(map(lambda file_name:h5py.File(file_name,'r'), L))
    h5f = h5py.File('promoter_hic_COscore.${cellType}.h5','w')
    h5f.create_dataset('data',data=np.zeros((dnase_cell[0]['data'].shape[0], len(L), dnase_cell[0]['data'].shape[1])))
    load_rows=10000
    for j in range(len(L)):
        for i in range(math.ceil(dnase_cell[0]['data'].shape[0]/load_rows)):
            h5f['data'][load_rows*i:load_rows*(i+1),j,:] = dnase_cell[j]['data'][load_rows*i:load_rows*(i+1),:]
    h5f.close()
    for i in range(len(dnase_cell)):
        dnase_cell[i].close()
    "
    gzip promoter_hic_COscore.${cellType}.h5
    """ 
}
 
ch_hic_COscore_features_out_enh
    .join(ch_hic_COscore_features_out_pr, by:0)
    .set{ ch_hic_COscore_features_out }

ch_hic_augment2_split
    .combine(ch_regElement_DNAseq3)
    .combine(ch_enhancer_bed_prep4)
    .combine(ch_promoter_bed_prep4)
    .set{ ch_hic_DNA_seq }

// Combining Data : step 2.1
// combine PCHi-C interactions to get DNA sequence for each element in the list
process COMBINE_PCHIC_DNA_SEQ {
    label 'bigCpuMem'

    input:
    tuple val(cellType), path(hic_aug), path(enhancer_DNAseq), path(promoter_DNAseq), path(enhancer_bed), path(enhancer_bg), path(promoter_bed), path(promoter_bg)  from ch_hic_DNA_seq

    output:
    tuple val(cellType), path("enhancer_hic.${cellType}.*.DNA_seq.h5.gz"), path("promoter_hic.${cellType}.*.DNA_seq.h5.gz") into ch_hic_DNA_seq_features_out

    script:  
    """
    python -c "import preptools as pt
    import numpy as np
    import pandas as pd
    import h5py
    import os
    reshapeDNA = lambda x:x.reshape(-1,x.shape[-1]//4,4)
    num = '$hic_aug'.split('.')[-2]
    pt.SelectElements_in_TrainingData('$hic_aug', \
            reshapeDNA(np.load('$promoter_DNAseq', allow_pickle=True)['sequence']), \
            reshapeDNA(np.load('$enhancer_DNAseq', allow_pickle=True)['sequence']), \
            pd.read_csv('$promoter_bed', delimiter='\t', names=$params.promoter_headers).name, \
            pd.read_csv('$enhancer_bed', delimiter='\t', names=$params.enhancer_headers).name, \
            f'promoter_hic.${cellType}.{num}.DNA_seq.npz', \
            f'enhancer_hic.${cellType}.{num}.DNA_seq.npz', \
            'sequence', \
            transform = lambda x:x.reshape(-1,x.shape[-1]*x.shape[-2]), \
            skip_assert = True)
    A = np.load(f'enhancer_hic.${cellType}.{num}.DNA_seq.npz', allow_pickle=True)['sequence']
    with h5py.File(f'enhancer_hic.${cellType}.{num}.DNA_seq.h5', 'w') as h5f:
        h5f.create_dataset('data', data=A)
    A = np.load(f'promoter_hic.${cellType}.{num}.DNA_seq.npz', allow_pickle=True)['sequence']
    with h5py.File(f'promoter_hic.${cellType}.{num}.DNA_seq.h5', 'w') as h5f:
        h5f.create_dataset('data', data=A)
    os.system(f'rm enhancer_hic.${cellType}.{num}.DNA_seq.npz')
    os.system(f'rm promoter_hic.${cellType}.{num}.DNA_seq.npz')
    os.system(f'gzip enhancer_hic.${cellType}.{num}.DNA_seq.h5')
    os.system(f'gzip promoter_hic.${cellType}.{num}.DNA_seq.h5')
    "
    """ 
}

ch_hic_DNA_seq_features_out
    .groupTuple(by:0)
    .into { ch_hic_enhDNA_seq_features_out_1; ch_hic_prDNA_seq_features_out_1 }

// Combining Data : step 2.2a
process COMBINE_PCHIC_OUT_ENHANCER {
    label 'bigmem'
    publishDir "${params.outdir}/pchic", mode: params.publish_dir_mode

    input:
    tuple val(cellType), path(enhancer_rep_dat), path(promoter_rep_dat) from ch_hic_enhDNA_seq_features_out_1

    output:
    tuple val(cellType), path("enhancer_hic.${cellType}.DNA_seq.h5.gz") into ch_hic_enhDNA_seq_features_out_1_

    script:
    """
    gzip -d $enhancer_rep_dat
    python -c "import numpy as np; 
    import h5py
    from natsort import natsorted
    L=natsorted('${enhancer_rep_dat}'.split())
    h5f_list = list(map(lambda file_name:h5py.File(file_name,'r'), L) )
    arr_list = list(map(lambda f:f['data'], h5f_list))
    h5f = h5py.File('enhancer_hic.${cellType}.DNA_seq.h5','w')
    h5f.create_dataset('data',data=np.zeros((np.sum(list(map(lambda i:h5f_list[i]['data'].shape[0],range(len(L))))),h5f_list[0]['data'].shape[1])))
    row=0;
    for i in range(len(L)):
        s = h5f_list[i]['data'].shape
        h5f['data'][row:row+s[0],:]=h5f_list[i]['data']
        row += s[0]
    list(map(lambda f:f.close(), h5f_list))
    h5f.close()"
    gzip enhancer_hic.${cellType}.*.DNA_seq.h5
    gzip enhancer_hic.${cellType}.DNA_seq.h5
    """
}

// Combining Data : step 2.2b
process COMBINE_PCHIC_OUT_PROMOTER {
    label 'bigmem'
    publishDir "${params.outdir}/pchic", mode: params.publish_dir_mode

    input:
    tuple val(cellType), path(enhancer_rep_dat), path(promoter_rep_dat) from ch_hic_prDNA_seq_features_out_1

    output:
    tuple val(cellType), path("promoter_hic.${cellType}.DNA_seq.h5.gz") into ch_hic_prDNA_seq_features_out_1_

    script:
    """
    gzip -d $promoter_rep_dat
    python -c "import numpy as np; 
    import h5py
    from natsort import natsorted
    L=natsorted('${promoter_rep_dat}'.split())
    h5f_list = list(map(lambda file_name:h5py.File(file_name,'r'), L) )
    arr_list = list(map(lambda f:f['data'], h5f_list))
    h5f = h5py.File('promoter_hic.${cellType}.DNA_seq.h5','w')
    h5f.create_dataset('data',data=np.zeros((np.sum(list(map(lambda i:h5f_list[i]['data'].shape[0],range(len(L))))),h5f_list[0]['data'].shape[1])))
    row=0;
    for i in range(len(L)):
        s = h5f_list[i]['data'].shape
        h5f['data'][row:row+s[0],:]=h5f_list[i]['data']
        row += s[0]
    list(map(lambda f:f.close(), h5f_list))
    h5f.close()"
    gzip promoter_hic.${cellType}.*.DNA_seq.h5
    gzip promoter_hic.${cellType}.DNA_seq.h5
    """
}

ch_hic_COscore_features_out
    .into{ch_hic_COscore_features_out_1; ch_hic_COscore_features_out_2}

// Combining Data : step 3.1a
// separate the data in ch_hic_DNA_seq_features_out_, ch_hic_COscore_features_out into data points
process SAVE_ENH_DNA_SEQ {
    label 'bigmem'

    input:
    tuple val(cellType), path(enhancer_hic_DNAseq) from ch_hic_enhDNA_seq_features_out_1_

    output:
    tuple val(cellType), path("enhancer_hic.${cellType}.DNA_seq.final.h5.gz") into ch_hic_features_enhDnaSeq

    script:  
    """
    python -c "import numpy as np
    import h5py
    import math
    NUM_SEQ=4
    shape1 = (-1, 1, $params.enhancer_window, NUM_SEQ)
    with h5py.File('$enhancer_hic_DNAseq','r') as enhHicDnaseq_h5:
        s=enhHicDnaseq_h5['data'].shape
        with h5py.File('enhancer_hic.${cellType}.DNA_seq.final.h5', 'w') as h5f:
            h5f.create_dataset('data',(s[0],shape1[1],shape1[3],shape1[2]), dtype=enhHicDnaseq_h5['data'].dtype)
            load_rows = 10000
            for i in range(math.ceil(s[0]/load_rows)):
                h5f['data'][load_rows*i:load_rows*(i+1),:,:,:] = np.array(enhHicDnaseq_h5['data'][load_rows*i:load_rows*(i+1),:]).reshape(shape1).transpose(0, 1, 3, 2)                
    "
    gzip enhancer_hic.${cellType}.DNA_seq.final.h5
    """
}

// Combining Data : step 3.1b
process SAVE_PR_DNA_SEQ {
    label 'bigmem'

    input:
    tuple val(cellType), path(promoter_hic_DNAseq) from ch_hic_prDNA_seq_features_out_1_

    output:
    tuple val(cellType), path("promoter_hic.${cellType}.DNA_seq.final.h5.gz") into ch_hic_features_prDnaSeq

    script:  
    """
    python -c "import numpy as np
    import h5py
    import math
    NUM_SEQ=4
    shape2 = (-1, 1, $params.promoter_window, NUM_SEQ)
    with h5py.File('$promoter_hic_DNAseq','r') as prHicDnaseq_h5:
        s=prHicDnaseq_h5['data'].shape
        with h5py.File(f'promoter_hic.${cellType}.DNA_seq.final.h5', 'w') as h5f:
            h5f.create_dataset('data', (s[0],shape2[1],shape2[3],shape2[2]), dtype=prHicDnaseq_h5['data'].dtype)
            load_rows = 10000
            for i in range(math.ceil(s[0]/load_rows)):
                h5f['data'][load_rows*i:load_rows*(i+1),:,:,:] = np.array(prHicDnaseq_h5['data'][load_rows*i:load_rows*(i+1),:]).reshape(shape2).transpose(0, 1, 3, 2)
    "
    gzip promoter_hic.${cellType}.DNA_seq.final.h5
    """
}

// Combining Data : step 3.1c
process SAVE_ENH_CO_SCORE {
    label 'bigmem'

    input:
    tuple val(cellType), path(enhancer_hic_COscore), path(promoter_hic_COscore) from ch_hic_COscore_features_out_1

    output:
    tuple val(cellType), path("enhancer_hic.${cellType}.CO_score.final.h5.gz") into ch_hic_features_enhCoScore

    script:  
    """
    python -c "import numpy as np
    import h5py
    import math
    with h5py.File('$enhancer_hic_COscore','r') as enhHicCO_h5:
        Tregion1_expr = np.array(enhHicCO_h5['data'])
        NUM_REP = Tregion1_expr.shape[1]
        shape1 = (-1, 1, NUM_REP, $params.enhancer_window)
        Tregion1_expr.reshape(shape1)
        with h5py.File(f'enhancer_hic.${cellType}.CO_score.final.h5', 'w') as h5f:
            h5f.create_dataset('data', data=Tregion1_expr)
    "
    gzip enhancer_hic.${cellType}.CO_score.final.h5
    """
}

// Combining Data : step 3.1d
process SAVE_PR_CO_SCORE {
    label 'bigmem'

    input:
    tuple val(cellType), path(enhancer_hic_COscore), path(promoter_hic_COscore) from ch_hic_COscore_features_out_2

    output:
    tuple val(cellType), path("promoter_hic.${cellType}.CO_score.final.h5.gz") into ch_hic_features_prCoScore

    script:  
    """
    python -c "import numpy as np
    import h5py
    import math
    with h5py.File('$promoter_hic_COscore','r') as prHicCO_h5:
        Tregion2_expr = np.array(prHicCO_h5['data'])
        NUM_REP = Tregion2_expr.shape[1]
        shape2 = (-1, 1, NUM_REP, $params.promoter_window)
        Tregion2_expr.reshape(shape2)
        with h5py.File(f'promoter_hic.${cellType}.CO_score.final.h5', 'w') as h5f:
            h5f.create_dataset('data', data=Tregion2_expr)
    "
    gzip promoter_hic.${cellType}.CO_score.final.h5
    """
}


ch_hic_features_enhCoScore
    .join(ch_hic_features_prCoScore, by:0)
    .join(ch_hic_features_enhDnaSeq, by:0)
    .join(ch_hic_features_prDnaSeq, by:0)
    .join( ch_hic_augment3, by:0 )
    .set{ ch_hic_features_out }
    // .view()


// Combining Data : step 3.2
// separate the data in ch_hic_DNA_seq_features_out_, ch_hic_COscore_features_out into data points
process SEPARATE_DATA {
    label 'bigmem'

    input:
    tuple val(cellType), path(enhancer_hic_CO), path(promoter_hic_CO), path(enhancer_hic_DNAseq), path(promoter_hic_DNAseq), path(hic_aug) from ch_hic_features_out

    output:
    tuple val(cellType), path("input_feature_part_*.h5.gz") into ch_hic_features_split_out mode flatten

    script:  
    """
    gzip -d $enhancer_hic_CO $promoter_hic_CO $enhancer_hic_DNAseq $promoter_hic_DNAseq
    python -c "import pandas as pd
    import h5py
    import math

    Tlabel = pd.read_csv('$hic_aug')['label']

    h5f_enhDNAseq = h5py.File('$enhancer_hic_DNAseq','r')
    Tregion1_seq = h5f_enhDNAseq['data']

    h5f_prDNAseq = h5py.File('$promoter_hic_DNAseq','r')
    Tregion2_seq = h5f_prDNAseq['data']


    ## load data: DNase
    h5f_enhDNase = h5py.File('$enhancer_hic_CO','r')
    Tregion1_expr = h5f_enhDNase['data']

    h5f_prDNase = h5py.File('$promoter_hic_CO','r')
    Tregion2_expr = h5f_prDNase['data']



    NUM = Tlabel.shape[0]

    step=math.ceil(NUM/$params.sepdata_split)
    for i in range(0,NUM,step):
        idx=i//step
        print(idx, flush=True)
        with h5py.File(f'input_feature_part_{idx}.h5','w') as h5f:
            h5f.create_dataset('enhseq',data=Tregion1_seq[i:min(NUM,i+step),:,:,:])
            h5f.create_dataset('prseq',data=Tregion2_seq[i:min(NUM,i+step),:,:,:])
            h5f.create_dataset('enhdnase',data=Tregion1_expr[i:min(NUM,i+step),:,:])
            h5f.create_dataset('prdnase',data=Tregion2_expr[i:min(NUM,i+step),:,:])
            h5f.create_dataset('label',data=Tlabel[i:min(NUM,i+step)])
            h5f.create_dataset('index',data=i)
    h5f_enhDNAseq.close()
    h5f_prDNAseq.close()
    h5f_enhDNase.close()
    h5f_prDNase.close()"
    gzip input_feature_part_*.h5
    """ 
}

// ch_hic_features_split_out_ = ch_hic_features_split_out.take( params.dev ? params.dev_lim_tar : -1 )


// Combining Data : step 3.2
// separate the data in ch_hic_DNA_seq_features_out_, ch_hic_COscore_features_out into data points
process CONVERT_TAR_XZ {
    label 'bigCpuMem'
    errorStrategy 'finish'

    input:
    tuple val(cellType), path(input_feature_part) from ch_hic_features_split_out

    output:
    tuple val(cellType), path("features.${cellType}.*.tar.xz") into ch_features_xz

    script:  
    """
    gzip -d input_feature_part
    python -c "import os
    import multiprocessing
    import h5py
    import numpy as np
    h5f = h5py.File('$input_feature_part','r')
    def npz_save(out_npz, enh_seq_data, pr_seq_data, enh_dnase_data, pr_dnase_data, label_data):
        np.savez(out_npz, enh_seq = enh_seq_data , pr_seq = pr_seq_data , enh_dnase = enh_dnase_data , pr_dnase = pr_dnase_data, label = label_data)
        return 0
    with multiprocessing.Pool() as pool:
        i=h5f['index'][()]
        outdir=f'data_{i}'
        os.makedirs(outdir,exist_ok=True)
        pool.starmap(npz_save, zip(map(lambda index:f'./{outdir}/{index}.npz', range(i,i+h5f['label'].shape[0])), h5f['enhseq'], h5f['prseq'], h5f['enhdnase'], h5f['prdnase'], h5f['label']))
    h5f.close()
    os.system(f'tar cJvfh features.${cellType}.{i}.tar.xz data_{i}')"
    """ 
}

ch_hic_features_xz_out=ch_features_xz.groupTuple(by:0)

// Combining Data : step 3.4
// separate the data in ch_hic_DNA_seq_features_out_, ch_hic_COscore_features_out into data points
process COMBINE_DATA_TAR {
    label 'bigmem'
    stageInMode 'copy'
    publishDir "${params.outdir}/pchic", mode: params.publish_dir_mode

    input:
    tuple val(cellType), path(input_xz) from ch_hic_features_xz_out

    output:
    tuple val(cellType), path("features.${cellType}.tar.xz") into ch_deept_features

    script:  
    """
    xzcat $input_xz | xz -c > features.${cellType}.tar.xz 
    """ 
    // uncompress command
    // tar --ignore-zeros -xJvf features.${cellType}.tar.xz
}


ch_deept_features.view()
















