nextflow.preview.dsl=2
include { PREPROCESS_PROMOTER; PREPROCESS_ENHANCER } from "${params.baseDir}/code/NN-feature-prep/modules/pr-enh-prep"

ch_promoter_bed = Channel.fromPath( params.promoter_bedfile )
ch_enhancer_bed = Channel.fromPath( params.enhancer_bedfile )

ch_input_fasta = Channel.fromPath(params.species_genome_fasta)

// promoter genomic sequence computation 
process GENOMIC_SEQUENCE_PROMOTER {
    storeDir "${params.store_dir}"

    input:
    tuple path(promoter), path(promoter_bg)
    path(fasta)

    output:
    file("promoter_DNAseq.npz.gz")

    script:  
    """
    python -c "import preptools as pt; \
    pt.generate_elem_fa('$fasta', \
        '$promoter', \
        out_fa = 'promoter.fa' )"

    python -c "import preptools as pt; \
    pt.fasta_to_onehot('promoter.fa', \
        outp='promoter_DNAseq.npz')"
    gzip promoter_DNAseq.npz
    """
}



// enhancer genomic sequence computation 
process GENOMIC_SEQUENCE_ENHANCER {
    storeDir "${params.store_dir}"

    input:
    tuple path(enhancer), path(enhancer_bg)
    path(fasta)

    output:
    file("enhancer_DNAseq.npz.gz")

    script:  
    """
    python -c "import preptools as pt; \
    pt.generate_elem_fa('$fasta', \
        '$enhancer', \
        out_fa = 'enhancer.fa' )"

    python -c "import preptools as pt; \
    pt.fasta_to_onehot('enhancer.fa', \
        outp='enhancer_DNAseq.npz')"
    gzip enhancer_DNAseq.npz
    """
}



workflow prep_gen_seq{
    ch_promoter_bed_prep = PREPROCESS_PROMOTER( ch_promoter_bed )
    ch_enhancer_bed_prep = PREPROCESS_ENHANCER( ch_enhancer_bed )

    ch_pr_genseq = GENOMIC_SEQUENCE_PROMOTER(ch_promoter_bed_prep, ch_input_fasta)
    ch_enh_genseq = GENOMIC_SEQUENCE_ENHANCER(ch_enhancer_bed_prep, ch_input_fasta)
    emit: ch_enh_genseq.combine(ch_pr_genseq)
}

workflow{
    ch_genseq = prep_gen_seq()
    ch_genseq.view()
}
// nextflow  -C DeepTact-input-preprocessing/Genomic_seq.config -C promoter-enhancer-preprocesssing/nextflow.config run DeepTact-input-preprocessing/Genomic_seq_reg_elem.nf -profile slurm -w "/fastscratch/agarwa/nf-tmp/work" -with-timeline -resume 












