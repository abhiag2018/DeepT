nextflow.preview.dsl=2

if (params.species == "hg" ){
    ch_input_fasta = Channel.fromPath(params.human_genome_fasta)
}
else if (params.species == "mm" ){
    ch_input_fasta = Channel.fromPath(params.mouse_genome_fasta)
}
else{
    println "species (: $params.species) should be hg or mm; "
    exit 1
}

ch_enhancer_bed_prep = Channel.fromPath("$params.store_dir/enhancer_${params.species}.bed").combine(Channel.fromPath("$params.store_dir/enhancer_bg_${params.species}.bed"))
ch_promoter_bed_prep = Channel.fromPath("$params.store_dir/promoter_${params.species}.bed").combine(Channel.fromPath("$params.store_dir/promoter_bg_${params.species}.bed"))

// promoter genomic sequence computation 
process GENOMIC_SEQUENCE_PROMOTER {
    storeDir "${params.save_dir}"

    input:
    tuple path(promoter), path(promoter_bg)
    path(fasta)

    output:
    file("promoter_DNAseq.${fasta.baseName}.npz.gz")

    script:  
    """
    python -c "import preptools as pt; \
    pt.generate_elem_fa('$fasta', \
        '$promoter', \
        out_fa = 'promoter.fa' )"

    python -c "import preptools as pt; \
    pt.fasta_to_onehot('promoter.fa', \
        outp='promoter_DNAseq.${fasta.baseName}.npz')"
    gzip promoter_DNAseq.${fasta.baseName}.npz
    """
}



// enhancer genomic sequence computation 
process GENOMIC_SEQUENCE_ENHANCER {
    storeDir "${params.save_dir}"

    input:
    tuple path(enhancer), path(enhancer_bg)
    path(fasta)

    output:
    file("enhancer_DNAseq.${fasta.baseName}.npz.gz")

    script:  
    """
    python -c "import preptools as pt; \
    pt.generate_elem_fa('$fasta', \
        '$enhancer', \
        out_fa = 'enhancer.fa' )"

    python -c "import preptools as pt; \
    pt.fasta_to_onehot('enhancer.fa', \
        outp='enhancer_DNAseq.${fasta.baseName}.npz')"
    gzip enhancer_DNAseq.${fasta.baseName}.npz
    """
}



workflow prep_gen_seq{
    ch_pr_genseq = GENOMIC_SEQUENCE_PROMOTER(ch_promoter_bed_prep, ch_input_fasta)
    ch_enh_genseq = GENOMIC_SEQUENCE_ENHANCER(ch_enhancer_bed_prep, ch_input_fasta)
    emit: ch_enh_genseq.combine(ch_pr_genseq)
}

workflow{
    ch_genseq = prep_gen_seq()
    ch_genseq.view()
}
// nextflow  -C DeepTact-input-preprocessing/Genomic_seq.config -C promoter-enhancer-preprocesssing/nextflow.config run DeepTact-input-preprocessing/Genomic_seq_reg_elem.nf -profile slurm -w "/fastscratch/agarwa/nf-tmp/work" -with-timeline -resume 












