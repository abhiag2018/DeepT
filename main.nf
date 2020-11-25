#! /projects/li-lab/agarwa/conda_envs/cube/bin/nextflow

ch_input = Channel.fromPath(params.input)

ch_input
  .splitCsv(header:true, sep:',')
  .map { row -> [ row.celltype, row.repetition, file(row.bam, checkIfExists: true) ]  }
  .set {ch_chrom_open_score}

process GEN_BAM_INDEX {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    tuple val(cell), val(rep), path(bamfile) from ch_chrom_open_score

    output:
    path '*.bam.bai' into result

    script:  // This script is bundled with the pipeline, in nf-core/atacseq/bin/
    """
      samtools index $bamfile ${cell}.${rep}.bam.bai
    """
}

process GEN_BAM_INDEX {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    tuple val(cell), val(rep), path(bamfile) from ch_chrom_open_score

    output:
    path '*.bam.bai' into result

    script:  // This script is bundled with the pipeline, in nf-core/atacseq/bin/
    """
      samtools index $bamfile ${cell}.${rep}.bam.bai
    """
}


// result.view { it.trim() }

//out_bed_ch
  //.collectFile(name: 'blast_output_combined.txt', storeDir: params.dataDir)