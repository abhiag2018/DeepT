nextflow.preview.dsl=2

if (params.species == "hg" ){
  all_chrom = "$params.all_chrom_hg"
  promoter_bedfile = "$params.promoter_bedfile_hg"
  enhancer_bedfile = "$params.enhancer_bedfile_hg"
  promoter_headers = "$params.promoter_headers_hg"
  enhancer_headers = "$params.enhancer_headers"
  insertRGB="True"
}
else if (params.species == "mm" ){
  all_chrom = "$params.all_chrom_mm"
  promoter_bedfile = "$params.promoter_bedfile_mm"
  enhancer_bedfile = "$params.enhancer_bedfile_mm"
  promoter_headers = "$params.promoter_headers_mm"
  enhancer_headers = "$params.enhancer_headers"
  insertRGB="False"
}
else{
    println "species (: $params.species) should be hg or mm; "
    exit 1
}

ch_promoter_bed = Channel.fromPath( promoter_bedfile )
ch_enhancer_bed = Channel.fromPath( enhancer_bedfile )

// promoter pre-processing
process PREPROCESS_PROMOTER {
    storeDir "${params.save_dir}"

    input:
    path(input)

    output:
    // stdout result
    tuple file("promoter_${params.species}.bed"), file("promoter_bg_${params.species}.bed")

    script:  
    """
    python -c "from preptools import process_promoter_bed; \
    process_promoter_bed('$input', \
        'promoter_${params.species}.bed', \
        $promoter_headers, \
        $all_chrom, \
        window = $params.promoter_window+$params.augment_length, \
        insertRGB=$insertRGB)"

    python -c "from preptools import process_promoter_bed; \
    process_promoter_bed('$input', \
        'promoter_bg_${params.species}.bed', \
        $promoter_headers, \
        $all_chrom, \
        window = $params.bgWindow, \
        insertRGB=$insertRGB)"
    """
    // func( allfield_bed, bg_path, headers, window=bg_window)
}

// enhancer pre-processing
process PREPROCESS_ENHANCER {
    storeDir "${params.save_dir}"

    input:
    path(input)

    output:
    // stdout result
    tuple file("enhancer_${params.species}.bed"), file("enhancer_bg_${params.species}.bed")


    script:  
    """
    python -c "from preptools import process_enhancer_bed; \
    process_enhancer_bed('$input', \
        'enhancer_${params.species}.bed', \
        $enhancer_headers, \
        window = $params.enhancer_window+$params.augment_length)"

    python -c "from preptools import process_enhancer_bed; \
    process_enhancer_bed('$input', \
        'enhancer_bg_${params.species}.bed', \
        $enhancer_headers, \
        window = $params.bgWindow)"
    """
}

workflow promoter_bed{
    ch_out = PREPROCESS_PROMOTER( ch_promoter_bed )
    emit: ch_out
}
workflow enhancer_bed{
    ch_out = PREPROCESS_ENHANCER( ch_enhancer_bed )
    emit: ch_out 
}

workflow {
    promoter_bed().view()
    enhancer_bed().view()
}