
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
    memory '10 GB'

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

