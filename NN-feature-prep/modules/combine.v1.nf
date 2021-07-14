
// Combining Data : step 5.0
// split HiC csv files for parallel preprocessing
// output : hic.aug.${cellType}.{0-$hic_split_combine}.csv
// time : 1min
process SPLIT_HIC_AUG{
    tag "hic.aug.${cellType}.{split}.csv"
    memory '2 GB'

    input:
    tuple val(cellType), path(hic_aug)

    output:
    tuple val(cellType), path("hic.aug.${cellType}.*.csv")

    script:
    """
    num=`wc -l $hic_aug | awk '{print \$1}'`
    num_part=\$((num/7000))
    python -c "import preptools as pt; \
    pt.splitCSV('$hic_aug', \
        [], \
        readArgs = {}, \
        writeArgs = {'header':True}, prefix='hic.aug.${cellType}.', \
        split_num= max(min(\$num_part, $params.hic_split_combine),1), \
        suffix='.csv')"
    """
}


// Combining Data : step 5.1.1
// combine PCHi-C interactions to get CO score for each element in the hic list; and for each technical repetition .bam 
// output : enhancer_hic.${cellType}.{0-19}.rep${rep}.h5.gz
//          promoter_hic.${cellType}.{0-19}.rep${rep}.h5.gz
// time : 1h20m
process COMBINE_PCHIC_CO_SCORE {
    tag "promoter/enhancer_hic.${cellType}.{split}.rep${rep}.h5.gz"
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }
    // afterScript "rm *.npz"

    input:
    tuple val(cellType), val(rep), path(enhancer_COscore), path(enhancer_bed), path(enhancer_bg), path(promoter_COscore), path(promoter_bed), path(promoter_bg), path(hic_aug)

    output:
    tuple val(cellType), val(rep), path("enhancer_hic.${cellType}.*.rep${rep}.h5.gz"), path("promoter_hic.${cellType}.*.rep${rep}.h5.gz")

    script:  
    """
    python -c "import preptools as pt; \
    import pandas as pd; \
    import numpy as np; \
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
    import os; \
    num='$hic_aug'.split('.')[-2]; \
    A = np.load(f'promoter_hic.${cellType}.{num}.rep${rep}.npz', allow_pickle=True)['expr']
    with h5py.File(f'promoter_hic.${cellType}.{num}.rep${rep}.h5', 'w') as h5f:
        h5f.create_dataset('data', data=A)"
    num=`echo $hic_aug | awk '{split(\$0,a,"."); print a[4]}'`
    gzip enhancer_hic.${cellType}.\${num}.rep${rep}.h5
    gzip promoter_hic.${cellType}.\${num}.rep${rep}.h5
    """ 
}

// Combining Data : step 5.1.2a
// combine the CO Score features for hic list parts for enhancers
// output : enhancer_hic.combined.${cellType}.rep${rep}.h5.gz
// time : 42m
process COMBINE_PCHIC_CO_SCORE_ENH {
    tag "enhancer_hic.combined.${cellType}.rep${rep}.h5.gz"
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }

    input:
    tuple val(cellType), val(rep), path(enhancer_rep_dat), path(promoter_rep_dat)


    output:
    tuple val(cellType), val(rep), path("enhancer_hic.combined.${cellType}.rep${rep}.h5.gz")

    script:  
    """
    gzip -df $enhancer_rep_dat
    python -c "import numpy as np; 
    import h5py
    from natsort import natsorted
    L = natsorted('$enhancer_rep_dat'.split())
    L = [l[:-3] for l in L]
    h5f_list = list(map(lambda file_name:h5py.File(file_name,'r'), L) )
    h5f = h5py.File('enhancer_hic.combined.${cellType}.rep${rep}.h5','w')
    h5f.create_dataset('data',data=np.zeros((np.sum(list(map(lambda i:h5f_list[i]['data'].shape[0],range(len(L))))),h5f_list[0]['data'].shape[1])))
    row=0;
    for i in range(len(L)):
        s = h5f_list[i]['data'].shape
        h5f['data'][row:row+s[0],:]=h5f_list[i]['data']
        row += s[0]
    list(map(lambda f:f.close(), h5f_list))
    h5f.close()"
    old_files="$enhancer_rep_dat.baseName"
    if [ "\${old_files:0:1}" == "[" ]; then
        NM="\${old_files:1:\${#old_files}-2}"
        rm `echo "\${NM//, / }"`
    else
        rm \$old_files
    fi
    gzip enhancer_hic.combined.${cellType}.rep${rep}.h5
    """ 
}

// Combining Data : step 5.1.2b
// combine the CO Score features for hic list parts for promoter
// output : promoter_hic.combined.${cellType}.rep${rep}.h5.gz
// time : 22m
process COMBINE_PCHIC_CO_SCORE_PR {
    tag "promoter_hic.combined.${cellType}.rep${rep}.h5.gz"
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }


    input:
    tuple val(cellType), val(rep), path(enhancer_rep_dat), path(promoter_rep_dat) 


    output:
    tuple val(cellType), val(rep), path("promoter_hic.combined.${cellType}.rep${rep}.h5.gz") 

    script:  
    """
    gzip -df $promoter_rep_dat
    python -c "import numpy as np; 
    import h5py
    from natsort import natsorted
    L = natsorted('$promoter_rep_dat'.split())
    L = [l[:-3] for l in L]
    h5f_list = list(map(lambda file_name:h5py.File(file_name,'r'), L) )
    h5f = h5py.File('promoter_hic.combined.${cellType}.rep${rep}.h5','w')
    h5f.create_dataset('data',data=np.zeros((np.sum(list(map(lambda i:h5f_list[i]['data'].shape[0],range(len(L))))),h5f_list[0]['data'].shape[1])))
    row=0;
    for i in range(len(L)):
        s = h5f_list[i]['data'].shape
        h5f['data'][row:row+s[0],:]=h5f_list[i]['data']
        row += s[0]
    list(map(lambda f:f.close(), h5f_list))
    h5f.close()"
    old_files="$promoter_rep_dat.baseName"
    if [ "\${old_files:0:1}" == "[" ]; then
        NM="\${old_files:1:\${#old_files}-2}"
        rm `echo "\${NM//, / }"`
    else
        rm \$old_files
    fi
    gzip promoter_hic.combined.${cellType}.rep${rep}.h5
    """ 
}


// Combining Data : step 5.1.3a
// combine the CO Score features for the combined hic list for enhancers, 
//      across different .bam repetitions
// output : enhancer_hic.combined.${cellType}.COscore.h5.gz
// time : 13h18m
process COMBINE_CO_SCORE_REPS_ENH {
    tag "enhancer_hic.combined.${cellType}.COscore.h5.gz"

    input:
    tuple val(cellType), val(reps), path(enhancer_COscore_reps), path(promoter_COscore_reps) 

    output:
    // stdout result
    tuple val(cellType), path("enhancer_hic.combined.${cellType}.COscore.h5.gz")

    script:  
    """
    gzip -df $enhancer_COscore_reps
    python -c "from natsort import natsorted
    import numpy as np
    import h5py
    import math
    L = natsorted('$enhancer_COscore_reps'.split())
    L = [l[:-3] for l in L]
    dnase_cell = list(map(lambda file_name:h5py.File(file_name,'r'), L))
    h5f = h5py.File('enhancer_hic.combined.${cellType}.COscore.h5','w')
    h5f.create_dataset('data',data=np.zeros((dnase_cell[0]['data'].shape[0], len(L), dnase_cell[0]['data'].shape[1])))
    load_rows=10000
    for j in range(len(L)):
        for i in range(math.ceil(dnase_cell[0]['data'].shape[0]/load_rows)):
            h5f['data'][load_rows*i:load_rows*(i+1),j,:] = dnase_cell[j]['data'][load_rows*i:load_rows*(i+1),:]
    h5f.close()
    for i in range(len(dnase_cell)):
        dnase_cell[i].close()
    "
    old_files="$enhancer_COscore_reps.baseName"
    if [ "\${old_files:0:1}" == "[" ]; then
        NM="\${old_files:1:\${#old_files}-2}"
        rm `echo "\${NM//, / }"`
    else
        rm \$old_files
    fi
    gzip enhancer_hic.combined.${cellType}.COscore.h5
    """ 
}


// Combining Data : step 5.1.3b
// combine the CO Score features for the combined hic list for promoters, 
//      across different .bam repetitions
// output : promoter_hic.combined.${cellType}.COscore.h5.gz
// time : 10h18m
process COMBINE_CO_SCORE_REPS_PR {
    tag "promoter_hic.combined.${cellType}.COscore.h5.gz"

    input:
    tuple val(cellType), val(reps), path(enhancer_COscore_reps), path(promoter_COscore_reps) 

    output:
    // stdout result
    tuple val(cellType), path("promoter_hic.combined.${cellType}.COscore.h5.gz") 

    script:  
    """ 
    gzip -df $promoter_COscore_reps
    python -c "from natsort import natsorted
    import numpy as np
    import h5py
    import math
    L = natsorted('$promoter_COscore_reps'.split())
    L = [l[:-3] for l in L]
    dnase_cell = list(map(lambda file_name:h5py.File(file_name,'r'), L))
    h5f = h5py.File('promoter_hic.combined.${cellType}.COscore.h5','w')
    h5f.create_dataset('data',data=np.zeros((dnase_cell[0]['data'].shape[0], len(L), dnase_cell[0]['data'].shape[1])))
    load_rows=10000
    for j in range(len(L)):
        for i in range(math.ceil(dnase_cell[0]['data'].shape[0]/load_rows)):
            h5f['data'][load_rows*i:load_rows*(i+1),j,:] = dnase_cell[j]['data'][load_rows*i:load_rows*(i+1),:]
    h5f.close()
    for i in range(len(dnase_cell)):
        dnase_cell[i].close()
    "
    old_files="$promoter_COscore_reps.baseName"
    if [ "\${old_files:0:1}" == "[" ]; then
        NM="\${old_files:1:\${#old_files}-2}"
        rm `echo "\${NM//, / }"`
    else
        rm \$old_files
    fi
    gzip promoter_hic.combined.${cellType}.COscore.h5
    """ 
}
 
// Combining Data : step 5.2.1
// combine the DNA sequence for combined hic list parts for promoter & enhancer
// output : enhancer_hic.${cellType}.{0-19}.DNA_seq.h5.gz
//          promoter_hic.${cellType}.{0-19}.DNA_seq.h5.gz
// time : 1h32m
process COMBINE_PCHIC_DNA_SEQ {
    tag "promoter/enhancer_hic.${cellType}.{split}.DNA_seq.h5.gz"
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }

    input:
    tuple val(cellType), path(hic_aug), path(enhancer_DNAseq), path(promoter_DNAseq), path(enhancer_bed), path(enhancer_bg), path(promoter_bed), path(promoter_bg)

    output:
    tuple val(cellType), path("enhancer_hic.${cellType}.*.DNA_seq.h5.gz"), path("promoter_hic.${cellType}.*.DNA_seq.h5.gz") 

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
    "
    num=`echo $hic_aug | awk '{split(\$0,a,"."); print a[4]}'`
    rm enhancer_hic.${cellType}.\${num}.DNA_seq.npz
    rm promoter_hic.${cellType}.\${num}.DNA_seq.npz
    gzip enhancer_hic.${cellType}.\${num}.DNA_seq.h5
    gzip promoter_hic.${cellType}.\${num}.DNA_seq.h5
    """ 
}

// Combining Data : step 5.2.2a
// combine the DNA sequence across the hic list parts for enhancers
// output : enhancer_hic.combined.${cellType}.DNA_seq.h5.gz
// time : 6h7m
process COMBINE_PCHIC_OUT_ENHANCER {
    tag "enhancer_hic.combined.${cellType}.DNA_seq.h5.gz"
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }


    input:
    tuple val(cellType), path(enhancer_rep_dat), path(promoter_rep_dat) 

    output:
    tuple val(cellType), path("enhancer_hic.combined.${cellType}.DNA_seq.h5.gz") 

    script:
    """
    gzip -df $enhancer_rep_dat
    python -c "import numpy as np; 
    import h5py
    from natsort import natsorted
    L=natsorted('${enhancer_rep_dat}'.split())
    L = [l[:-3] for l in L]
    h5f_list = list(map(lambda file_name:h5py.File(file_name,'r'), L) )
    arr_list = list(map(lambda f:f['data'], h5f_list))
    h5f = h5py.File('enhancer_hic.combined.${cellType}.DNA_seq.h5','w')
    h5f.create_dataset('data',data=np.zeros((np.sum(list(map(lambda i:h5f_list[i]['data'].shape[0],range(len(L))))),h5f_list[0]['data'].shape[1])))
    row=0;
    for i in range(len(L)):
        s = h5f_list[i]['data'].shape
        h5f['data'][row:row+s[0],:]=h5f_list[i]['data']
        row += s[0]
    list(map(lambda f:f.close(), h5f_list))
    h5f.close()"
    old_files="$enhancer_rep_dat.baseName"
    if [ "\${old_files:0:1}" == "[" ]; then
        NM="\${old_files:1:\${#old_files}-2}"
        rm `echo "\${NM//, / }"`
    else
        rm \$old_files
    fi
    gzip enhancer_hic.combined.${cellType}.DNA_seq.h5
    """
}


// Combining Data : step 5.2.2b
// combine the DNA sequence across the hic list parts for promoters
// output : promoter_hic.combined.${cellType}.DNA_seq.h5.gz
// time : 3h
process COMBINE_PCHIC_OUT_PROMOTER {
    tag "promoter_hic.combined.${cellType}.DNA_seq.h5.gz"
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }


    input:
    tuple val(cellType), path(enhancer_rep_dat), path(promoter_rep_dat) 

    output:
    tuple val(cellType), path("promoter_hic.combined.${cellType}.DNA_seq.h5.gz") 

    script:
    """
    gzip -df $promoter_rep_dat
    python -c "import numpy as np; 
    import h5py
    from natsort import natsorted
    L=natsorted('${promoter_rep_dat}'.split())
    L = [l[:-3] for l in L]
    h5f_list = list(map(lambda file_name:h5py.File(file_name,'r'), L) )
    arr_list = list(map(lambda f:f['data'], h5f_list))
    h5f = h5py.File('promoter_hic.combined.${cellType}.DNA_seq.h5','w')
    h5f.create_dataset('data',data=np.zeros((np.sum(list(map(lambda i:h5f_list[i]['data'].shape[0],range(len(L))))),h5f_list[0]['data'].shape[1])))
    row=0;
    for i in range(len(L)):
        s = h5f_list[i]['data'].shape
        h5f['data'][row:row+s[0],:]=h5f_list[i]['data']
        row += s[0]
    list(map(lambda f:f.close(), h5f_list))
    h5f.close()"
    old_files="$promoter_rep_dat.baseName"
    if [ "\${old_files:0:1}" == "[" ]; then
        NM="\${old_files:1:\${#old_files}-2}"
        rm `echo "\${NM//, / }"`
    else
        rm \$old_files
    fi
    gzip promoter_hic.combined.${cellType}.DNA_seq.h5
    """
}


// Combining Data : developer process for unzipping files
process UNZIP {

    input:
    tuple val(cellType), path(file_gz)

    output:
    tuple val(cellType), path(filename)

    script:  
    filename = file_gz.baseName
    """
    gzip -df $file_gz
    """
}

// Combining Data : step 5.3
// separate the data into NPZ files.. into $sepdata_split chunks for easy processing in next step
// output : data_chr.${cellType}.tar.gz
// time: 1h40m
process SEPARATE_DATA {
    tag "data_chr.${cellType}.tar.gz"
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
    import multiprocessing as mp
    import numpy as np

    hicInteractions = pd.read_csv('$hic_aug')
    Tlabel = hicInteractions['label']
    Tindex = hicInteractions['index']
    Tenh_coscore = h5py.File('$enhancer_hic_CO_gz.baseName','r')['data']
    Tenh_dnaseq = h5py.File('$enhancer_hic_DNAseq_gz.baseName','r')['data']
    Tpr_coscore = h5py.File('$promoter_hic_CO_gz.baseName','r')['data']
    Tpr_dnaseq = h5py.File('$promoter_hic_DNAseq_gz.baseName','r')['data']

    def npz_save(enh_coscore_, enh_dnaseq_, pr_coscore_, pr_dnaseq_, label_, index_): 
        print(f'data/{index_}.npz',flush=True)
        np.savez(f'data/{index_}.npz', enh_seq = enh_dnaseq_.reshape((1,$params.enhancer_window,4)).transpose(0,2,1) , 
            pr_seq = pr_dnaseq_.reshape((1,$params.promoter_window,4)).transpose(0,2,1) , 
            enh_dnase = enh_coscore_[np.newaxis,:,:] , 
            pr_dnase = pr_coscore_[np.newaxis,:,:] , 
            label = label_,
            index = index_)

    os.makedirs('data')
    print(f'{mp.cpu_count()} CPUs visible', flush=True)
    with mp.Pool(processes=None) as pool:
        pool.starmap(npz_save, zip(Tenh_coscore, Tenh_dnaseq, Tpr_coscore, Tpr_dnaseq, Tlabel, Tindex ))
    "
    # rm $enhancer_hic_CO_gz.baseName $promoter_hic_CO_gz.baseName $enhancer_hic_DNAseq_gz.baseName $promoter_hic_DNAseq_gz.baseName
    tar -czvf data_chr.${cellType}.tar.gz data --remove-files
    """
}
// ch_hic_features_split_out_ = ch_hic_features_split_out.take( params.dev ? params.dev_lim_tar : -1 )

