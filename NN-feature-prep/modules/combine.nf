
// Combining Data : step 5.0
process SPLIT_HIC_AUG{
    input:
    tuple val(cellType), path(hic_aug)

    output:
    tuple val(cellType), path("hic.aug.${cellType}.*.csv")

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


// Combining Data : step 5.1.1
// combine PCHi-C interactions to get CO score for each element in the list
process COMBINE_PCHIC_CO_SCORE {
    label 'bigmem'
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
        h5f.create_dataset('data', data=A)
    os.system(f'gzip enhancer_hic.${cellType}.{num}.rep${rep}.h5')
    os.system(f'gzip promoter_hic.${cellType}.{num}.rep${rep}.h5')"
    """ 
}

// Combining Data : step 5.1.2a
// combine PCHi-C interactions to get CO score for each element in the list
process COMBINE_PCHIC_CO_SCORE_ENH {
    label 'bigmem'
    input:
    tuple val(cellType), val(rep), path(enhancer_rep_dat), path(promoter_rep_dat)


    output:
    tuple val(cellType), val(rep), path("enhancer_hic.${cellType}.rep${rep}.h5.gz")

    script:  
    """
    gzip -df $enhancer_rep_dat
    python -c "import numpy as np; 
    import h5py
    from natsort import natsorted
    L = natsorted('$enhancer_rep_dat'.split())
    L = [l[:-3] for l in L]
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

// Combining Data : step 5.1.2b
// combine PCHi-C interactions to get CO score for each element in the list
process COMBINE_PCHIC_CO_SCORE_PR {
    label 'bigmem'
    input:
    tuple val(cellType), val(rep), path(enhancer_rep_dat), path(promoter_rep_dat) 


    output:
    tuple val(cellType), val(rep), path("promoter_hic.${cellType}.rep${rep}.h5.gz") 

    script:  
    """
    gzip -df $promoter_rep_dat
    python -c "import numpy as np; 
    import h5py
    from natsort import natsorted
    L = natsorted('$promoter_rep_dat'.split())
    L = [l[:-3] for l in L]
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


// Combining Data : step 5.1.3
// combine PCHi-C interactions to get CO score for each element in the list
process COMBINE_CO_SCORE_REPS_ENH {
    label 'bigmem'

    input:
    tuple val(cellType), val(reps), path(enhancer_COscore_reps), path(promoter_COscore_reps) 

    output:
    // stdout result
    tuple val(cellType), path("enhancer_hic_COscore.${cellType}.h5.gz")

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

// Combining Data : step 5.1.3
// combine PCHi-C interactions to get CO score for each element in the list
process COMBINE_CO_SCORE_REPS_PR {
    label 'bigmem'

    input:
    tuple val(cellType), val(reps), path(enhancer_COscore_reps), path(promoter_COscore_reps) 

    output:
    // stdout result
    tuple val(cellType), path("promoter_hic_COscore.${cellType}.h5.gz") 

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
 
// Combining Data : step 5.2.1
// combine PCHi-C interactions to get DNA sequence for each element in the list
process COMBINE_PCHIC_DNA_SEQ {
    label 'bigCpuMem'

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
    os.system(f'rm enhancer_hic.${cellType}.{num}.DNA_seq.npz')
    os.system(f'rm promoter_hic.${cellType}.{num}.DNA_seq.npz')
    os.system(f'gzip enhancer_hic.${cellType}.{num}.DNA_seq.h5')
    os.system(f'gzip promoter_hic.${cellType}.{num}.DNA_seq.h5')
    "
    """ 
}

// Combining Data : step 5.2.2a
process COMBINE_PCHIC_OUT_ENHANCER {
    label 'bigmem'

    input:
    tuple val(cellType), path(enhancer_rep_dat), path(promoter_rep_dat) 

    output:
    tuple val(cellType), path("enhancer_hic.${cellType}.DNA_seq.h5.gz") 

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

// Combining Data : step 5.2.2b
process COMBINE_PCHIC_OUT_PROMOTER {
    label 'bigmem'

    input:
    tuple val(cellType), path(enhancer_rep_dat), path(promoter_rep_dat) 

    output:
    tuple val(cellType), path("promoter_hic.${cellType}.DNA_seq.h5.gz") 

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
process SAVE_ENH_DNA_SEQ {
    label 'bigmem'

    input:
    tuple val(cellType), path(enhancer_hic_DNAseq)
    output:
    tuple val(cellType), path("enhancer_hic.${cellType}.DNA_seq.final.h5.gz")

    script:  
    """
    gzip -df $enhancer_hic_DNAseq
    python -c "import numpy as np
    import h5py
    import math
    NUM_SEQ=4
    shape1 = (-1, 1, $params.enhancer_window, NUM_SEQ)
    with h5py.File('$enhancer_hic_DNAseq'[:-3],'r') as enhHicDnaseq_h5:
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

// Combining Data : step 5.3.1b
process SAVE_PR_DNA_SEQ {
    label 'bigmem'

    input:
    tuple val(cellType), path(promoter_hic_DNAseq)

    output:
    tuple val(cellType), path("promoter_hic.${cellType}.DNA_seq.final.h5.gz")

    script:  
    """
    gzip -df $promoter_hic_DNAseq
    python -c "import numpy as np
    import h5py
    import math
    NUM_SEQ=4
    shape2 = (-1, 1, $params.promoter_window, NUM_SEQ)
    with h5py.File('$promoter_hic_DNAseq'[:-3],'r') as prHicDnaseq_h5:
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

// Combining Data : step 5.3.1c
process SAVE_ENH_CO_SCORE {
    label 'bigmem'

    input:
    tuple val(cellType), path(enhancer_hic_COscore), path(promoter_hic_COscore)

    output:
    tuple val(cellType), path("enhancer_hic.${cellType}.CO_score.final.h5.gz") 

    script:  
    """
    gzip -df $enhancer_hic_COscore
    python -c "import numpy as np
    import h5py
    import math
    import os
    with h5py.File('$enhancer_hic_COscore'[:-3],'r') as enhHicCO_h5:
        Tregion1_expr = np.array(enhHicCO_h5['data'])
        NUM_REP = Tregion1_expr.shape[1]
        shape1 = (-1, 1, NUM_REP, $params.enhancer_window)
        Tregion1_expr.reshape(shape1)
        with h5py.File(f'enhancer_hic.${cellType}.CO_score.final.h5', 'w') as h5f:
            h5f.create_dataset('data', data=Tregion1_expr)
    os.system('rm $enhancer_hic_COscore'[:-3])
    "
    gzip enhancer_hic.${cellType}.CO_score.final.h5
    """
}

// Combining Data : step 5.3.1d
process SAVE_PR_CO_SCORE {
    label 'bigmem'

    input:
    tuple val(cellType), path(enhancer_hic_COscore), path(promoter_hic_COscore) 

    output:
    tuple val(cellType), path("promoter_hic.${cellType}.CO_score.final.h5.gz") 

    script:  
    """
    gzip -df $promoter_hic_COscore
    python -c "import numpy as np
    import h5py
    import math
    import os
    with h5py.File('$promoter_hic_COscore'[:-3],'r') as prHicCO_h5:
        Tregion2_expr = np.array(prHicCO_h5['data'])
        NUM_REP = Tregion2_expr.shape[1]
        shape2 = (-1, 1, NUM_REP, $params.promoter_window)
        Tregion2_expr.reshape(shape2)
        with h5py.File(f'promoter_hic.${cellType}.CO_score.final.h5', 'w') as h5f:
            h5f.create_dataset('data', data=Tregion2_expr)
    os.system('rm $promoter_hic_COscore'[:-3])
    "
    gzip promoter_hic.${cellType}.CO_score.final.h5
    """
}

// Combining Data : step 5.3.2
// separate the data in ch_hic_DNA_seq_features_out_, ch_hic_COscore_features_out into data points
process SEPARATE_DATA {
    label 'bigmem'

    input:
    tuple val(cellType), path(enhancer_hic_CO_gz), path(promoter_hic_CO_gz), path(enhancer_hic_DNAseq_gz), path(promoter_hic_DNAseq_gz), path(hic_aug) 

    output:
    tuple val(cellType), path("input_feature_part_*.h5.gz")

    script:  
    enhancer_hic_CO = enhancer_hic_CO_gz.baseName
    promoter_hic_CO = promoter_hic_CO_gz.baseName
    enhancer_hic_DNAseq = enhancer_hic_DNAseq_gz.baseName
    promoter_hic_DNAseq = promoter_hic_DNAseq_gz.baseName
    """
    gzip -df $enhancer_hic_CO_gz $promoter_hic_CO_gz $enhancer_hic_DNAseq_gz $promoter_hic_DNAseq_gz
    python -c "import pandas as pd
    import h5py
    import math
    import numpy as np

    Tlabel = pd.read_csv('$hic_aug')['label']

    h5f_enhDNAseq = h5py.File('$enhancer_hic_DNAseq','r')
    Tregion1_seq = np.expand_dims(np.array(h5f_enhDNAseq['data']),0)

    h5f_prDNAseq = h5py.File('$promoter_hic_DNAseq','r')
    Tregion2_seq = np.expand_dims(np.array(h5f_prDNAseq['data']),0)


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
    rm $enhancer_hic_CO $promoter_hic_CO $enhancer_hic_DNAseq $promoter_hic_DNAseq
    """ 
}

// ch_hic_features_split_out_ = ch_hic_features_split_out.take( params.dev ? params.dev_lim_tar : -1 )


// Combining Data : step 5.3.2
// separate the data in ch_hic_DNA_seq_features_out_, ch_hic_COscore_features_out into data points
process CONVERT_TAR_XZ {
    label 'bigCpuMem'
    errorStrategy 'finish'

    input:
    tuple val(cellType), path(input_feature_part_gz)

    output:
    tuple val(cellType), path("features.${cellType}.*.tar.xz") 

    script:  
    input_feature_part = input_feature_part_gz.baseName
    """
    gzip -df $input_feature_part_gz
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
    rm $input_feature_part
    """ 
}


// Combining Data : step 5.3.4
// separate the data in ch_hic_DNA_seq_features_out_, ch_hic_COscore_features_out into data points
process COMBINE_DATA_TAR {
    label 'bigmem'
    stageInMode 'copy'
    storeDir "${params.store_dir}"

    input:
    tuple val(cellType), path(input_xz)

    output:
    tuple val(cellType), path("features.${cellType}.tar.xz") 

    script:  
    """
    xzcat $input_xz | xz -c > features.${cellType}.tar.xz 
    """ 
    // uncompress command
    // tar --ignore-zeros -xJvf features.${cellType}.tar.xz
}
