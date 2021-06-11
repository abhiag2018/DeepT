#!/usr/bin/env python
import sys
import pandas as pd
import h5py
import math
import numpy as np

hic_aug = sys.argv[1]
enhancer_hic_DNAseq = sys.argv[2]
promoter_hic_DNAseq = sys.argv[3]
enhancer_hic_CO = sys.argv[4]
promoter_hic_CO = sys.argv[5]
sepdata_split = int(sys.argv[6])
# print(enhancer_hic_DNAseq, promoter_hic_DNAseq, enhancer_hic_CO, promoter_hic_CO, sepdata_split)

Tlabel = pd.read_csv(hic_aug)['label']

print(enhancer_hic_DNAseq)
h5f_enhDNAseq = h5py.File(enhancer_hic_DNAseq,'r')
Tregion1_seq = h5f_enhDNAseq['data']

h5f_prDNAseq = h5py.File(promoter_hic_DNAseq,'r')
Tregion2_seq = h5f_prDNAseq['data']


## load data: DNase
h5f_enhDNase = h5py.File(enhancer_hic_CO,'r')
Tregion1_expr = h5f_enhDNase['data']

h5f_prDNase = h5py.File(promoter_hic_CO,'r')
Tregion2_expr = h5f_prDNase['data']


NUM = Tlabel.shape[0]

step=math.ceil(NUM/sepdata_split)

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
h5f_prDNase.close()