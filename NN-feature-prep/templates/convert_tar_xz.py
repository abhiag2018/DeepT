#!/usr/bin/env python
import os, sys
import multiprocessing
import h5py
import numpy as np

celltype=sys.argv[1]
input_feature_part=sys.argv[2]

h5f = h5py.File(input_feature_part,'r')

i=h5f['index'][()]
outdir=f'data_{i}'
os.makedirs(outdir,exist_ok=True)
# pool.starmap(npz_save, zip(map(lambda index:f'./{outdir}/{index}.npz', range(i,i+h5f['label'].shape[0])), h5f['enhseq'], h5f['prseq'], h5f['enhdnase'], h5f['prdnase'], h5f['label']))
for index in range(h5f['label'].shape[0]):
    np.savez(f'./{outdir}/{index}.npz', enh_seq = h5f['enhseq'][index,:,:,:] , 
        pr_seq = h5f['prseq'][index,:,:,:] , 
        enh_dnase = h5f['enhdnase'][index,:,:][np.newaxis,:,:] , 
        pr_dnase = h5f['prdnase'][index,:,:][np.newaxis,:,:] , 
        label = h5f['label'][index])

h5f.close()
print(i)
sys.stdout.flush()
sys.exit(0)
# os.system(f'tar cJvfh features.cellType.{i}.tar.xz data_{i}')