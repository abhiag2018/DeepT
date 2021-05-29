import os
import multiprocessing
import hickle as hkl
import h5py

import numpy as np

cellDir = "/projects/li-lab/agarwa/CUBE/DeepTact/dataset/DeepTact_tmp.1/TrainingData"
cell = 'nCD4'; cell_len = 685700; init = 0
# cell = 'tB'; cell_len = 670140; init = 685700
# cell = 'Mon'; cell_len = 591840; init = 670140+685700

def copy_links(idx, cell = cell, init=init, cellDir = cellDir):
    out = os.link(cellDir+'/'+cell+'/P-E/data/'+str(idx)+".npz", cellDir+'/COMB/P-E/data/'+str(init+idx)+".npz")
    return out


# copy data
# a_pool = multiprocessing.Pool(128)
# result = a_pool.map(copy_links, range(cell_len))


# ## copy hkl to h5py
# for cell in prep.DnaseCells.keys():
#     with h5py.File(f"{cellDir}/{cell}/P-E/test.hdf5",'w') as D:
#         D['idx'] = hkl.load(f"{cellDir}/{cell}/P-E/test.hkl")
#     with h5py.File(f"{cellDir}/{cell}/P-E/train.hdf5",'w') as D:
#         D['idx'] = hkl.load(f"{cellDir}/{cell}/P-E/train.hkl")

## train and test list
with h5py.File(f"{cellDir}/COMB/P-E/test.h5py", 'w') as D:
    D.create_dataset("idx", shape = (0,), dtype='i', maxshape=(None, ))

with h5py.File(f"{cellDir}/COMB/P-E/test.h5py", 'a') as D:
    sum_cell_len = 0
    for cell,cell_len in [('nCD4',685700),('tB',670140),('Mon',591840)]:
        cell = hkl.load(f"{cellDir}/{cell}/P-E/test.hkl")
        D['idx'].resize((D['idx'].shape[0]+cell.shape[0]), axis = 0)
        D['idx'][-cell.shape[0]:] = sum_cell_len+cell
        sum_cell_len +=  cell_len

with h5py.File(f"{cellDir}/COMB/P-E/test.h5py", 'r') as D:
    hkl.dump(np.array(D['idx'], np.int),f"{cellDir}/COMB/P-E/test.hkl")


with h5py.File(f"{cellDir}/COMB/P-E/train.h5py", 'w') as D:
    D.create_dataset("idx", shape = (0,), dtype='i', maxshape=(None, ))

with h5py.File(f"{cellDir}/COMB/P-E/train.h5py", 'a') as D:
    sum_cell_len = 0
    for cell,cell_len in [('nCD4',685700),('tB',670140),('Mon',591840)]:
        cell = hkl.load(f"{cellDir}/{cell}/P-E/train.hkl")
        D['idx'].resize((D['idx'].shape[0]+cell.shape[0]), axis = 0)
        D['idx'][-cell.shape[0]:] = sum_cell_len+cell
        sum_cell_len +=  cell_len

with h5py.File(f"{cellDir}/COMB/P-E/train.h5py", 'r') as D:
    hkl.dump(np.array(D['idx'], np.int),f"{cellDir}/COMB/P-E/train.hkl")

