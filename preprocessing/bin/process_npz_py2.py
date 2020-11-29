# import colored_traceback.always
import sys
import h5py

import numpy as np


tmpBaseDir = "/projects/li-lab/agarwa/CUBE/DeepTact/dataset/DeepTact_tmp.1"
tmpBaseDir_hic = tmpBaseDir+'/'+"TrainingData"


def h5py_to_npz_py2(in_path, out_path):
    """
    convert hdf5 file to python2 readable .npz file
    """
    hdf_var = h5py.File(in_path, 'r')
    out_dict = {}
    for k in hdf_var.keys():
        out_dict[k] = hdf_var[k]
    np.savez(out_path, **out_dict)
    return 0
    


taskType = sys.argv[1]
cell = sys.argv[2]

assert taskType in ['prDNase', 'enhDNase' ,'prDNA','enhDNA']

if taskType == 'prDNA':
    inpath = "{_dir}/promoterDNA_{_cell}.hdf5".format(_dir=tmpBaseDir, _cell=cell)
    outpath = "{_dir}/{_cell}/P-E/promoter_Seq.npz".format(_dir=tmpBaseDir_hic, _cell=cell)
    h5py_to_npz_py2(inpath, outpath)
elif taskType == 'enhDNA':
    inpath = "{_dir}/enhancerDNA_{_cell}.hdf5".format(_dir=tmpBaseDir, _cell=cell)
    outpath = "{_dir}/{_cell}/P-E/enhancer_Seq.npz".format(_dir=tmpBaseDir_hic, _cell=cell) 
    h5py_to_npz_py2(inpath, outpath)
elif taskType == 'prDNase':
    inpath = "{_dir}/{_cell}/P-E/promoter_DNase.hdf5".format(_dir=tmpBaseDir_hic, _cell=cell)
    outpath = "{_dir}/{_cell}/P-E/promoter_DNase.npz".format(_dir=tmpBaseDir_hic, _cell=cell)
    h5py_to_npz_py2(inpath, outpath)
elif taskType == 'enhDNase':
    inpath = "{_dir}/{_cell}/P-E/enhancer_DNase.hdf5".format(_dir=tmpBaseDir_hic, _cell=cell)
    outpath = "{_dir}/{_cell}/P-E/enhancer_DNase.npz".format(_dir=tmpBaseDir_hic, _cell=cell)
    h5py_to_npz_py2(inpath, outpath)



# if taskType == 'prDNase':
#     file = h5py.File("{D}/pr_Seq.hdf5".format(D=base_train_py3), 'r')
#     np.savez("{D}/promoter_Seq.npz".format(D=base_train), **{'sequence':file['DNA'], 'label':file['label']})

# elif taskType == 'enhDNase':
#     file = h5py.File("{D}/enh_Seq.hdf5".format(D=base_train_py3), 'r')
#     np.savez("{D}/enhancer_Seq.npz".format(D=base_train), **{'sequence':file['DNA'], 'label':file['label']})


# elif taskType == 'prDNA':
#     file = h5py.File("{D}/pr_DNase.hdf5".format(D=base_train_py3), 'r')
#     np.savez("{D}/promoter_DNase.npz".format(D=base_train), **{'expr':file['DNase']})

# elif taskType == 'enhDNA':
#     file = h5py.File("{D}/enh_DNase.hdf5".format(D=base_train_py3), 'r')
#     np.savez("{D}/enhancer_DNase.npz".format(D=base_train), **{'expr':file['DNase']})

