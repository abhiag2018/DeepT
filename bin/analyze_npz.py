import colored_traceback.always
import sys, os, argparse
import swifter
import itertools
import pickle
import h5py

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

import scipy.stats as sps
import numpy as np
import pandas as pd


import preptools as pt
import preprocessing as prep


args = pt.process_inputArgs(input_parse=sys.argv[1:], argType={'file_index':None, 'nTasks':None, 'taskType':None , 'cellType': None } )


assert args.taskType in ['prDNA','enhDNA','stats','plotDNase','tohdf5']

base_train = "/projects/li-lab/agarwa/CUBE/DeepTact/dataset/DeepTact_tmp.1/train_Mon/P-E"
base_train_py3 = "/projects/li-lab/agarwa/CUBE/DeepTact/dataset/DeepTact_tmp.1/train_Mon/P-E/python3"


def onehotDNA_as_pdseries(hdfpath):
    path_list = hdfpath.split("/")
    out_path = '/'.join(path_list[:-1])+'/'+path_list[-1].split('.')[0]+'_pd.hdf5'

    def vectoATCG(vec):
        if sum(vec)>1:
            raise NameError('multiple nucleotides!')
        elif vec[0]==1:
            return ('A')
        elif vec[1]==1:
            return ('T')
        elif vec[2]==1:
            return ('C')
        elif vec[3]==1:
            return ('G')
        return 'N'
    convert_to_seq = lambda onehot_arr: pd.Series(list(onehot_arr)).apply(vectoATCG).str.cat()

    onehot_DNA_pd =  pd.Series(list(h5py.File(hdfpath, 'r')['sequence']))
    onehot_DNA_pd_arr = onehot_DNA_pd.swifter.apply(lambda x:np.array(x).reshape(x.shape[0]//4,4))
    DNA_seq = onehot_DNA_pd_arr.swifter.apply(convert_to_seq)
    
    store = pd.HDFStore(out_path)
    store['sequence'] = DNA_seq
    store.close()
    # pd.read_hdf(out_path)
    return 0


def get_nucle_count(hdfpath, sample_n_frac = None):
    path_list = hdfpath.split("/")
    out_path = '/'.join(path_list[:-1])+'/'+path_list[-1]+'_count.hdf5'

    one_hot_seq_pd = pd.Series(h5py.File(hdfpath, 'r')['sequence'])
    if sample_n_frac:
        if type(sample_n_frac) == int:
            one_hot_seq_pd = one_hot_seq_pd.sample(n=sample_n_frac)
        else:
            one_hot_seq_pd = one_hot_seq_pd.sample(frac=sample_n_frac)

    transform = lambda x:x.reshape(x.shape[-1]//4,4)

    onehot_reshaped = one_hot_seq_pd.apply(lambda x:transform(np.array(x)))

    first = lambda x: x[0] if len(x)>0 else 'N'
    def vectoATCG(vec):
        if sum(vec)>1:
            raise NameError('multiple nucleotides!')
        return first(list(itertools.compress('ATCG', [x==1 for x in vec])))
    convert_to_seq = lambda onehot_arr: pd.Series(list(onehot_arr)).apply(vectoATCG).str.cat()

    onehot_seq = onehot_reshaped.apply(convert_to_seq)
    store = pd.HDFStore(out_path)
    store['A'] = onehot_seq.apply(lambda x:x.count('A'))
    store['T'] = onehot_seq.apply(lambda x:x.count('T'))
    store['C'] = onehot_seq.apply(lambda x:x.count('C'))
    store['G'] = onehot_seq.apply(lambda x:x.count('G'))
    
    # pd.read_hdf(out_pathA)
    return (onehot_seq.apply(lambda x:x.count('A')),onehot_seq.apply(lambda x:x.count('T')),
         onehot_seq.apply(lambda x:x.count('C')),onehot_seq.apply(lambda x:x.count('G')))

def plot_DNase(cell,bamf):
    label_cell = pd.read_csv(f"{prep.tmpBaseDir}/TrainingData/{cell}/P-E/pairs_train_augment.csv").label
    # bamfs_cell = list(itertools.chain.from_iterable(prep.DnaseCells[cell]))
    # for bamf in bamfs_cell:
    bamf_dir = prep.dnaseTmpDir(bamf)
    print(f"{cell}, {os.path.basename(bamf_dir)}..",end="")

    enh_DNase = np.nan_to_num(np.stack(np.load(f"{bamf_dir}/enhancer_Train_Augment.npz", allow_pickle=True)['expr']))
    pr_DNase = np.nan_to_num(np.stack(np.load(f"{bamf_dir}/promoter_Train_Augment.npz", allow_pickle=True)['expr']))

    enh_DNase_agg = np.sum(enh_DNase,axis=1) 
    pr_DNase_agg = np.sum(pr_DNase,axis=1)

    # fs=15
    # fig, axes = plt.subplots(1,1,figsize=(5,5))

    # axes.scatter(enh_DNase_agg[label_cell==1],pr_DNase_agg[label_cell==1],color='blue',alpha=1.0,sizes=(10,))
    # axes.scatter(enh_DNase_agg[label_cell==0],pr_DNase_agg[label_cell==0],color='red',alpha=.25,sizes=(5,),marker='o', facecolors='none')

    # axes.set_ylabel('promoter score',fontsize=fs)
    # axes.set_xlabel('enhancer score',fontsize=fs)
    # fig.suptitle('Chromatin Openness Score',fontsize=fs)
    # plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    # custom_lines = (matplotlib.lines.Line2D([0], [0], color='blue', lw=8,alpha=1), matplotlib.lines.Line2D([0], [0], color='red', lw=8,alpha=1))
    # vpt = fig.legend(custom_lines, ('Positive Label', 'Neg Label'),bbox_to_anchor=(.95,.85),fontsize = fs)

    # plt.show()
    # plt.savefig(f'{prep.tmpBaseDir}/cell_{cell}_{os.path.basename(bamf_dir)}.png')
    # plt.close('all')

    pos_data_x = enh_DNase_agg[label_cell==1]
    neg_data_x = enh_DNase_agg[label_cell==0]
    pos_data_y = pr_DNase_agg[label_cell==1]
    neg_data_y = pr_DNase_agg[label_cell==0]

    ## calculate axes limits
    pos_data = np.stack([pos_data_x, pos_data_y])
    neg_data = np.stack([neg_data_x, neg_data_y])
    kde = sps.gaussian_kde(pos_data)# estimate kernel density of data

    xx, yy = np.meshgrid(
        np.linspace(-25000, 200000, 500),
        np.linspace(-25000, 200000, 500))

    z = kde.pdf([xx.ravel(), yy.ravel()]).reshape(xx.shape)# calculate probability density on these points

    # proportion of points above the 10%, i.e. approx the third contour line
    zi = z > np.max(z) * 0.1


    xxx = xx[0][np.logical_and(xx[0]>min(xx[zi]), xx[0]<max(xx[zi]))]; xlen = len(xxx)
    yyy = yy[:,0][np.logical_and(yy[:,0]>min(yy[zi]), yy[:,0]<max(yy[zi]))]; ylen = len(yyy)
    xxx = np.stack([xxx]*ylen)
    yyy = np.stack([yyy]*xlen).T

    z = sps.gaussian_kde(pos_data).pdf([xxx.ravel(), yyy.ravel()]).reshape(xxx.shape)
    plt.pcolormesh(xxx, yyy, z, cmap='Reds')
    plt.xlabel('enhancer')
    plt.ylabel('promoter')
    plt.title('Pos interactions')
    plt.savefig(f'{prep.tmpBaseDir}/cell_{cell}_pos_{os.path.basename(bamf_dir)}.png')

    z = sps.gaussian_kde(neg_data).pdf([xxx.ravel(), yyy.ravel()]).reshape(xxx.shape)
    plt.pcolormesh(xxx, yyy, z, cmap='Reds')
    plt.xlabel('enhancer')
    plt.ylabel('promoter')
    plt.title('Neg interactions')
    plt.savefig(f'{prep.tmpBaseDir}/cell_{cell}_neg_{os.path.basename(bamf_dir)}.png')
    
    # g = sns.JointGrid(x=pos_data_x, y=pos_data_y, hue=[1]*len(pos_data_x), xlim = [min(xx[zi]), max(xx[zi])],ylim = [min(yy[zi]), max(yy[zi])])
    # g.plot_joint(sns.kdeplot, cmap='Reds', shade=True, alpha=0.7 )
    # g.plot_marginals(sns.kdeplot, color='red', shade=True)
    # plt.savefig(f'{prep.tmpBaseDir}/cell_{cell}_pos_{os.path.basename(bamf_dir)}.png')

    # g = sns.JointGrid(x=neg_data_x, y=neg_data_y, hue=[1]*len(neg_data_x), xlim = [min(xx[zi]), max(xx[zi])],ylim = [min(yy[zi]), max(yy[zi])])
    # g.plot_joint(sns.kdeplot, cmap='Reds', shade=True, alpha=0.7 )
    # g.plot_marginals(sns.kdeplot, color='red', shade=True)
    # plt.savefig(f'{prep.tmpBaseDir}/cell_{cell}_neg_{os.path.basename(bamf_dir)}.png')

    plt.close('all')


def plot_DNA():
    prDNA_hdf_pd = OrderedDict()
    enhDNA_hdf_pd = OrderedDict()
    prDNA_count = OrderedDict()
    enhDNA_count = OrderedDict()
    counts_nuc = list('ATCGN')
    for cell in cell_list:
        prDNA_hdfpath = f"{prep.tmpBaseDir}/promoterDNA_{cell}_pd.hdf5"
        enhDNA_hdfpath = f"{prep.tmpBaseDir}/enhancerDNA_{cell}_pd.hdf5"
        prDNA_hdf_pd[cell] = pd.read_hdf(prDNA_hdfpath)
        enhDNA_hdf_pd[cell] = pd.read_hdf(enhDNA_hdfpath)
        prDNA_count[cell] = OrderedDict()
        enhDNA_count[cell] = OrderedDict()
        for nuc in counts_nuc:
            prDNA_count[cell][nuc] = prDNA_hdf_pd[cell].swifter.apply(lambda x:x.count(nuc))
            enhDNA_count[cell][nuc] = enhDNA_hdf_pd[cell].swifter.apply(lambda x:x.count(nuc))
        enhDNA_count[cell] = pd.DataFrame(enhDNA_count[cell])


    violinplot_options = {'showmeans' : True, 'showmedians'  :  True, 'quantiles'  :  [[0,.25,.5,.75,1.0],[0,.25,.5,.75,1.0]]}


    cell_list = list(prep.DnaseCells.keys())
    _flatten = lambda x:list(itertools.chain.from_iterable(x))
    fs=40
    custom_xlim = (0.5, 2.5)
    custom_ylim = (0, 1200)
    ncols = 5
    for cell in list(cell_list):
        label_cell = pd.read_csv(f"{prep.tmpBaseDir}/TrainingData/{cell}/P-E/pairs_train_augment.csv").label
        fig, axes = plt.subplots(2,ncols-2,figsize = (ncols*10,(ncols-2)*10))

        axes= _flatten(axes)
        vpt = fig.suptitle(f'cell: {cell}', fontsize=fs*2)
        print(cell,end="")
        for ax, nuc in zip((axes),counts_nuc):
            print(".",end="",flush=True)
            data = [list(prDNA_count[cell][nuc][label_cell==1]), list(enhDNA_count[cell][nuc][label_cell==1])]
            vp1 = ax.violinplot(data, **violinplot_options)
            for pc in vp1['bodies']:
                c1 = pc.get_facecolor()
                vpt = pc.set_edgecolor('black')
                vpt = pc.set_alpha(0.7)
            data = [list(prDNA_count[cell][nuc][label_cell==0]), list(enhDNA_count[cell][nuc][label_cell==0])]
            vp2 = ax.violinplot(data, **violinplot_options)
            for pc in vp2['bodies']:
                c2 = pc.get_facecolor()
                vpt = pc.set_edgecolor('black')
                vpt = pc.set_alpha(0.7)

            labels = ['promoters','enhancers']
            vpt = ax.tick_params(
                axis='y',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=True,      # ticks along the bottom edge are off
                top=True,         # ticks along the top edge are off
                labelbottom=False,
                length = 10,
                width=2,
                direction='out',
                labelsize = fs)
            vpt = ax.set(xlim=custom_xlim,ylim=custom_ylim, xticks= np.arange(1, len(labels) + 1),yticks=np.arange(0,1200,100))
    #         vpt = ax.set_xticks(np.arange(1, len(labels) + 1))
            vpt = ax.set_xticklabels(labels,fontsize=fs)
    #         vpt = ax.set_yticks(np.arange(0,1200,100),minor=True)
            vpt = ax.set_yticklabels(np.arange(0,1200,100),fontsize=fs)
            vpt = ax.set_ylabel('nucleotide count',fontsize=fs)
            vpt = ax.set_title(f'nucleotude {nuc}', fontsize=fs)  


        custom_lines = (Line2D([0], [0], color=c1[0], lw=8,alpha=1), Line2D([0], [0], color=c2[0], lw=8,alpha=1))
        vpt = fig.legend(custom_lines, ('Positive Label', 'Neg Label'),bbox_to_anchor=(.95,.25),fontsize = fs*1.5)

        vpt = plt.setp(axes[1].get_yticklabels(), visible=False)
        vpt = plt.setp(axes[2].get_yticklabels(), visible=False)
        vpt = plt.setp(axes[4].get_yticklabels(), visible=False)
        vpt = plt.setp(axes[5].get_yticklabels(), visible=False)
        
        vpt = plt.setp(axes[3].get_xticklabels(), visible=True)
        vpt = plt.setp(axes[4].get_xticklabels(), visible=True)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        fig.delaxes(axes[5])
        plt.show()
        plt.savefig(f'{prep.tmpBaseDir}/cell_{cell}.png')

        


if args.taskType == 'stats':
    cell = args.cellType

#     trainLabel = pd.read_csv(f"{prep.tmpBaseDir_hic}/{cell}/P-E/pairs.csv")['label']
#     p = sum(trainLabel); tot = len(trainLabel)
#     print(f"\n{cell}; pairs:: pos, neg: {p}/{tot}")
#     assert p==tot//2
    
#     trainLabel = pd.read_csv(f"{prep.tmpBaseDir_hic}/{cell}/P-E/pairs_test.csv")['label']
#     p = sum(trainLabel); tot = len(trainLabel)
#     print(f"{cell}; test:: pos, neg: {p}/{tot}")

    trainLabel = pd.read_csv(f"{prep.tmpBaseDir_hic}/{cell}/P-E/pairs_train_augment.csv")['label']
    p = sum(trainLabel); tot = len(trainLabel)
    print(f"{cell}; augment:: pos: {p}/{tot}")

#     trainLabel = pd.read_csv(f"{prep.tmpBaseDir_hic}/{cell}/P-E/pairs_train.csv")['label']
#     p = sum(trainLabel); tot = len(trainLabel)
#     print(f"_train:: pos, neg: {p}/{tot}")


    pr_dnase = np.load(f"{prep.tmpBaseDir_hic}/{cell}/P-E/promoter_DNase.npz",allow_pickle=True)['expr']
    pr_dna =  np.load(f"{prep.tmpBaseDir_hic}/{cell}/P-E/promoter_Seq.npz",allow_pickle=True)['sequence']
    print(cell,"promoter DNase",pr_dnase.shape,pr_dnase.dtype,flush=True)
    print('\t\t DNA',pr_dna.shape,pr_dna.dtype,flush=True)
    
    enh_dnase = np.load(f"{prep.tmpBaseDir_hic}/{cell}/P-E/enhancer_DNase.npz",allow_pickle=True)['expr']
    enh_dna =  np.load(f"{prep.tmpBaseDir_hic}/{cell}/P-E/enhancer_Seq.npz",allow_pickle=True)['sequence']
    print(cell,"enhancer DNase",enh_dnase.shape,enh_dnase.dtype,flush=True)
    print('\t\t DNA',enh_dna.shape,enh_dna.dtype,flush=True)

elif args.taskType == 'plotDNase':
    cell=args.cellType
    bamfs_cell = pt.flatten(prep.DnaseCells[cell])

    tasklist = bamfs_cell

    obj_func = lambda task:plot_DNase(cell,task)

    args.nTasks = min(args.nTasks,len(tasklist))
    print(f"gen DNase plots..", end=" ", flush=True)
    pt.distribute_task(task_list = tasklist, 
        nTasks = int(args.nTasks), file_index=int(args.file_index), 
        func=obj_func, num_tasks=len(tasklist),
        dry_run=False)
    print(".", flush=True)

elif args.taskType == 'tohdf5':
    cell_list = [args.cellType]

    tasklist = [f"{prep.tmpBaseDir}/promoterDNA_{cell}.hdf5" for cell in cell_list]
    tasklist += [f"{prep.tmpBaseDir}/enhancerDNA_{cell}.hdf5" for cell in cell_list]

    args.nTasks = min(args.nTasks,len(tasklist))
    print(f"convert .npz files to pd series hdf5 {args.file_index}/{len(tasklist)})..", end=" ", flush=True)
    pt.distribute_task(task_list = tasklist, 
        nTasks = int(args.nTasks), file_index=int(args.file_index), 
        func=onehotDNA_as_pdseries, num_tasks=len(tasklist),
        dry_run=False)
    print(".", flush=True)







