# Tools/Resources

UCSC : http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/  
BEDOPS : https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/wig2bed.html  
GENOME base conversion : https://genome.ucsc.edu/cgi-bin/hgLiftOver

# Conda Packages

1. Install Required Packages for running Neural Network
env name : dt1

		conda create -p dt1 python=2.7  
		conda install biopython=1.70 --yes  
		conda install tensorflow=0.12 --yes  
		conda install pandas=0.20.1 --yes  
		pip install keras==1.2.0  
		conda install hickle  
		conda install scikit-learn=0.18.2  
		pip install theano==0.9.0  

2. Change backend to theano in "\~/.keras/keras.json"


# Input Download

- Promoter Capture Hi-C from 17 cell types from Javierre et al [ref 18]

- DNase-Seq data from ENCODE : 
		python encode_downloader.py --file-types bam:alignments --ignore-unpublished --assemblies hg19 --dir ../dataset/Dnase-Seq --max-download 26 ../dataset/Dnase-Seq/experiments_dnase.txt

- Permissive enhancers from [FANTOM5](https://fantom.gsc.riken.jp/5/) [ref 21]: Download all files in [weblink](https://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/) : 
		wget -i /projects/li-lab/agarwa/CUBE/DeepTact/dataset/enhancers_fantom5/dl_list.txt

- Promoter bed files from [UCSC gene browser](https://genome.ucsc.edu/cgi-bin/hgTables) with options :
		1. clade : Mammal
		2. genome : Human
		3. assembly : GRCh37/hg19 Feb 2009
		4. group : Genes and Gene Predictions
		5. track : Ensembl Genes
		6. table : ensGene
		7. region : genome
		8. output format : all fields from selected table
		9. output file : <filename>

		Format info : click on describe table schema button

- Human hg19_fasta files from [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz)

- HiC ENCODE Data : 

		python encode_downloader.py --file-types fastq --ignore-unpublished --dir ../dataset/HiC --max-download 26 ../dataset/table_HiC.txt


# Input Data Format

Located at demo/P-E in .npz file format. For each of promoter and enhancer as \<elem\>
- `elem_DNase.npz` : contains variable `"expr"`
- `elem_Seq.npz` : contains variables `"label"`, `"sequence"`


## DNA Sequence Input
	pr = np.load('promoter_Seq.npz')
	enh = np.load('enhancer_Seq.npz')

- `pr['label']` : array of size *R* of indicators denoting whether the current pair is DRE-interaction or not
- `enh['label'] == pr['label']`
- `pr['sequence']` : array of size *R x 4000*
- `pr['sequence'].reshape(*R*,1000,4)` gives the appropriate one hot encoded tensor containing the DNA Sequence of each element. Similarly, `enh['sequence'].reshape(*R*,2000,4)` gives the appropriate one hot encoded tensor.

## DNase-Seq Input

	pr = np.load('promoter_DNase.npz')
	enh = np.load('promoter_DNase.npz')

- `pr['expr']` : array of size *R x (1000 \* rep)* where rep is the number of biological replicates for the data.
- `pr['expr'].reshape(R,1000,rep)` :  gives the appropriate DNase-Seq data tensor. 

Similarly for `enh['expr']`.

