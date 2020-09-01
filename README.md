# Input Data Pre-Processing 

- Input parameters & files : `parameters.py`  
- Pre-processing output names : `preprocessing.py`  
- Pre-process **DNA-Seq** data : `process_fasta.py`  
- Pre-process **DNase-Seq** data : `process_Dnase.py`
- Pre-process **PCHiC** data : `process_PCHiC.py`

# Pre-Processing Script 


1. Create pre-processing directories. Generate main and background .bed files for promoters and enhancers : `scriptMain.sh prep`

2. Generate promoter and enhancer one hot encoded DNA-Sequence : `scriptMain.sh DNA`

3. Generate regulatory elements' DNase-Seq data : 
	1. split .bam file for bedtools intersect : `scriptMain.sh splitBam 300`
	2. create promoter and enhancer window .bed files:  `scriptMain.sh Win 300`
	3. Generate *bedtools intersect* tasklist: `scriptMain.sh TaskList`
	4. Perform *bedtools intersect*: `scriptMain.sh Intersect 300`
	5. Combine *bedtools intersect* output. Create tasklist: `scriptMain.sh cTaskList`
	6. Perform *bedtools intersect* output combine: `scriptMain.sh Combine 300`
	7. Post Processing :  `scriptMain.sh PostProcess 30`

3. Generate Training Data:
	1. HiCMatch: `scriptMain.sh hicMatch 13`
	2. Generate Labels: `scriptMain.sh hicLabels`


# BAM Pre-Processing

## Bamtools Split

split .bam files before preprocessing :

	bamtools split -in ENCFF138MTV.bam -reference


## Bedtools Intersect vs Bedtools Coverage

Common Options:  
`-f` : counts intersections as fraction of input A  
`-F` : counts intersections as fraction of input B  
`-e` : require that the minimum fraction be satisfied for A _OR_ B

**Coverage**
\[[Doc](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html)\]  
bedtool coverage : coverage of features in file B on the features in file A
	usage : bedtool coverage [options] `-a` [input file A] `-b` [input file B]
`-counts` :  only report counts of intersection. (Restricted by `-f `and `-r`)

	bedtools coverage -counts -e -f 0.5 -F 0.5  -a <bed-file> -b <bam-file>  > <out-file>"

**Intersect**  
`-c`: For each entry in A, report the number of hits in B 

	bedtools intersect -c -e -f 0.5 -F 0.5 -a <bed-file> -b <bam-file>  > <out-file>


### Conclusion
For bed file with \~200,000 lines and 852MB bam file, the memory and time requirements are :
bedtools intersect : 80s CPU-time, 23.15GB memory utilized
bedtools coverage : 82s CPU-time, 22.55GB memory utilized
 
So they are similar in both time and memory. But, if there is enough memory **bedtools intersect is better**. And the difference in memory requirements is not too large anyways.
 







