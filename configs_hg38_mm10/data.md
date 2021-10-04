# Data Download

dataDir="/projects/li-lab/agarwa/CUBE/DeepTact/dataset"

## Promoter List

hg38 : https://genome.ucsc.edu/cgi-bin/hgTables

	clade : Mammal
	genome : human
	assembly : Dec 2013 (hg38)
	group : genes and gene predictions
	track : NCBI RefSeq
	table : RefSeq All (ncbiRefSeq)
	output format : all fields from selected table

local file: $dataDir/hg38_ncbiRefSeq.bed

mm10 :  https://genome.ucsc.edu/cgi-bin/hgTables


	clade : Mammal
	genome : mouse
	assembly : Dec 2011 (mm10)
	group : genes and gene predictions
	track : NCBI RefSeq
	table : RefSeq All (ncbiRefSeq)
	output format : all fields from selected table

local file: $dataDir/mm10_ncbiRefSeq.bed


### GTF file

https://genome.ucsc.edu/cgi-bin/hgTables

output format : GTF

## Enhancer List

hg38 : https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer/F5.hg38.enhancers.bed.gz

local file: $dataDir/enhancers_fantom5/F5.hg38.enhancers.bed

mm10 : https://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/enhancer/F5.mm10.enhancers.bed.gz

local file: $dataDir/enhancers_fantom5/F5.mm10.enhancers.bed

## Genome Sequence

hg38 : https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/latest/hg38.fa.gz

mm10 : https://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/

