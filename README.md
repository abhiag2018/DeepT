# Input Data Pre-Processing 

- `parameters.py`   : Input parameters & files
- `preprocessing.py`   : Intermediate/Output path
- `process_fasta.py`   : Pre-process **DNA-Seq** data for all regulatory elements
- `process_Dnase.py` : Pre-process **DNase-Seq** data for all regulatory elements
- `process_PCHiC.py` : 
	- Pre-process **PCHiC** data to generate training labels.
	- Select input features for training data.
- `DeepTact_0.py` : DeepTact code ported to python 3.

# Pre-Processing with SLURM


1. Create pre-processing directories. Generate main and background .bed files for promoters and enhancers : `scriptMain.sh prep`

2. Generateone hot encoded DNA-Sequence for all elements in  promoter and enhancer lists. Each element is of length *extended-window* : `scriptMain.sh DNA`

3. Generate regulatory elements' DNase-Seq data : 

	1. `scriptMain genIndex`
	2. `scriptMain bamToArray`
	3. `scriptMain pgenProfile`
	4. `scriptMain egenProfile`

4. Generate Training .csv files with augmented training data :
	1. HiCMatch: `scriptMain.sh hicMatch -n 13`
	2. Run `scriptMain.sh hicLabels` for steps 1 to 7 in script. Complete step 6 running python script seperately.

5. Generate the required DNase and DNA input features by selecting the appropriate elements (and window) from the features generated in steps 2 and 3. 
	i. `scriptMain selectDNA` for tasks `split`, `run` and `combine`.
	ii. `scriptMain -n 25 selectDNase` and `scriptMain -n 2 combineDNaseReps`. `scriptMain.sh -c <CELL>  -r <NUM_REP> sepData`  to seperate the data points into seperate files for processing with tensorflow dataset generator.

6.
	- Separate the data points into training and validation. Then generate bootstrapped dataset indices from training data : `DeepTact_0.split_train_val_bootstrap(val=0.2)`. 
	-  Save the index list for validation and training sets at `CELL+'/'+TYPE+'/bagData/val.hkl` and  `CELL+'/'+TYPE+'/bagData/train_{i}.hkl`, respectively.


7. Train DeepTact model : `scriptGpU.sh -a 0-<NUM_ENSEMBL-1> <CELL> <NUM_REP>`

