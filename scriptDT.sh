#!/bin/bash
#SBATCH --array=0-0
#SBATCH --output=./slurm-log/slurm-%A_%a.out
#SBATCH --error=./slurm-log/slurm-%A.out
#SBATCH --qos batch 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500G
#SBATCH --time 12:00:00

num_tasks=$SLURM_ARRAY_TASK_COUNT

optionsPrep=('prepBam' 'prepPr' 'prepEnh' 'prepHiC')

optionsDNase=('genIndex' 'bamToArray' 'pgenProfile' 'egenProfile')

optionsDNA=('pDNA' 'eDNA' 'selectDNA')

optionsPCHiC=('hicMatch' 'hicLabels')

options=("${optionsPrep[@]}" "${optionsDNase[@]}" "${optionsDNA[@]}" "${optionsPCHiC[@]}")
# echo ${options[@]}
# echo ${#options[@]}


helpFunction()
{
	echo -e """usage: scriptDT.sh <TaskType> \n\t\tin ${options[@]} \n\n"""
	exit 0
}


option=""
if [ $# -lt 1 ]
  then
	helpFunction
else
	option=$1
	if [[ ! ${options[@]} == *$option* ]]; then
	    helpFunction
	fi
fi

# echo "SLURM_CPUS_ON_NODE=$SLURM_CPUS_ON_NODE"
# echo "SLURM_MEM_PER_CPU=$SLURM_MEM_PER_CPU"
# echo "SLURM_MEM_PER_NODE=$SLURM_MEM_PER_NODE"

# echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"
# echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"

#####
source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

if [[ ${optionsPrep[@]} == *$option* ]]; then
    echo  "run preprocessing.py"
    python preprocessing.py --nTasks ${num_tasks} --taskType ${option}  --file_index=$SLURM_ARRAY_TASK_ID
elif [[ ${optionsDNA[2]} == *$option* ]]; then
	echo  "run process_fasta.py"
	cellType=$2
	python process_fasta.py --nTasks ${num_tasks} --taskType ${option}  --file_index=$SLURM_ARRAY_TASK_ID --cellType=$cellType
elif [[ ${optionsDNA[@]} == *$option* ]]; then
	echo  "run process_fasta.py"
	python process_fasta.py --nTasks ${num_tasks} --taskType ${option}  --file_index=$SLURM_ARRAY_TASK_ID	
elif [[ ${optionsPCHiC[@]} == *$option* ]]; then
	echo  "run process_PCHiC.py"
	cellType=$2
	python process_PCHiC.py --nTasks ${num_tasks} --taskType ${option}  --file_index=$SLURM_ARRAY_TASK_ID --cellType=$cellType
else
	echo  "run process_Dnase.py"
	python process_Dnase.py --nTasks ${num_tasks} --taskType ${option}  --file_index=$SLURM_ARRAY_TASK_ID
fi


# sbatch -J prepBam scriptDT.sh prepBam
# sbatch -J prepPr scriptDT.sh prepPr
# sbatch -J prepEnh scriptDT.sh prepEnh
# sbatch -J prepHiC scriptDT.sh prepHiC
# # sbatch -J splitBam -a 0-30 scriptDT.sh splitBam

# sbatch -J genIndex -a 0-30 scriptDT.sh genIndex
# sbatch -J bamToArray -a 0-700 scriptDT.sh bamToArray
# sbatch -J pgenProfile -a 0-30 scriptDT.sh pgenProfile
# sbatch -J egenProfile -a 0-30 scriptDT.sh egenProfile

# sbatch -J pDNA scriptDT.sh pDNA
# sbatch -J eDNA scriptDT.sh eDNA
# sbatch -J selectDNA scriptDT.sh selectDNA tB
# sbatch -J selectDNA scriptDT.sh selectDNA tCD4
# sbatch -J selectDNA scriptDT.sh selectDNA nCD4
# sbatch -J selectDNA scriptDT.sh selectDNA FoeT
# sbatch -J selectDNA scriptDT.sh selectDNA Mon
# sbatch -J selectDNA scriptDT.sh selectDNA tCD8


###

# sbatch -J hicMatch -a 0-13  scriptDT.sh hicMatch tB
# sbatch -J hicMatch -a 0-13  scriptDT.sh hicMatch tCD4
# sbatch -J hicMatch -a 0-13  scriptDT.sh hicMatch nCD4
# sbatch -J hicMatch -a 0-13  scriptDT.sh hicMatch FoeT
# sbatch -J hicMatch -a 0-13  scriptDT.sh hicMatch Mon
# sbatch -J hicMatch -a 0-13  scriptDT.sh hicMatch tCD8

# sbatch -J hicLabels scriptDT.sh hicLabels tB
# sbatch -J hicLabels scriptDT.sh hicLabels tCD4
# sbatch -J hicLabels scriptDT.sh hicLabels nCD4
# sbatch -J hicLabels scriptDT.sh hicLabels FoeT
# sbatch -J hicLabels scriptDT.sh hicLabels Mon
# sbatch -J hicLabels scriptDT.sh hicLabels tCD8









