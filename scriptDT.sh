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

optionsPrep=('splitBam' 'prep' 'prepPr' 'prepEnh')

optionsProc=()
for _t in 'Win' 'TaskList' 'Intersect' 'cTaskList' 'Combine' 'PostProcess'
do
	for t in '' 'p' 'e'
	do
		optionsProc+=(${t}${_t})
	done
done

optionsDNA=('DNA' 'pDNA' 'eDNA')

options=("${optionsPrep[@]}" "${optionsProc[@]}" "${optionsDNA[@]}")
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

# python process_PCHiC.py
# python download_DNase.py

if [[ ${optionsPrep[@]} == *$option* ]]; then
    echo  "run preprocessing.py"
    python preprocessing.py --nTasks ${num_tasks} --taskType ${option}  --file_index=$SLURM_ARRAY_TASK_ID
elif [[ ${optionsDNA[@]} == *$option* ]]; then
	echo  "run process_fasta.py"
	python process_fasta.py --nTasks ${num_tasks} --taskType ${option}  --file_index=$SLURM_ARRAY_TASK_ID	
else
	echo  "run process_Dnase.py"
	python process_Dnase.py --nTasks ${num_tasks} --taskType ${option}  --file_index=$SLURM_ARRAY_TASK_ID
fi


# sbatch -J prepPr scriptDT.sh prepPr
# sbatch -J prepEnh scriptDT.sh prepEnh
# # sbatch -J splitBam -a 0-300 scriptDT.sh splitBam
# sbatch -J pWin -a 0-300 scriptDT.sh pWin
# sbatch -J eWin -a 0-300 scriptDT.sh eWin

# sbatch -J pTaskList scriptDT.sh pTaskList
# sbatch -J eTaskList scriptDT.sh eTaskList
# sbatch -J pIntersect -a 0-300 scriptDT.sh pIntersect
# sbatch -J eIntersect -a 0-300 scriptDT.sh eIntersect

# sbatch -J pcTaskList scriptDT.sh pcTaskList
# sbatch -J ecTaskList scriptDT.sh ecTaskList

# sbatch -J pCombine -a 0-300 scriptDT.sh pCombine
# sbatch -J eCombine -a 0-300 scriptDT.sh eCombine
# sbatch -J pPostProcess -a 0-30 scriptDT.sh pPostProcess
# sbatch -J ePostProcess -a 0-30 scriptDT.sh ePostProcess


###

# sbatch -J pDNA scriptDT.sh pDNA
# sbatch -J eDNA scriptDT.sh eDNA








