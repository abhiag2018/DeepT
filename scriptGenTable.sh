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

options=('HiC')

# options=("${optionsPrep[@]}" "${optionsProc[@]}" "${optionsDNA[@]}" "${optionsHiC[@]}")
# echo ${options[@]}
# echo ${#options[@]}


helpFunction()
{
	echo -e """usage: scriptGenTable.sh <TaskType> \n\t\tin ${options[@]} \n\n"""
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


if [[ ${options[@]} == *$option* ]]; then
	echo  "run HiC.py"
	python HiC.py --nTasks ${num_tasks} --taskType ${option}  --file_index=$SLURM_ARRAY_TASK_ID	
fi

# sbatch -J HiC -a 0-300 run_script.sh HiC













