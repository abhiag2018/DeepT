#!/bin/bash
#SBATCH --array=0-0
#SBATCH --output=./slurm-log/slurm-%A_%a.out
#SBATCH --qos batch 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500G
#SBATCH --time 70:00:00

num_tasks=$SLURM_ARRAY_TASK_COUNT

optionsPrep=('splitBam' 'prep' 'prepPr' 'prepEnh')

optionsProc=()
for _t in 'Win' 'Intersect' 'Combine' 'PostProcess' 'TaskList'
do
	for t in '' 'p' 'e'
	do
		optionsProc+=(${t}${_t})
	done
done

# options=("${optionsPrep[@]}")
# options+=("${optionsProc[@]}")
options=("${optionsPrep[@]}" "${optionsProc[@]}")
# echo ${options[@]}
# echo ${#options[@]}


helpFunction()
{
	echo -e """usage: run_script.sh <TaskType> \n\t\tin ${options[@]} \n\n"""
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
else
	echo  "run process_Dnase.py"
	python process_Dnase.py --nTasks ${num_tasks} --taskType ${option}  --file_index=$SLURM_ARRAY_TASK_ID
fi


# sbatch run_script.sh prepPr
# sbatch run_script.sh prepEnh
# # sbatch -a 0-300 run_script.sh splitBam
# sbatch -a 0-300 run_script.sh pWin
# sbatch -a 0-300 run_script.sh eWin
# sbatch run_script.sh pTaskList
# sbatch run_script.sh eTaskList
# sbatch -J pIntersect -a 0-300 run_script.sh pIntersect
# sbatch -J eIntersect -a 0-300 run_script.sh eIntersect
# sbatch -a 0-300 run_script.sh pCombine
# sbatch -a 0-300 run_script.sh eCombine
# sbatch -a 0-300 run_script.sh pPostProcess
# sbatch -a 0-300 run_script.sh ePostProcess