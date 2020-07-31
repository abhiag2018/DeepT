#!/bin/bash
#SBATCH --array=0-0
#SBATCH --output=./slurm-log/slurm-%A_%a.out
#SBATCH --qos batch 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500G
#SBATCH --time 70:00:00

num_tasks=$SLURM_ARRAY_TASK_COUNT

helpFunction()
{
	printf """usage: run_script.sh <TaskType> \n\t\tin prep|prepPr|prepEnh|pWin|pIntersect|pPostProcess|eWin|eIntersect|ePostProcess \n\n"""
	exit 1
}


option=""
if [ $# -lt 1 ]
  then
	helpFunction
else
	option=$1
	if [[ ! "$option" =~ ^(splitBam|prep|prepPr|prepEnh|pWin|pIntersect|pPostProcess|eWin|eIntersect|ePostProcess) ]]; then
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

if [[ "$option" =~ (splitBam|prep|prepPr|prepEnh) ]]; then
    echo  "run preprocessing.py"
    python preprocessing.py --nTasks ${num_tasks} --taskType ${option}  --file_index=$SLURM_ARRAY_TASK_ID
else
	echo  "run process_Dnase.py"
	python process_Dnase.py --nTasks ${num_tasks} --taskType ${option}  --file_index=$SLURM_ARRAY_TASK_ID
	# python process_Dnase.py --taskType ${option}  --nTasks=${num_tasks}  --file_index=$SLURM_ARRAY_TASK_ID
fi


# sbatch run_script.sh prep
# sbatch -a 0-300 run_script.sh splitBam
# sbatch -a 0-300 run_script.sh pWin
# sbatch -a 0-300 run_script.sh eWin
# sbatch -a 0-300 run_script.sh pIntersect
# sbatch -a 0-300 run_script.sh eIntersect
# sbatch -a 0-300 run_script.sh pPostProcess
# sbatch -a 0-300 run_script.sh ePostProcess