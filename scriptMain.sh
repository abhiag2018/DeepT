#!/bin/bash


optionsPrep1=('prep' 'prepPr' 'prepEnh' 'prepHiC')
optionsPrep2=('splitBam')
optionsPrep=("${optionsPrep1[@]}" "${optionsPrep2[@]}")

optionsProc1=('Win' 'TaskList' 'Intersect' 'cTaskList' 'Combine' 'PostProcess')
optionsProc2=('' 'p' 'e')
optionsProc=()
for _t in ${optionsProc1[@]}
do
	for t in ${optionsProc2[@]}
	do
		optionsProc+=(${t}${_t})
	done
done

optionsDNA=('DNA' 'pDNA' 'eDNA')

optionsPCHiC=('hicMatch','hicLabels')
cellList=('tB' 'tCD4' 'nCD4' 'FoeT' 'Mon' 'tCD8')


options=("${optionsPrep[@]}" "${optionsProc[@]}" "${optionsDNA[@]}" "${optionsPCHiC[@]}")

helpFunction()
{
	echo -e """usage: scriptMain.sh <TaskType> \n\t\tin ${options[@]} \n\n"""
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
	if [ $# -gt 1 ]; then
		sbatchopt="-a 0-${2} "
	else
		sbatchopt=""
	fi
fi



if [[ ${optionsPrep1[0]} == *$option* ]]; then
	for opt in ${optionsPrep1[@]:1}
	do
		sbatch -J ${opt} ${sbatchopt}scriptDT.sh ${opt}
	done
elif [[ ${optionsPrep2[0]} == *$option* ]]; then
	sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option}
elif [[ ${optionsProc1[@]} == *$option* ]]; then
	for _t in ${optionsProc2[@]:1}
	do
		opt=(${_t}${option})
		sbatch -J ${opt} ${sbatchopt}scriptDT.sh ${opt}
	done
elif [[ ${optionsDNA[0]} == *$option* ]]; then
	for opt in ${optionsDNA[@]:1}
	do
		sbatch -J ${opt} ${sbatchopt}scriptDT.sh ${opt}
	done
elif [[ ${optionsPCHiC[@]} == *$option* ]]; then
	for cell in ${cellList[@]}
	do
		sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option} ${cell}
	done
fi





