#!/bin/bash


optionsPrep=('prep' 'prepBam' 'prepPr' 'prepEnh' 'prepHiC')

optionsDNase=('DNase' 'genIndex' 'bamToArray' 'pgenProfile' 'egenProfile')
numArrayDNase=(30 700 30 30)

optionsDNA=('DNA' 'pDNA' 'eDNA')

optionsPCHiC=('hicMatch','hicLabels')
cellList=('tB' 'tCD4' 'nCD4' 'FoeT' 'Mon' 'tCD8')


options=("${optionsPrep[@]}" "${optionsDNase[@]}" "${optionsDNA[@]}" "${optionsPCHiC[@]}")

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


if [[ ${optionsDNase[@]:1} == *$option* && $# -lt 2 ]]; then
	tmp=("${optionsDNase[@]:1}")
	for i in "${!tmp[@]}"; do
	   if [[ "${tmp[$i]}" = "${option}" ]]; then
		   	sbatchopt="-a 0-${numArrayDNase[$i]} "
	   fi
	done
fi


if [[ ${optionsPrep[0]} == *$option* ]]; then
	for opt in ${optionsPrep[@]:1}
	do
		sbatch -J ${opt} ${sbatchopt}scriptDT.sh ${opt}
	done
elif [[ ${optionsPrep[@]:1} == *$option* ]]; then
	sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option}
elif [[ ${optionsDNase[0]} == *$option* ]]; then
	tmp=("${optionsDNase[@]:1}")
	for i in "${!tmp[@]}"; do
	   	sbatchopt="-a 0-${numArrayDNase[$i]} "
		sbatch -J ${tmp[$i]} ${sbatchopt}scriptDT.sh ${tmp[$i]}
	done
elif [[ ${optionsDNase[@]:1} == *$option* ]]; then
	# echo "sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option}"
	sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option}
elif [[ ${optionsDNA[0]} == *$option* ]]; then
	for opt in ${optionsDNA[@]:1}
	do
		sbatch -J ${opt} ${sbatchopt}scriptDT.sh ${opt}
	done
elif [[ ${optionsDNA[@]:1} == *$option* ]]; then
	sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option}
elif [[ ${optionsPCHiC[@]} == *$option* ]]; then
	for cell in ${cellList[@]}
	do
		sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option} ${cell}
	done
else
	sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option}
fi




# scriptMain.sh prep
# scriptMain.sh DNase
# scriptMain.sh DNA


# scriptMain.sh hicLabels 30



