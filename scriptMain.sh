#!/bin/bash


optionsPrep=('splitBam' 'prepBam' 'prepPr' 'prepEnh' 'prepHiC')

optionsDNase=('genIndex' 'bamToArray' 'pgenProfile' 'egenProfile')
numArrayDNase=(30 700 30 30)

optionsDNA=('pDNA' 'eDNA')

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


if [[ ${optionsDNase[@]} == *$option* && $# -lt 2 ]]; then
	for i in "${!optionsDNase[@]}"; do
	   if [[ "${optionsDNase[$i]}" = "${option}" ]]; then
		   	sbatchopt="-a 0-${numArrayDNase[$i]} "
	   fi
	done
fi


if [[ ${optionsPrep[@]} == *$option* ]]; then
	sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option}
elif [[ ${optionsDNase[@]:1} == *$option* ]]; then
	sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option}
elif [[ ${optionsDNA[@]} == *$option* ]]; then
	sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option}
elif [[ ${optionsPCHiC[@]} == *$option* ]]; then
	for cell in ${cellList[@]}
	do
		sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option} ${cell}
	done
fi



## preprocessing
# scriptMain.sh splitBam
# scriptMain.sh prepBam
# scriptMain.sh prepPr
# scriptMain.sh prepEnh
# scriptMain.sh prepHiC

## DNA
# scriptMain.sh pDNA
# scriptMain.sh eDNA

## DNase
# scriptMain.sh genIndex
# scriptMain.sh bamToArray
# scriptMain.sh pgenProfile
# scriptMain.sh egenProfile


## Training Data
# scriptMain.sh hicMatch 13 : 21 min
# scriptMain.sh hicLabels


## 



