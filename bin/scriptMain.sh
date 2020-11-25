#!/bin/sh

#abort on error
set -e

optionsPrep=('splitBam' 'prepBam' 'prepPr' 'prepEnh' 'prepHiC')

optionsDNase=('genIndex' 'bamToArray' 'pgenProfile' 'egenProfile')
# numArrayDNase=(30 700 30 30)
numArrayDNase=(30 100 30 30)

optionsDNA=('pDNA' 'eDNA')

optionsPCHiC=('hicMatch' 'hicLabels' 'selectDNase' 'combineDNaseReps' 'selectDNA' 'sepData')
cellList=('tB' 'tCD4' 'nCD4' 'FoeT' 'Mon' 'tCD8')

optionsPy2=('bootstrap')

optionsNPZ=('stats' 'tohdf5' 'plotDNase')

function usage
{
    echo "usage: scriptMain.sh [-n AN_ARG || -c SOME_MORE_ARGS || -h] task_type"
    echo "   ";
    echo "  -n | --num-tasks    : number of sub-tasks";
    echo "  -c | --cell         : cell type";
    echo "  -r | --num-rep      : num of rep of DNase";
    echo "  -h | --help         : This message";
    exit 0
}

function parse_args
{
  # positional args
  args=()

  # named args
  while [ "$1" != "" ]; do
      case "$1" in
          -n | --num-tasks )       num_tasks="$2";      shift;;
          -c | --cell )                 cell="$2";      shift;;
          -r | --num-rep )              num_rep="$2";      shift;;
          -h | --help )                 usage;                   exit;; # quit and show usage
          * )                           args+=("$1")             # if no match, add it to the positional args
      esac
      shift # move to next kv pair
  done

  # restore positional args
  set -- "${args[@]}"

  # set positionals to vars
  option="${args[0]}"

  # only one positional argument required
  if [[ -z "$num_rep" ]]; then
      num_rep=1;
  fi
  if [[ $# -gt 1 ]]; then
    usage
  fi

}


function run
{

  parse_args "$@"

  if [[ ${optionsDNase[@]} == *$option* && -z "$num_tasks" ]]; then
    for i in "${!optionsDNase[@]}"; do
       if [[ "${optionsDNase[$i]}" = "${option}" ]]; then
          sbatchopt="-a 0-${numArrayDNase[$i]} "
       fi
    done
  elif [[ -n "$num_tasks" ]]; then
    sbatchopt="-a 0-${num_tasks} "
  fi

  if [[ -n "$cell" ]]; then
    cellList=(${cell})
  fi

  if [[ ${optionsPrep[@]} == *$option* ]]; then
    echo "sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option}"
    sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option}
  elif [[ ${optionsDNase[@]:1} == *$option* ]]; then
    echo "sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option}"
    sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option}
  elif [[ ${optionsDNA[@]} == *$option* ]]; then
    echo "sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option}"
    sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option}
  elif [[ ${optionsPCHiC[@]} == *$option* ]]; then
    for cell in ${cellList[@]}
    do
      echo "sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option} ${cell} ${num_rep}"
      sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option} ${cell} ${num_rep}
    done
  elif  [[ ${optionsPy2[@]} == *$option* ]]; then
    for idx in $(seq 0 $num_tasks)
    do
      echo "sbatch -J ${option} scriptPy2 ${idx}"
      sbatch -J ${option} scriptPy2 ${idx}
    done
  elif  [[ ${optionsNPZ[@]} == *$option* ]]; then
    for cell in ${cellList[@]}
    do
      echo "sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option} ${cell}"
      sbatch -J ${option} ${sbatchopt}scriptDT.sh ${option} ${cell}
    done
  fi
}



run "$@";


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
# scriptMain.sh bamToArray : 10 min
# scriptMain.sh pgenProfile : 10 min
# scriptMain.sh egenProfile : 10 min


## Training Data
# scriptMain.sh.sh hicMatch -n 13 : 21 min
# scriptMain.sh.sh hicLabels : 9 hours


## Select Input Features for Training Data : DNA
# scriptMain.sh selectDNA; task = split : 3 min
# scriptMain.sh selectDNA -n 100; task = run : 25 min
# scriptMain.sh selectDNA; task = combine ; 8 min

## Select Input Features for Training Data : DNase
# scriptMain.sh -n 25 selectDNase : 21+ hours
# scriptMain.sh -n 2 combineDNaseReps : 5 min
# scriptMain.sh -c <CELL> sepData : 1 hours


## Clean-Up
# find . -type f -regex ".*_[0-9]*\.npz" | rm

# DeepTact train Demo : Mon : ENCFF295OEK
# python BootStrap : 3.25 hours

## to hdf5
# scriptMain.sh tohdf5 -n 2 [-c FoeT]
# scriptMain.sh plotDNase -n 11 [-c FoeT]



