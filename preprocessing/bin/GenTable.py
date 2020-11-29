import sys, os
import glob
import preptools as pt

baseHiC = "/projects/li-lab/agarwa/CUBE/DeepTact/dataset/HiC"
tmp = "ENCSR046XXF/released/reads/fastq/rep1/pair1"
# fastqF = [f"{baseHiC}/{tmp}/ENCFF047VMB.part-01.fastq",f"{baseHiC}/{tmp}/ENCFF047VMB.part-02.fastq",]
fastqF = glob.glob(f"{baseHiC}/**/*.fastq.gz",recursive=True)




def find_lines(fastqF,outp):
    cmdL = lambda fn,outp:os.system(f"echo $(zcat {fn} | wc -l)/4|bc >> {outp}")
    with open(outp, 'a') as f:
        f.write(f"{os.path.basename(fastqF)}: ") 
    result = cmdL(fastqF,outp) 
    print(result,flush=True)
    return result



if __name__=="__main__":

    args = pt.process_inputArgs(input_parse=sys.argv[1:])

    outputf = f"{baseHiC}/reads.{args.file_index}.txt"

    pt.distribute_task(task_list = fastqF, nTasks = int(args.nTasks), file_index=int(args.file_index), func=lambda x:find_lines(x,outputf), dry_run=False)