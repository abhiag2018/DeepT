import re
# base_dir = "/projects/li-lab/agarwa/CUBE/DeepTact/dataset"
# for chrid in list(range(1,23))+list('XY'):
#     fasta_path = base_dir+f"/genome_seq/Homo_sapiens.GRCh37.75.dna.chromosome.{chrid}.fa"
#     new_fasta_path = base_dir+f"/genome_seq/Homo_sapiens.GRCh37.75.dna.chromosome.chr{chrid}.fa"
#     nf = open(new_fasta_path,'w')
#     with open(fasta_path,'r') as ff:
#         for line in ff:
#             if re.match(">",line):
#                 tmp=nf.write(re.sub(r">(\d+)", r">chr\1", line))
#             else:
#                 tmp=nf.write(line);
#     nf.close()



if __name__=="__main__":
    import preprocessing as prep
    baseDir = prep.baseDataDir
    hg19 = prep.hg19
    p = prep.promoter
    e = prep.enhancer
    bgwin = prep.bgWindow

    bamfs = prep.bamfiles
    interesctOpt = prep.intersectOptions

    reRun = prep.reRun

    #promoter
    pt.generate_elem_fa(f"{baseDir}/{hg19}", f"{baseDir}/{p['bed-path']}", out_fa = f"{baseDir}/{p['fa-out']}")
    prName,prLoc,prSeq,prOhc = pt.fasta_to_onehot(f"{baseDir}/{p['fa-out']}")
    

    #enhancer
    pt.generate_elem_fa(f"{baseDir}/{hg19}", f"{baseDir}/{e['bed-path']}", out_fa = f"{baseDir}/{e['fa-out']}")
