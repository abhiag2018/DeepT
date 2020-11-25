import re, os
import subprocess
import glob
import itertools
import httplib2
import urllib3
from urllib.parse import unquote
from urllib.request import urlretrieve
from bs4 import BeautifulSoup, SoupStrainer


base_dataDir = "/projects/li-lab/agarwa/CUBE/DeepTact/dataset"

def download_experiments(exprLi):
    http = httplib2.Http()
    exprLi = {'tB':exprLi[:5+1],'tCD4':exprLi[6:9],'nCD4':exprLi[9:12],'FoeT':exprLi[12:23],'Mon':exprLi[23:28],'tCD8':exprLi[28:]}
    for cell,expr_cell in itertools.islice(exprLi.items(),0,None):
        for expr in expr_cell:
            status, response = http.request(f'https://www.encodeproject.org/experiments/{expr}/')
            dirname = f"{base_dataDir}/Dnase-Seq/{cell}/{expr}"
            for link in BeautifulSoup(response, parse_only=SoupStrainer('a'),features="html.parser"):
                if link.has_attr('href') and re.match('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi\?acc=GSM',link['href']):
                    try:
                        status1, response1 = http.request(link['href'])
                    except BrokenPipeError:
                        print("broken link : ",link['href'])
                        continue
                    for link1 in BeautifulSoup(response1, parse_only=SoupStrainer('a'),features="html.parser"):
                        if link1.has_attr('href') and (re.match('ftp:.*bed%2Egz',link1['href']) or re.match('ftp:.*wig%2Egz',link1['href']) or re.match('ftp:.*narrowPeak%2Egz',link1['href']) or re.match('ftp:.*broadPeak%2Egz',link1['href']) or re.match('ftp:.*bigWig',link1['href'])):
                            fpath = f"{dirname}/{os.path.basename(unquote(link1['href']))}"
                            if not os.path.exists(fpath):
                                os.makedirs(dirname,exist_ok=True)
                                print("downloading ",link1['href'],".. ",end="")
                                # urlretrieve(link1['href'], fpath)
                                print(" .")
            print()

def encode_exp_biosample_summary(exprLi):
    http = httplib2.Http()
    exprLi = {'tB':exprLi[:5+1],'tCD4':exprLi[6:9],'nCD4':exprLi[9:12],'FoeT':exprLi[12:23],'Mon':exprLi[23:28],'tCD8':exprLi[28:]}
    for cell,expr_cell in itertools.islice(exprLi.items(),0,None):
        for expr in expr_cell:
            status, response = http.request(f'https://www.encodeproject.org/experiments/{expr}/')
            dirname = f"{base_dataDir}/Dnase-Seq/{cell}/{expr}"
            for link in BeautifulSoup(response, parse_only=SoupStrainer('a'),features="html.parser"):
                if link.has_attr('href') and re.match('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi\?acc=GSM',link['href']):
                    try:
                        status1, response1 = http.request(link['href'])
                    except BrokenPipeError:
                        print("broken link : ",link['href'])
                        continue
                    for link1 in BeautifulSoup(response1, parse_only=SoupStrainer('a'),features="html.parser"):
                        if link1.has_attr('href') and (re.match('ftp:.*bed%2Egz',link1['href']) or re.match('ftp:.*wig%2Egz',link1['href']) or re.match('ftp:.*narrowPeak%2Egz',link1['href']) or re.match('ftp:.*broadPeak%2Egz',link1['href']) or re.match('ftp:.*bigWig',link1['href'])):
                            fpath = f"{dirname}/{os.path.basename(unquote(link1['href']))}"
                            if not os.path.exists(fpath):
                                os.makedirs(dirname,exist_ok=True)
                                print("downloading ",link1['href'],".. ",end="")
                                # urlretrieve(link1['href'], fpath)
                                print(" .")
            print()

if __name__ == "__main__":
    # exprLi = ['ENCSR000EMJ','ENCSR000EMJ','ENCSR247IUJ','ENCSR646BJT','ENCSR891VOV','ENCSR643GHI',
    #     'ENCSR425MRQ','ENCSR569ATD','ENCSR000EML','ENCSR000EML','ENCSR000EMM','ENCSR000EMM','ENCSR164PQJ',
    #     'ENCSR342NCX','ENCSR541UPY','ENCSR541UPY','ENCSR544KDB','ENCSR682PKY','ENCSR777USE','ENCSR838IPF',
    #     'ENCSR888UPQ','ENCSR917VCP','ENCSR962EAP','ENCSR000ELE','ENCSR000ELE','ENCSR674JIL','ENCSR000EPK',
    #     'ENCSR695AUY','ENCSR840TVG','ENCSR609DDQ','ENCSR842VTJ']
    # download_experiments(exprLi)

    with open('DeepTact/dataset/ATAC-Seq/experiments_atac.txt', "r") as text_file:
        exp_atac = text_file.readlines()
        exp_atac = [exp[:-1] for exp in exp_atac]
        exp_atac

    # for file in glob.glob(f"{base_dataDir}/Dnase-Seq/**/**/*.bigWig"):
    #     print("converting", os.path.basename(file),"..",end=" ")
    #     cmd = f"bigWigToBedGraph {file} {file[:-6]}bedGraph"
    #     result = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
    #     print(".")
    #     print(result,end="\n\n")











