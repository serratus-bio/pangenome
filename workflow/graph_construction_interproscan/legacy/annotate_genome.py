import xmltodict
import pickle as pkl
import tqdm
from pathlib import Path
import pandas as pd
import os
import warnings
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq

from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
from typing import Tuple,Dict,Union,List
# from typing import Union
Entrez.email='zfevan.deng@mail.utoronto.ca'

def fetch_seq(f:str)->str:
    gb:dict=pkl.load(open(f,'rb'))
    seq:str=gb['GBSet']['GBSeq']['GBSeq_sequence']
    return seq

translate= lambda x:Seq(x).translate()._data.decode()

def hextranslate(g:str)->List[str]:
    '''
    genome: input nt seq
    output: list of sense trans*3 + antisense trans*3
    '''
    o=[]
    genome=Seq(g)
    genome_r:Seq=genome.reverse_complement()
    for i in [0,1,2]:
        o.append(genome[i:].translate()._data.decode())
        o.append(genome_r[i:].translate()._data.decode())
    return o

get_hex = lambda f: hextranslate(fetch_seq(f))
routine_dict=Dict[str,Union[str,int,float]]
transannot=['0s','0a','1s','1a','2s','2a']
transannot_indice={j:i for i,j in enumerate(transannot)}

def get_valid_seg(translist:List[str],name:str='xx',thresh:int=100)->List[Tuple[routine_dict,str]]:
    output=[]
    for i,annot in zip(translist,transannot):
        b=e=0
        for seg in i.split('*'):
            e=b+len(seg)
            if len(seg)>thresh:
                head={'name':name,
                      'transannot':annot,
                      'prob':b,
                      'proe':e}
                output.append((head,seg))
            b=e+1
    return output

def seg_to_fasta(seg:List[Tuple[routine_dict,str]])->str:
    o=[]
    for s in seg:
       head='> ' + ','.join([f'{k}={v}' for k,v in s[0].items()])
       o.extend([head,s[1]])
    return '\n'.join(o)

def fetch_fasta(idir:str,odir:str,thresh:int=100):
    odict={}
    idir:Path=Path(idir)
    odir:Path=Path(odir)
    odir.mkdir(exist_ok=True)
    with open(odir.with_suffix('.err').absolute(),'w') as err:
        for i in idir.iterdir():
            k=i.stem
            try:
                odict[k]=get_hex(i)
                genseq=fetch_seq(i)
                #TODO print -> logging
                print(f'> {k}\n{genseq}',file=open(odir/(k+':genome.fasta'),'w'))
                segs=get_valid_seg(hextranslate(genseq),name=k,thresh=thresh)
                # print(len(segs))
                print(seg_to_fasta(segs),file=open(odir/(k+':segs.fasta'),'w'))
            except Exception as e:
                err.write(f'{k}\t{e}\n')

def gen_sub_sys(idir:str,odir:str,
            scandir:str,
            workdir='.',
            tasknum=5,
            conda_dir="/home/hugheslab1/zfdeng/miniconda3",
            conda_env_name='interproscan',
            interproscan_bin="/home/hugheslab1/zfdeng/pangenome/interproscan/interproscan-5.65-97.0/interproscan.sh"):
    
    if workdir=='.':
        workdir=os.getcwd()
    head=f'''#!/usr/bin/env bash
    source "{conda_dir}/etc/profile.d/conda.sh"
    conda init
    conda activate {conda_env_name}

    workdir={workdir}
    cd {workdir}
    '''.replace("    ","")
    # single_line='''source {interproscan_bin} -appl Pfam,AntiFam,ProSiteProfiles  -i "genome_fasta/{stem}.fasta" -f JSON -b "scan_result/{stem}" -t {type} -cpu 4; cd $workdir'''
    single_line='''source {interproscan_bin} -appl Pfam,AntiFam,ProSiteProfiles  -i "{infile}" -f JSON -b "{ofile}" -t {type} -cpu 4; cd $workdir'''
    files=[]
    for i in range(tasknum):
        files.append([head])
    m=0
    idir:Path=Path(idir)
    odir:Path=Path(odir)
    scandir:Path=Path(scandir)
    odir.mkdir(exist_ok=True)
    scandir.mkdir(exist_ok=True)
    
    for fasta in idir.iterdir():
        stem=fasta.stem
        o=scandir/(stem+'.json')
        if not o.is_file():
            t='n' if stem.endswith('genome') else 'p'
            
            line=single_line.format(interproscan_bin=interproscan_bin,
                            infile=fasta.absolute(),
                            ofile=o.with_suffix('').absolute(),
                            type=t)
            _=(m%tasknum)
            files[_].append(line)
            m+=1
    for i in range(tasknum):
        print('\n'.join(files[i]),file=open(odir/f'subtask{i}.sh','w'))
        
    print(f"""for i in `ls {odir.absolute()}/*sh`
    do
        submitjob -w 24 -m 8 -c 4 source $i
    done""".replace("    ",""),
    file=open(odir.with_suffix('.sub.sh'),'w'))

# %%
def test():
    idir='../test/genbank_mata'
    odir='../test/genbank_fasta'
    fetch_fasta(idir,odir)
    
    gen_sub_sys(idir='../test/genbank_fasta',
            odir='../test/pbs_scripts',
            scandir='../test/scan_result')