import pandas as pd
from typing import Dict,List
from warnings import warn
virus_list=pd.read_csv('_data/VMR_MSL38_v2.csv')
from multiprocessing import Pool
from subprocess import run
from pathlib import Path

def hhblits(infile:Path,cpu:int=2,threshold:int=50):
    # infile='tmp/AAbV||GBBW01007738#2s#11629#11886.fasta'
    o=run(['hhblits','-i',infile,'-d','hhs-db/pfam','-cpu',str(cpu)],capture_output=True)#,'-o','hahaha.hhr'
    skip_tag=1
    output=[]
    for i in o.stdout.decode().split('\n'):
        if skip_tag:
            if i.startswith(' No Hit'):
                skip_tag=0
                # print(','.join([i[:4],i[4:34],i[34:40],i[40:48],i[48:56],i[56:63],i[63:69],i[69:74],i[74:85],i[85:94],i[94:]]))
        else:
            if len(i)>90 and float(i[34:40].strip())>threshold:
                # print(','.join([i[:4],i[4:34],i[34:40],i[40:48],i[48:56],i[56:63],i[63:69],i[69:74],i[74:85],i[85:94],i[94:]]))
                output.append([i[:4],i[4:34],i[34:40],i[40:48],i[48:56],i[56:63],i[63:69],i[69:74],i[74:85],i[85:94],i[94:]])
    return (infile,output)

from tempfile import TemporaryDirectory
import os
from warnings import warn
def scan_hhalign(infile:Path,cpu:int=2):
    infile,output=hhblits(infile,cpu,0)
    if infile.with_suffix('.hhr').is_file():
        os.remove(infile.with_suffix('.hhr'))
    else:
        warn(f'Failed run:{infile}')
    return infile.stem,output


if 1:
    import tqdm
    from glob import glob
    import sys
    import pickle as pkl

    all_hmm=[Path(f'/home/hugheslab1/zfdeng/pangengraph_2/hhs-a3m/{i}.a3m') for i in pkl.load(open('rna_virus_domains.pkl','rb'))] #demo_4_domains

    qbar=tqdm.tqdm(total=len(all_hmm))
    def callback(i):
        qbar.update()
    
    res=[]
    # infile=used_genome[0]
    pool=Pool(processes=11,maxtasksperchild=100)
    for infile in all_hmm:
        r=pool.apply_async(scan_hhalign,(infile,),callback=callback)
        res.append(r)
    pool.close()
    pool.join() 
    
    import pickle as pkl
    _=[i.get() for i in res]
    pkl.dump(_,open(f'hmm_align-rna_virus_domains.pkl','wb'))
    
    