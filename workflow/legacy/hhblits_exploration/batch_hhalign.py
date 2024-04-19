'''
Inspect profiles
'''
import pandas as pd
from typing import Dict,List
from warnings import warn
virus_list=pd.read_csv('_data/VMR_MSL38_v2.csv')
from multiprocessing import Pool
from subprocess import run
from pathlib import Path

# ref_genome_dict={
#     'zika':'ZIKV||AY632535',
#     'ebola':'EBOV||AF086833',
#     'measles':'MeV||AB016162',
#     'covid':'SARS-CoV-2||MN908947'
# }
# virus_meta={'covid':('Order','Nidovirales'),
#             'zika':('Family','Flaviviridae'),
#             'ebola':('Family','Filoviridae'),
#             'measles':('Family','Paramyxoviridae'),
#             }

# #--- ---
# def is_multiple_access(access:Dict[str,str])->bool:
#     ''' 
#     helper for `get_genbank_id`'s output 
#     `False`: one seg; `True`:c segs 
#     '''
#     if '_' in access:
#         return False
#     else:
#         return True
    
# def get_file_stem(virus_series:pd.Series,
#     access_col='true_access',name_col='true_name')->List[str]:
#     '''
#     v_list['true_access']=nido_list['Virus GENBANK accession'].apply(get_genbank_id)
#     nido_list['true_name']=nido_list['Virus name abbreviation(s)'].apply(get_correct_name)
#     '''
#     _=virus_series
#     access:Dict[str,str]= _[access_col]
#     if not is_multiple_access(_[access_col]):
#         return [f"{_[name_col]}||{access['_']}"]
#     else:
#         return [f"{_[name_col]}|{k}|{v}" for k,v in access.items()]
    

# def robust_split(s:str,split:str)->List[str]:
#     '''
#     split and remove blankspace
#     '''
#     return [i.strip() for i in s.split(split)]

# def get_genbank_id(accession:str)->Dict[str,str]:
#     '''
#     accession: entry in VMR
#     '''
#     # TODO compatible with '/'
#     assert isinstance(accession,str),f'accession: {accession} is not a str!'
#     o={}
#     holder_token=0
#     if ';' in accession:
#         for sub_a in robust_split(accession,';'):
#             if ':' in sub_a:
#                 subk,subv=robust_split(sub_a,':')
#                 o[subk]=subv
#             else:
#                 o[f'{holder_token}']=sub_a
#                 holder_token+=1
#     else:
#         if ":" in accession:
#             subk,subv=robust_split(accession,':')
#             o[subk]=subv
#         else:
#             o['_']=accession
#     return o
#     # else:
#     #     return {'nan':'nan'}

# def get_correct_name(name:str)->str:
#     """
#     some names block contains multiple name sep by ';'
#     only use the first one
#     warning: this name is not unique!
#     """
#     assert isinstance(name,str),f'valid input: {name}'
#     if ';' in name:
#         warn(f'multiple names: {name},use the first one')
#         return robust_split(name,';')[0]
#     else:
#         return name
    
    
# # if 0:
# used_genome=[]
# for k,v in virus_meta.items():
#     sub_v_list:pd.DataFrame=virus_list[virus_list[v[0]]==v[1]].copy(deep=True)
#     sub_v_list['true_access']=sub_v_list['Virus GENBANK accession'].apply(get_genbank_id)
#     sub_v_list['true_name']=sub_v_list['Virus name abbreviation(s)'].apply(get_correct_name)

#     for _,s in sub_v_list.iterrows():
#         for stem in get_file_stem(s):
#             p=Path(f'_data/genome_fasta/{stem}:segs.fasta')
#             if p.is_file():
#                 # print(stem)
#                 used_genome.append(p)
#             else:
#                 print(p)

# #--- ---
# def parse_fasta_name(fa_na:str)->Dict[str,str]:
#     fa_na=fa_na[2:-1]
#     return {i.split('=')[0]:i.split('=')[1] for i in fa_na.split(',')}

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
    
# def mpi_scan_hhblits(infile:Path,dir:Path=Path('tmp'),processes=8,cpu:int=2):
#     pool1=Pool(processes=processes)
#     # c=0
#     f='placeholder'
#     res=[]
#     for i in open(infile,'r').readlines():
#         if i.startswith('>'):
#             if f!='placeholder':
#                 f.close()
#                 r=pool1.apply_async(hhblits,(dir/name_stem,))
#                 res.append(r)
#             # name_dict=parse_fasta_name(i)
#             # name_stem='#'.join(name_dict.values())+'.fasta'
#             f=open(dir/name_stem,'w')
#         if f!='placeholder':
#             f.write(i)
#     if f!='placeholder':
#         f.close()
#         r=pool1.apply_async(hhblits,(dir/name_stem,cpu))
#         res.append(r)
#     pool1.close()
#     pool1.join()
#     return [i.get() for i in res]

# def scan_hhblits(infile:Path,dir:Path=Path('tmp'),cpu:int=2):
#     # c=0
#     f='placeholder'
#     res=[]
#     for i in open(infile,'r').readlines():
#         if i.startswith('>'):
#             if f!='placeholder':
#                 f.close()
#                 r=hhblits(dir/name_stem,cpu)
#                 res.append(r)
#             name_dict=parse_fasta_name(i)
#             name_stem='#'.join(name_dict.values())+'.fasta'
#             f=open(dir/name_stem,'w')
#         if f!='placeholder':
#             f.write(i)
#     if f!='placeholder':
#         f.close()
#         r=hhblits(dir/name_stem,cpu)
#         res.append(r)
#     return res
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

# def tmp(i):
#     '''
#     scan used genomes
#     '''
#     with TemporaryDirectory() as t:  
#         tmpdir=Path(t)
#         return (i,scan_hhblits(i,tmpdir,4))
if 1:
    import tqdm
    from glob import glob
    import sys
    b,e=int(sys.argv[1]),int(sys.argv[2])
    all_hmm=[Path(i) for i in glob('/home/hugheslab1/zfdeng/pangengraph_2/hhs-a3m/*.a3m')][b:e]
    # print(len(all_segs))
    # print(len(used_genome))
    # used_genome=[i for i in all_segs if i not in used_genome]
    # print(len(used_genome))
    # import sys;sys.exit(0)
    # all_segs=[i for i in glob('_data/genome_fasta/*segs.fasta')]
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
    pkl.dump(_,open(f'hmm_align-{b}-{e}.pkl','wb'))
    
    