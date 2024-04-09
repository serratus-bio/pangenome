# %%
from glob import glob
from pathlib import Path
import logging
import pickle as pkl

# %%
def get_pfam_type(sto_file:str):
    'sto file: from https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.seed.gz'
    for l in open(sto_file,'r') .readlines():
        if l.startswith('#=GF TP'):
            return l.split()[-1]
    else:
        logging.warn(f'no TP line found: {sto_file}')
        return 'NULL'
    
def dump_pfam_type():
    pfam_type_dict={}
    for i in glob('/home/hugheslab1/zfdeng/pangengraph/pfam_self_compile/sto/*.sto'):
        stem=Path(i).stem.split('.')[0]
        pfam_type_dict[stem]=get_pfam_type(i)
    pkl.dump(pfam_type_dict,open('data/pfam_type_dict.pkl','wb'))
    
# %%
