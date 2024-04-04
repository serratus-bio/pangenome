# %%
from pathlib import Path
import wget
import pyhmmer
from glob import glob
from tqdm import tqdm
import pickle as pkl
import pandas as pd
def get_mlength(h:str)->int:
    with pyhmmer.plan7.HMMFile(h) as hmm_file:
        hmm = hmm_file.read()
    return hmm.M

def pull_hmms():
    for i in tqdm([Path(i).stem.split('.')[0] for i in glob('/home/hugheslab1/zfdeng/pangengraph/pfam_self_compile/sto/*')]):
        wget.download(
                f'https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{i}?annotation=hmm',
                f'/home/hugheslab1/zfdeng/pangengraph/pfam/hmm/{i}.hmm.gz'
            )

def get_mlength_list():
    mlength_list=[]
    for i in glob("/home/hugheslab1/zfdeng/pangengraph/pfam/hmm/*gz"):#to_dataframe(db.cypher_query(''''''))['accession']:
        m=get_mlength(i)
        mlength_list.append(
            {'accession':Path(i).stem.split('.')[0],
            'std_length':m}
        )
    pkl.dump(mlength_list,open('pfam36_profile_length.pkl','wb'))

