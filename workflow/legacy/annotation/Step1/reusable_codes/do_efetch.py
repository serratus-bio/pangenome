from Bio import SeqIO
from pathlib import Path
from Bio import Entrez
from warnings import warn
import xmltodict
import pickle as pkl
import tqdm
from pathlib import Path
import pandas as pd
from typing import List,Dict

Entrez.email='zfevan.deng@mail.utoronto.ca'

def fetch_parse(i:str):
    '''
    i: genbank accession
    '''
    xml_file=Entrez.efetch(db='nuccore',id=i,rettype="gb", retmode="xml")
    data_dict = xmltodict.parse(xml_file.read())
    return data_dict

def robust_split(s:str,split:str)->List[str]:
    '''
    split and remove blankspace
    '''
    return [i.strip() for i in s.split(split)]

def get_genbank_id(accession:str)->Dict[str,str]:
    '''
    accession: entry in VMR
    '''
    assert isinstance(accession,str),f'accession: {accession} is not a str!'
    o={}
    holder_token=0
    if ';' in accession:
        for sub_a in robust_split(accession,';'):
            if ':' in sub_a:
                subk,subv=robust_split(sub_a,':')
                o[subk]=subv
            else:
                o[f'{holder_token}']=sub_a
                holder_token+=1
    else:
        if ":" in accession:
            subk,subv=robust_split(accession,':')
            o[subk]=subv
        else:
            o['_']=accession
    return o

def get_correct_name(name:str)->str:
    """
    some names block contains multiple name sep by ';'
    only use the first one
    warning: this name is not unique!
    """
    assert isinstance(name,str),f'valid input: {name}'
    if ';' in name:
        warn(f'multiple names: {name},use the first one')
        return robust_split(name,';')[0]
    else:
        return name

def is_multiple_access(access:Dict[str,str])->bool:
    ''' 
    helper for `get_genbank_id`'s output 
    `False`: one seg; `True`:c segs 
    '''
    if '_' in access:
        return False
    else:
        return True
    
def get_file_stem(virus_series:pd.Series,
    access_col='true_access',name_col='true_name')->List[str]:
    '''
    v_list['true_access']=nido_list['Virus GENBANK accession'].apply(get_genbank_id)
    nido_list['true_name']=nido_list['Virus name abbreviation(s)'].apply(get_correct_name)
    '''
    _=virus_series
    access:Dict[str,str]= _[access_col]
    if not is_multiple_access(_[access_col]):
        return [f"{_[name_col]}||{access['_']}"]
    else:
        return [f"{_[name_col]}|{k}|{v}" for k,v in access.items()]
    