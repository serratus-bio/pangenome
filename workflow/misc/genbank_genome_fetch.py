import xmltodict
import pickle as pkl
from pathlib import Path
import pandas as pd
import warnings
from Bio import Entrez
from Bio.Seq import Seq
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)
from typing import Tuple,Dict,Union,List

# from typing import Union
Entrez.email='zfevan.deng@mail.utoronto.ca'
def fetch_parse(i:str):
    xml_file=Entrez.efetch(db='nuccore',id=i,rettype="gb", retmode="xml")
    data_dict = xmltodict.parse(xml_file.read())
    return data_dict

def fetch_genbank(odir:str,in_df:pd.DataFrame,
           get_name=lambda s:f'{s["Name"]}||{s["Genbank"]}',
           get_accession=lambda s:s["Genbank"]):
    odir:Path=Path(odir)
    odir.mkdir(exist_ok=True)
    o_dict={}
    with open(odir.with_suffix('.err').absolute(),'w') as err:
        for _,s in in_df.iterrows():
            try:
                k=get_name(s)
                ofile=odir/f'{k}.pkl'
                if not ofile.is_file():
                    data_dict=fetch_parse(get_accession(s))
                    o_dict[k]=data_dict
                    with open(ofile,'wb') as f:
                        pkl.dump(data_dict,f)
            except Exception as e:
                err.write(f'{k}\t{e}\n')
                
                
                
# %%
def test():
    infile='../test/2024.01.31-accession_list_JS.csv'
    vmr=pd.read_csv(infile)
    odir='../test/genbank_mata'
    fetch_genbank(odir,vmr)

