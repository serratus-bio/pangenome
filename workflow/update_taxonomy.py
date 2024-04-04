# %%
from neomodel import config, db
from neomodel.integration.pandas import to_dataframe
from workflow.schema.models import *
import pandas as pd
config.DATABASE_URL = 'bolt://neo4j:WBrtpKCUW28e@44.206.130.87:7687'
import logging
from typing import Optional, Dict, List
import pickle as pkl

# %%
def robust_split(s:str,split:str)->List[str]:
    '''
    split and remove blankspace
    '''
    return [i.strip() for i in s.split(split)]
def get_genbank_id(accession:str)->Dict[str,str]:
    '''
    accession: entry in VMR
    '''
    # assert isinstance(accession,str),f'accession: {accession} is not a str!'
    if not isinstance(accession,str):
        logging.warning(f'accession: {accession} is not a str!')
        return {'_':'Null'}
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


# %%
def map_genid_taxonomy(entry:pd.Series):
    true_access:dict=entry['true_access']
    taxonomy:str=entry['taxonomy_str']
    odict={}
    for genid in true_access.values():
        if genid!='Null':
            odict[genid]=taxonomy
    return odict
def gen_genid_taxonomy_dict():
    virus_list=pd.read_csv('data/VMR_MSL38_v2.csv')
    # virus_list.columns
    virus_list['true_access']=virus_list['Virus GENBANK accession'].apply(get_genbank_id)
    virus_list['taxonomy_str']=virus_list[['Realm', 'Subrealm', 'Kingdom', 'Subkingdom',
        'Phylum', 'Subphylum', 'Class', 'Subclass', 'Order', 'Suborder',
        'Family', 'Subfamily', 'Genus', 'Subgenus', 'Species']].fillna('Null').apply(lambda x: ' ;'.join([str(i) for i in x]),axis=1)

    genid_taxonomy_dict={}
    for i in virus_list.apply(map_genid_taxonomy,axis=1):
        genid_taxonomy_dict.update(i)
    return genid_taxonomy_dict

# %%
def update_taxonomy(genid_taxonomy_dict:dict):
    genid_taxonomy_dicts=[dict(genid=k,taxo=v) for k,v in genid_taxonomy_dict.items()]
    q='''
        UNWIND $genid_taxonomy_dicts as genid_taxonomy
        OPTIONAL MATCH (g:Genome {accession:genid_taxonomy.genid})
        WITH g,genid_taxonomy
        WHERE g IS NOT NULL
        SET g.taxonomy = genid_taxonomy.taxo
        RETURN g
        '''
    t=to_dataframe(db.cypher_query(q,
            params={'genid_taxonomy_dicts':genid_taxonomy_dicts},
            resolve_objects=True))
    return t

# %%
if __name__=='__main__':
    genid_taxonomy_dict=gen_genid_taxonomy_dict()
    update_taxonomy(genid_taxonomy_dict)
    # pkl.dump(genid_taxonomy_dict,open('data/genid_taxonomy_dict.pkl','wb'))