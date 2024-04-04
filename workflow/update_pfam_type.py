# %%
from neomodel import config, db
from neomodel.integration.pandas import to_dataframe
from schema.models import *
import pandas as pd
config.DATABASE_URL = 'bolt://neo4j:WBrtpKCUW28e@44.206.130.87:7687'
import logging
from typing import Optional, Dict, List
import pickle as pkl

# %%


# %%
def update_pfam_type(pfam_type_dicts:List[dict]):
    q='''
            UNWIND $pfam_type_dicts as pfam_type_dict
            OPTIONAL MATCH (d:FuncDomainSet {accession:pfam_type_dict.accession})
            WITH d,pfam_type_dict
            WHERE d IS NOT NULL
            SET d.pfam_subtype = pfam_type_dict.pfam_type
            RETURN d
            '''
    t=to_dataframe(db.cypher_query(q,
            params={'pfam_type_dicts':pfam_type_dicts},
            resolve_objects=True))
    return t

if __name__=='__main__':
    q='''
        MATCH (domainset:FuncDomainSet)
        RETURN domainset.accession as accession
        '''
    accessions=to_dataframe(db.cypher_query(q,resolve_objects=True))['accession']
    pfam_type_dict:dict=pkl.load(open('../data/pfam_type_dict.pkl','rb'))
    pfam_type_dicts=[dict(accession=k,pfam_type=v) for k,v in pfam_type_dict.items() if k in accessions.to_list()]
    update_pfam_type(pfam_type_dicts)

# %%
