from neomodel import config, db
from neomodel.integration.pandas import to_dataframe
from schema.models import *
from typing import List
import logging
logger=logging.getLogger()
config.DATABASE_URL = 'bolt://neo4j:WBrtpKCUW28e@44.206.130.87:7687'

def update_hmm_length(mlength_list):
    '''
    mlength_list: pkl.load(open('pfam36_profile_length.pkl','rb'))
    '''
    
    q='''
    UNWIND $mlength_list as mlen
    OPTIONAL MATCH (s:FuncDomainSet {accession:mlen.accession})
    WITH s,mlen
    WHERE s IS NOT NULL
    SET s.std_length = mlen.std_length
    RETURN s
    '''
    t=to_dataframe(db.cypher_query(q,
        params={'mlength_list':mlength_list},
        resolve_objects=True))
    return t