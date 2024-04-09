# %%
from neomodel import config, db
from neomodel.integration.pandas import to_dataframe
from schema.models import *
from typing import List
config.DATABASE_URL = 'bolt://neo4j:WBrtpKCUW28e@44.206.130.87:7687'

def connect_domain_set():
    q='''
    MATCH (n:FuncDomainSet)-[x:hasDownstreamSet]->()-[:hasDownstreamSet]->(m:FuncDomainSet) 
    MERGE (n)-[y:hasNextDomainSet]->(m)
    SET y.degree=x.degree
    return n,m,y
    '''
    # MERGE (n)-[y:hasNextDomain {degree:x.degree}]->(m)
    return to_dataframe(db.cypher_query(q,resolve_objects=True))


def connect_domain():
    q='''
    MATCH (n:FuncDomain)-[:HasDownstream]->(:DomainLinkage)-[:HasDownstream]->(m:FuncDomain) 
    MERGE (n)-[y:hasNextDomain]->(m)
    SET y.begin_at=m.b-n.e
    WITH n,m,y,
    CASE 
        WHEN m.e > n.e THEN 0
        ELSE  1
    END as nested
    SET y.nested=nested
    RETURN n,y,m
    '''
    # MERGE (n)-[y:hasNextDomain {degree:x.degree}]->(m)
    return to_dataframe(db.cypher_query(q,resolve_objects=True))

