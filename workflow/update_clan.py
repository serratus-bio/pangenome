# %%
from neomodel import config, db
from neomodel.integration.pandas import to_dataframe
from schema.models import *
# from typing import List
# from tempfile import TemporaryDirectory
# from subprocess import run
from pathlib import Path
# import pandas as pd
config.DATABASE_URL = 'bolt://neo4j:WBrtpKCUW28e@44.206.130.87:7687'
import pickle as pkl
import networkx as nx
import logging

# %%
def update_clan_affiliation(clan_dict:dict,domain_accession:str,homology:float):
    '''
    `FuncDomainClan` could be created here,
    but `FuncDomainSet` must have existed.
    '''
    with db.transaction:
        clan=FuncDomainClan.get_or_create(clan_dict)[0]
        domainset=FuncDomainSet.nodes.get_or_none(accession=domain_accession)
        if domainset and clan:
            clan.hasclanentry.connect(domainset,properties=dict(linkscore=homology))
        else:
            logging.warning(f'invalid pair: {clan_dict} {domain_accession}')
    
def update_clanwise_homology(query_accession:str,target_accession:str,homology:float):
    with db.transaction:
        query:FuncDomainSet=FuncDomainSet.nodes.get_or_none(accession=query_accession)
        target:FuncDomainSet=FuncDomainSet.nodes.get_or_none(accession=target_accession)
        if query and target:
            query.hasclanhomology.connect(target,properties=dict(linkscore=homology))
        else:
            logging.warning(f'invalid pair: {query_accession} {target_accession}')


# %%
# accession='PF00869'

if __name__=='__main__':     
    clan_graph:nx.DiGraph=pkl.load(open('data/clan_graph.pkl','rb'))
    '''
    schema:
    {'label': 'hasMember', 'score': 0.004}
    {'label': 'hasLink', 'score': 0.0051}

    {'label': 'clan',
    'name': 'EGF',
    'description': 'Members of ...'}

    {'label': 'entry',
    'name': 'DSL', 'description': 'EGF'}
    '''
    from tqdm import tqdm
    q='''
            MATCH (domainset:FuncDomainSet)
            RETURN domainset.accession as accession
            '''
    accessions=to_dataframe(db.cypher_query(q,resolve_objects=True))['accession']
    for accession in tqdm(accessions):
        if accession not in clan_graph:
            logging.info(f'no clan: {accession}')
        else:
            for i in clan_graph.predecessors(accession):
                n=clan_graph.nodes[i]
                if n['label']=='clan':
                    clan_dict=dict(name=n['name'],source='Pfam',accession=i,annotation=n['description'])
                    homology=clan_graph[i][accession]['score']
                    update_clan_affiliation(clan_dict,accession,homology)
                elif n['label']=='entry':
                    homology=clan_graph[i][accession]['score']
                    update_clanwise_homology(i,accession,homology)
                    
            for i in clan_graph.neighbors(accession):
                homology=clan_graph[accession][i]['score']
                update_clanwise_homology(accession,i,homology)
                
    