
# acan=pd.read_csv('zika/zika-domains-acan.csv')
from models import *
from neomodel import config, db
import pandas as pd
from pathlib import Path
from Bio.Seq import Seq
from tqdm import tqdm
import logging 
from neomodel.integration.pandas import to_dataframe
logger=logging.getLogger()
config.DATABASE_URL = 'bolt://neo4j:WBrtpKCUW28e@44.206.130.87:7687'
import sys
domains=pd.read_csv(sys.argv[1]) #'zika/zika-domains.csv'
def get_seq(k:str,fasta_dir='CoreData/genome_fasta'):
    genome_dir=Path(fasta_dir)
    seq=open(genome_dir/f"{k}:genome.fasta").readlines()[-1].strip()
    return seq

def get_rc_seq(seq:str):
    return Seq(seq).reverse_complement()._data.decode()

def get_props(last_b:int,last_e:int,b:int,e:int,id:int,sub_id:int):
    assert b>last_b
    if b>=last_e:
        begin_at,nested=0,0
        this_lk_id=float(id)
        this_dm_id=float(id+1)
        linkb=last_e 
        
        id+=2
        sub_id=1
        iselemental=1
        
    else:
        begin_at= b-last_e
        nested = 0 if e>last_e else 1
        this_lk_id=float(f"{id}.{sub_id}")
        this_dm_id=float(f"{id}.{sub_id+1}")
        linkb=b
        
        sub_id+=2
        iselemental=0
    return (begin_at,nested,
            this_lk_id,this_dm_id,
            linkb,iselemental,
            id,sub_id
            )

def dump_genome(g:str,subd:pd.DataFrame):
    with db.transaction:
        subd.sort_values(by='start',inplace=True)
        if Genome.nodes.get_or_none(name=g) is None: 
            # transaction indicates none or all: one exists == all exists
            genome:Genome=Genome.get_or_create(dict(
                    name=g,
                    seq=get_seq(g),
                    source='GenBank',
                    accession=g.split('|')[-1]))[0]
            last_b,last_e,id,sub_id=0,0,0,1
            last_region:Optional[Region]=None
            last_domain_name='BEGIN'
            
            for _,s in subd.iterrows():
                b,e=s['start'],s['end']
                this_domain_name,annot=s['domain_annotation'].split(':')
                (begin_at,nested,this_lk_id,this_dm_id,
                    linkb,iselemental,id,sub_id)=get_props(
                    last_b,last_e,b,e,id,sub_id)
                    
                linkname=f'Linkage:{last_domain_name}:{this_domain_name}'
                linkage=DomainLinkage.create(dict(
                    name=f'{s["genome_name"]}:{linkname}',
                    b=linkb,e=b))[0]
                linkage.connect_fasta_shortcut(genome,this_lk_id,iselemental)
                linkage.connect_regionset(dict(
                    name=linkname))
                linkage.connect_last_region(last_region,begin_at,nested)
                last_region=linkage

                funcdomain_name=f'Funcdomain:{this_domain_name}'
                funcdomain=HmmFuncDomain.get_or_create(dict(
                    name=f'{s["genome_name"]}:{funcdomain_name}',
                    b=b,hmmstart=s["hmmStart"],
                    e=e,hmmend=s["hmmEnd"]))[0]
                funcdomain.save()
                funcdomain.connect_fasta_shortcut(genome,this_dm_id,iselemental)
                funcdomain.connect_regionset(dict(
                            name=this_domain_name,source='Pfam',
                            annotation=annot, 
                            accession=s['domain_accession']))
                funcdomain.connect_last_region(last_region,begin_at,nested)
                
                last_region=funcdomain
                last_b,last_e=b,e
                last_domain_name=this_domain_name
                    
            this_domain_name='END'
            b,e=s['genome_length'],-1
            (begin_at,nested,this_lk_id,this_dm_id,
                    linkb,iselemental,id,sub_id)=get_props(
                    last_b,last_e,b,e,id,sub_id)
            linkname=f'Linkage:{last_domain_name}:{this_domain_name}'
            linkage=DomainLinkage.create(dict(
                name=f'{s["genome_name"]}:{linkname}',
                b=linkb,e=b))[0]
            linkage.connect_fasta_shortcut(genome,this_lk_id,iselemental)
            linkage.connect_regionset(dict(
                name=linkname))
            linkage.connect_last_region(last_region,begin_at,nested)
            last_region=linkage
        else:
            logger.info(f'processed entry. this function is only used to dump new genome entries from *domain.csv')

def connect_domain_set():
    q='''
    MATCH (n:FuncDomainSet)-[x:hasDownstreamSet]->()-[:hasDownstreamSet]->(m:FuncDomainSet) 
    MERGE (n)-[y:hasNextDomain]->(m)
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
# from typing import Union,Optional
for g,subd in tqdm(domains.groupby('genome_name')):
    dump_genome(g,subd)
connect_domain_set()
connect_domain()
    
# from multiprocessing import Pool
# def t(tp):
#     g,subd=tp
#     dump_genome(g,subd)
# pool=Pool(processes=8)
# pool.map(t,list(domains.groupby('genome_name'))[:2])

# MATCH (n:FuncDomainSet)-[x:hasDownstreamSet]->()-[:hasDownstreamSet]->(m:FuncDomainSet) MERGE (n)-[:NextDomain {degree:x.degree}]->(m)