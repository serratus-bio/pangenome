'''
Temperal codes, should be merged into dump_genomes later.
'''
# %%
from neomodel import config, db
from neomodel.integration.pandas import to_dataframe
from schema.models import *
import pandas as pd
config.DATABASE_URL = 'bolt://neo4j:WBrtpKCUW28e@44.206.130.87:7687'
import logging
from typing import Tuple,Generator
import pickle as pkl
import json
from neomodel import NodeSet
from pathlib import Path
# %%
def iter_match(genome_dict:dict)->Generator[Tuple[dict,dict],None,None]:
    '''
    in a scan result json dict (of nt),
    iter through all orf-match pairs
    '''
    for orf in genome_dict['results'][0]['openReadingFrames']:
        for match in orf['protein']['matches']:
            yield (orf,match)

def parse_genome_dict(genome_dict:dict):
     genome_length=len(genome_dict['results'][0]['sequence'])
     annotations=[]
     for orf,match in iter_match(genome_dict):
          orf_info=orf['start'],orf['end'],orf['strand']
          match_info=(match['signature']['accession'],
                    f"{match['signature']['name']}:{match['signature']['description']}",
                    match['signature']['signatureLibraryRelease']['library'])
          if match_info[2]=='PFAM': #'PROSITE_PROFILES'
               for loc in match['locations']:
                    data=(loc['start'],
                         loc['end'],
                         loc['hmmStart'],
                         loc['hmmEnd'],
                         loc['evalue'],
                         loc['score'])
                    entry={}
                    entry['genome_name']=genome_dict['results'][0]['crossReferences'][0]['name']#SeqIO.read(infasta,'fasta').name
                    entry['genome_length']=genome_length
                    entry['domain_accession']=match_info[0]
                    
                    entry['strand']=orf_info[2]
                    if entry['strand']=='SENSE':
                         entry['start']=orf_info[0]+data[0]*3
                         entry['end']=orf_info[0]+data[1]*3
                    else:
                         entry['start']=orf_info[1]-data[1]*3
                         entry['end']=orf_info[1]-data[0]*3
                    entry['hmmStart']=data[2]
                    entry['hmmEnd']=data[3]
                    entry['evalue']=data[4]
                    entry['score']=data[5]
                    entry['domain_annotation']=match_info[1]
                    annotations.append(entry)
          
     domains=pd.DataFrame(annotations)
     return domains


# %%


def update_hitscore(annotation_json:str):
    domains=parse_genome_dict(json.load(open(annotation_json,'r')))
    for _,s in domains.iterrows():
        HmmFuncDomainnodes:NodeSet=HmmFuncDomain.nodes
        funcdomain_name=s['domain_annotation'].split(':')[0]
        funcdomain:HmmFuncDomain=HmmFuncDomainnodes.get_or_none(
                            name=f'{s["genome_name"]}:Funcdomain:{funcdomain_name}')
        if funcdomain:
            rel=funcdomain.memberof.relationship(funcdomain.memberof[0])
            rel.evalue=s['evalue']
            rel.score=s['score']
            rel.save()
        else:
            logging.warning(f'missing domain:{s["genome_name"]}:Funcdomain:{funcdomain_name}')


# %%     
if __name__=='__main__':
    from tqdm import tqdm
    logging.basicConfig(filename='../data/add_hitscores.log', 
                filemode='a',           
                level=logging.INFO,
                format='%(asctime)s - %(levelname)s - %(message)s')
    q='''
        MATCH (genome:Genome)
        RETURN genome.name as name
        '''
    names=to_dataframe(db.cypher_query(q,resolve_objects=True))['name']
    for name in tqdm(names):
        f=f'/home/hugheslab1/zfdeng/pangengraph/CoreData/scan_result/{name}:genome.json'
        if not Path(f).exists():
            logging.warn(f'no annotation file: {name}')
            continue
        try:
            update_hitscore(f)
            logging.info(f'updated: {name}')
        except:
            logging.warn(f'failed to update: {name}')
        
# %% 