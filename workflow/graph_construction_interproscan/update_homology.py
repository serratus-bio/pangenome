# %%
from neomodel import config, db
from neomodel.integration.pandas import to_dataframe
from schema.models import *
from typing import List
from tempfile import TemporaryDirectory
from subprocess import run
from pathlib import Path
import pandas as pd
config.DATABASE_URL = 'bolt://neo4j:WBrtpKCUW28e@44.206.130.87:7687'
# %%
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO import write
from typing import Optional
import logging
logging.basicConfig(filename='data/add_identity.log', 
                filemode='a',           
                level=logging.INFO,
                format='%(asctime)s - %(levelname)s - %(message)s')
def write_domain_seqs(accession:str,output:Optional[str]=None):
    q='''
    MATCH (region:FuncDomain)<-[r:HasMember]-(domainset:FuncDomainSet {accession:$accession})
    RETURN region
    '''
    output=output if output else accession+'.fasta'
    regions=to_dataframe(db.cypher_query(q,params={'accession':accession},resolve_objects=True))
    if len(regions)<=1:
        logging.warning(f'less than 2 regions:{accession}')
    regions['records']=regions['region'].apply(lambda x: SeqRecord(seq=Seq(x.seq),id=x.name,description=''))
    write(regions['records'],
        handle=open(output,'w'),format='fasta')
    return regions
    
#%%
def mmseq_distance_matrix(infasta:str,
    mmseqbin:str='~/miniconda3/envs/rdrp/bin/mmseqs',fo=lambda x:Path(x).absolute().with_suffix('.align.tsv'),
    thread:int=4):
    # fo=lambda x:Path(x).stem+'-align.tsv'
    with TemporaryDirectory() as temdir:
        oname=fo(infasta)
        infasta=Path(infasta).absolute()
        cmd=f'''
        {mmseqbin} createdb {infasta} db 
        {mmseqbin} createdb {infasta} query
        {mmseqbin} search query db aln tmp/ -a --search-type 3 --threads {thread}
        {mmseqbin} convertalis query db aln {oname} --format-mode 4
        '''
        print(cmd,file=open(f'{temdir}/scripts.sh','w'))
        r=run(['bash',f'{temdir}/scripts.sh'],cwd=temdir,capture_output=True)#{temdir}/ ,stdout=open(f'out','w'),stderr=open(f'err','w'),
        # print(r.stderr.decode())
        # print(r.stdout.decode())
        return pd.read_csv(oname,delimiter='\t')
    
# %%
def add_identity(distance_matrix:pd.DataFrame,ret=True)->Optional[pd.DataFrame]:
    distance_dicts=list(distance_matrix[distance_matrix['query']!=distance_matrix['target']][
        ['query','target','fident']].T.to_dict().values())
    if len(distance_dicts)>0:
        if ret:
            q='''
                UNWIND $distance_dicts as distance_dict
                MATCH (query:FuncDomain {name: distance_dict.query}),(target:FuncDomain {name: distance_dict.target})
                MERGE (query)-[identity:homologousTo {identity: distance_dict.fident}]->(target)
                RETURN query,identity,target
                '''
            identities=to_dataframe(db.cypher_query(q,params={'distance_dicts':distance_dicts},resolve_objects=True))
            return identities
        else:
            q='''
                UNWIND $distance_dicts as distance_dict
                MATCH (query:FuncDomain {name: distance_dict.query}),(target:FuncDomain {name: distance_dict.target})
                MERGE (query)-[identity:homologousTo {identity: distance_dict.fident}]->(target)
                '''
            db.cypher_query(q,params={'distance_dicts':distance_dicts})
            return 
          
#%%    
def add_identity_apoc(distance_matrix:pd.DataFrame,ret=True)->Optional[pd.DataFrame]:
    '''
    No Any Faster, but why?
    '''
    distance_dicts=list(distance_matrix[distance_matrix['query']!=distance_matrix['target']][
        ['query','target','fident']].T.to_dict().values())
    if len(distance_dicts)>0: 
        q='''
            CALL apoc.periodic.iterate(
            "UNWIND $distance_dicts AS distance_dict RETURN distance_dict",
            "MATCH (query:FuncDomain {name: distance_dict.query}),(target:FuncDomain {name: distance_dict.target})
             MERGE (query)-[identity:homologousTo {identity: distance_dict.fident}]->(target)
             RETURN query,identity,target",
            {parallel:False, params:{distance_dicts: $distance_dicts}}
            )
            '''
        identities=to_dataframe(db.cypher_query(q,params={'distance_dicts':distance_dicts},resolve_objects=True))
        return identities

#%%
if __name__=='__main__':
    import sys
    import tqdm
    logging.basicConfig(filename='data/add_identity.log', 
                filemode='w',           
                level=logging.INFO,
                format='%(asctime)s - %(levelname)s - %(message)s')
    logger=logging.getLogger()
    # for i in to_dataframe(db.cypher_query(q,resolve_objects=True))['accession'][:2]:
    q='''
        MATCH (domainset:FuncDomainSet)
        RETURN domainset.accession as accession
        '''
    accessions=to_dataframe(db.cypher_query(q,resolve_objects=True))['accession']
    for i in tqdm.tqdm(accessions):
        if 0:
            # SLOW! use cluster instead:
            #   for i in `ls tmp_fasta/*.fasta`; do 
            #   submitjob -w 24 -m 4 -c 4 conda run -p ~/miniconda3/envs/n4j python post_dump.py $i; done
            # distance_matrix=mmseq_distance_matrix(sys.argv[1])
            write_domain_seqs(accession=i,output=f'tmp_fasta/{i}.fasta')
            distance_matrix=mmseq_distance_matrix(f'tmp_fasta/{i}.fasta')
            
        if Path(f'tmp_fasta/{i}.align.tsv').exists():
            try:
                add_identity(pd.read_csv(f'tmp_fasta/{i}.align.tsv',delimiter='\t'))
                logging.info(f'`add_identity` finished for for {i}')
            except:
                logging.warning(f'`add_identity` fails for for {i}')

        else:
            logging.info(f'`mmseq_distance_matrix` has no output for for {i}')# %%
