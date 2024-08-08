# %%
from Bio import SeqIO,Seq
import pandas as pd
from subprocess import run
from tempfile import TemporaryDirectory
from pathlib import Path
import json
from typing import Tuple,List,Generator

from neomodel import config, db
from schema.models import *
import pandas as pd
from pathlib import Path
from Bio.Seq import Seq
from tqdm import tqdm
import logging 
from neomodel.integration.pandas import to_dataframe
logger=logging.getLogger()
config.DATABASE_URL = 'bolt://neo4j:YeLEtoRkhame@35.174.232.151:7687'
import sys
from glob import glob

logging.basicConfig(filename='dump.log', 
                    filemode='a',           
                    level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


# %%
def interproscan(infasta:str,
          fasta_type='n', # or 'q',
          conda_dir="~/miniconda3/envs/interproscan/",
          interproscan_bin='/home/hugheslab1/zfdeng/pangenome/interproscan/interproscan-5.65-97.0/interproscan.sh'
          )->dict:
     '''
     '''
     p_infasta=Path(infasta)
     with TemporaryDirectory() as tdir:
          cmd=["conda", "run", "-p", conda_dir,
               interproscan_bin,
               '-appl',
               'Pfam,AntiFam',
               '-i',
               p_infasta.absolute(),
               '-f',
               'JSON',
               '-b',
               f'{p_infasta.stem}',
               '-t',fasta_type,
               '-cpu','8']
          run(cmd,cwd=tdir,capture_output=True)
          genome_dict=json.load(open(f'{tdir}/{p_infasta.stem}.json'))
     return genome_dict
 
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
            id,sub_id)

def dump_single_genome(genome_dict:dict):
    '''
    draft version 
    very inefficient 
    `genome_dict`: from interproscan
    '''
    try:
        domains = parse_genome_dict(genome_dict)
        domains.sort_values(by='start',inplace=True)
        genome_seq=genome_dict['results'][0]['sequence']
        genome_name=genome_dict['results'][0]['crossReferences'][0]['name']
        genome_accession=genome_name.split('|')[-1]
        if Genome.nodes.get_or_none(name=genome_name) is None: 
       
            # domains.sort_values(by='start',inplace=True)
            genome:Genome=Genome.get_or_create(dict(
                    name=genome_name,
                    seq=genome_seq,
                    source='GenBank',
                    accession=genome_accession))[0]
            last_b,last_e,id,sub_id=0,0,0,1
            last_region:Optional[Region]=None
            last_domain_name='BEGIN'

            for _,s in domains.iterrows():
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
            logger.info(f'dumped: {genome_name}')
        else:
            logger.info(f'skipped: {genome_name}')
            pass
    except:
        logger.info(f'failed: {genome_name}')
        
# %%
if __name__=='__main__':
    # from tqdm import tqdm
    dumped=to_dataframe(db.cypher_query('''
        MATCH (g:Genome)
        RETURN g.name as genomes
        ''',resolve_objects=True))     
    valid_genome=[i.strip() for i in open('genome_classification/valid_genome','r').readlines()]

    for i in valid_genome:
        try:
            file=f'../../pangenome/CoreData/scan_result/{i}:genome.json'
            genome_dict=json.load(open(file))
        except:
            logger.info(f'not found: {i}')
        dump_single_genome(genome_dict)
    