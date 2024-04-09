# %%
from neomodel import config, db
from neomodel.integration.pandas import to_dataframe
from workflow.schema.models import *
from typing import List
from tempfile import TemporaryDirectory
from subprocess import run
from pathlib import Path
import pandas as pd
config.DATABASE_URL = 'bolt://neo4j:WBrtpKCUW28e@44.206.130.87:7687'
import logging
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO import write
from typing import Optional
import pickle as pkl

# logging.LogRecord

# %%
# def write_domain_seqs(accession:str,output:Optional[str]=None):
#     q='''
#     MATCH (region:FuncDomain)<-[r:HasMember]-(domainset:FuncDomainSet {accession:$accession})
#     RETURN region
#     '''
#     output=output if output else accession+'.fasta'
#     regions=to_dataframe(db.cypher_query(q,params={'accession':accession},resolve_objects=True))
#     regions['records']=regions['region'].apply(lambda x: SeqRecord(seq=Seq(x.seq),id=x.name,description=''))
#     write(regions['records'],
#         handle=open(output,'w'),format='fasta')
#     return regions
    
    
# #%%
# def mmseq_distance_matrix(infasta:str,
#     mmseqbin:str='~/miniconda3/envs/rdrp/bin/mmseqs',fo=lambda x:Path(x).absolute().with_suffix('.align.tsv'),
#     thread:int=4):
#     # fo=lambda x:Path(x).stem+'-align.tsv'
#     with TemporaryDirectory() as temdir:
#         oname=fo(infasta)
#         infasta=Path(infasta).absolute()
#         cmd=f'''
#         {mmseqbin} createdb {infasta} db 
#         {mmseqbin} createdb {infasta} query
#         {mmseqbin} search query db aln tmp/ -a --search-type 3 --threads {thread}
#         {mmseqbin} convertalis query db aln {oname} --format-mode 4
#         '''
#         print(cmd,file=open(f'{temdir}/scripts.sh','w'))
#         r=run(['bash',f'{temdir}/scripts.sh'],cwd=temdir,capture_output=True)#{temdir}/ ,stdout=open(f'out','w'),stderr=open(f'err','w'),
#         # print(r.stderr.decode())
#         # print(r.stdout.decode())
#         return pd.read_csv(oname,delimiter='\t')
    
# # distance_matrix=mmseq_distance_matrix('PF00271.fasta')
# # %%
# def add_identity(distance_matrix:pd.DataFrame,ret=True)->Optional[pd.DataFrame]:
#     distance_dicts=list(distance_matrix[distance_matrix['query']!=distance_matrix['target']][
#         ['query','target','fident']].T.to_dict().values())
#     if len(distance_dicts)>0:
#         if ret:
#             q='''
#                 UNWIND $distance_dicts as distance_dict
#                 MATCH (query:FuncDomain {name: distance_dict.query}),(target:FuncDomain {name: distance_dict.target})
#                 MERGE (query)-[identity:homologousTo {identity: distance_dict.fident}]->(target)
#                 RETURN query,identity,target
#                 '''
#             identities=to_dataframe(db.cypher_query(q,params={'distance_dicts':distance_dicts},resolve_objects=True))
#             return identities
#         else:
#             q='''
#                 UNWIND $distance_dicts as distance_dict
#                 MATCH (query:FuncDomain {name: distance_dict.query}),(target:FuncDomain {name: distance_dict.target})
#                 MERGE (query)-[identity:homologousTo {identity: distance_dict.fident}]->(target)
#                 '''
#             db.cypher_query(q,params={'distance_dicts':distance_dicts})
#             return 
          
# #%%    
# def add_identity_apoc(distance_matrix:pd.DataFrame,ret=True)->Optional[pd.DataFrame]:
#     '''
#     No Any Faster, but why?
#     '''
#     distance_dicts=list(distance_matrix[distance_matrix['query']!=distance_matrix['target']][
#         ['query','target','fident']].T.to_dict().values())
#     if len(distance_dicts)>0: 
#         q='''
#             CALL apoc.periodic.iterate(
#             "UNWIND $distance_dicts AS distance_dict RETURN distance_dict",
#             "MATCH (query:FuncDomain {name: distance_dict.query}),(target:FuncDomain {name: distance_dict.target})
#              MERGE (query)-[identity:homologousTo {identity: distance_dict.fident}]->(target)
#              RETURN query,identity,target",
#             {parallel:False, params:{distance_dicts: $distance_dicts}}
#             )
#             '''
#         identities=to_dataframe(db.cypher_query(q,params={'distance_dicts':distance_dicts},resolve_objects=True))
#         return identities
# i='PF00869'
# add_identity_apoc(pd.read_csv(f'tmp_fasta/{i}.align.tsv',delimiter='\t'))
#%%
# if __name__=='__main__':
#     import sys
#     import tqdm
#     logging.basicConfig(filename='data/add_identity.log', 
#                 filemode='w',           
#                 level=logging.INFO,
#                 format='%(asctime)s - %(levelname)s - %(message)s')
#     logger=logging.getLogger()
#     # for i in to_dataframe(db.cypher_query(q,resolve_objects=True))['accession'][:2]:
#     q='''
#         MATCH (domainset:FuncDomainSet)
#         RETURN domainset.accession as accession
#         '''
#     accessions=to_dataframe(db.cypher_query(q,resolve_objects=True))['accession']
#     for i in tqdm.tqdm(accessions):
#         if 0:
#             # SLOW! use cluster instead:
#             #   for i in `ls tmp_fasta/*.fasta`; do 
#             #   submitjob -w 24 -m 4 -c 4 conda run -p ~/miniconda3/envs/n4j python post_dump.py $i; done
#             # distance_matrix=mmseq_distance_matrix(sys.argv[1])
#             write_domain_seqs(accession=i,output=f'tmp_fasta/{i}.fasta')
#             distance_matrix=mmseq_distance_matrix(f'tmp_fasta/{i}.fasta')
            
#         if Path(f'tmp_fasta/{i}.align.tsv').exists():
#             try:
#                 add_identity(pd.read_csv(f'tmp_fasta/{i}.align.tsv',delimiter='\t'))
#                 logging.info(f'`add_identity` finished for for {i}')
#             except:
#                 logging.warning(f'`add_identity` fails for for {i}')
#             # try:
#             #     add_identity(pd.read_csv(f'tmp_fasta/{i}.align.tsv',delimiter='\t'))
#             #     logging.warning(f'`add_identity` fails for for {i}')
#             # except:
#             #     logging.warning(f'`add_identity` fails for for {i}')
#         else:
#             logging.info(f'`mmseq_distance_matrix` has no output for for {i}')
            
# import sys;sys.exit(0)


              
# get_type('/home/hugheslab1/zfdeng/pangengraph/pfam_self_compile/sto/PF00003.26.sto')      
# %%

# def robust_split(s:str,split:str)->List[str]:
#     '''
#     split and remove blankspace
#     '''
#     return [i.strip() for i in s.split(split)]
# def get_genbank_id(accession:str)->Dict[str,str]:
#     '''
#     accession: entry in VMR
#     '''
#     # assert isinstance(accession,str),f'accession: {accession} is not a str!'
#     if not isinstance(accession,str):
#         logging.warning(f'accession: {accession} is not a str!')
#         return {'_':'Null'}
#     o={}
#     holder_token=0
#     if ';' in accession:
#         for sub_a in robust_split(accession,';'):
#             if ':' in sub_a:
#                 subk,subv=robust_split(sub_a,':')
#                 o[subk]=subv
#             else:
#                 o[f'{holder_token}']=sub_a
#                 holder_token+=1
#     else:
#         if ":" in accession:
#             subk,subv=robust_split(accession,':')
#             o[subk]=subv
#         else:
#             o['_']=accession
#     return o


# # %%
# def map_genid_taxonomy(entry:pd.Series):
#     true_access:dict=entry['true_access']
#     taxonomy:str=entry['taxonomy_str']
#     odict={}
#     for genid in true_access.values():
#         if genid!='Null':
#             odict[genid]=taxonomy
#     return odict
# def gen_genid_taxonomy_dict():
#     virus_list=pd.read_csv('data/VMR_MSL38_v2.csv')
#     # virus_list.columns
#     virus_list['true_access']=virus_list['Virus GENBANK accession'].apply(get_genbank_id)
#     virus_list['taxonomy_str']=virus_list[['Realm', 'Subrealm', 'Kingdom', 'Subkingdom',
#         'Phylum', 'Subphylum', 'Class', 'Subclass', 'Order', 'Suborder',
#         'Family', 'Subfamily', 'Genus', 'Subgenus', 'Species']].fillna('Null').apply(lambda x: ' ;'.join([str(i) for i in x]),axis=1)

#     genid_taxonomy_dict={}
#     for i in virus_list.apply(map_genid_taxonomy,axis=1):
#         genid_taxonomy_dict.update(i)
#     return genid_taxonomy_dict

# genid_taxonomy_dict=gen_genid_taxonomy_dict()
# pkl.dump(genid_taxonomy_dict,open('data/genid_taxonomy_dict.pkl','wb'))
# # %%
# def update_taxonomy(genid_taxonomy_dict:dict):
#     genid_taxonomy_dicts=[dict(genid=k,taxo=v) for k,v in genid_taxonomy_dict.items()]
#     q='''
#         UNWIND $genid_taxonomy_dicts as genid_taxonomy
#         OPTIONAL MATCH (g:Genome {accession:genid_taxonomy.genid})
#         WITH g,genid_taxonomy
#         WHERE g IS NOT NULL
#         SET g.taxonomy = genid_taxonomy.taxo
#         RETURN g
#         '''
#     t=to_dataframe(db.cypher_query(q,
#             params={'genid_taxonomy_dicts':genid_taxonomy_dicts},
#             resolve_objects=True))
#     return t