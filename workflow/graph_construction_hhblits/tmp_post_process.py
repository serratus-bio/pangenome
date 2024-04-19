'''
this is a temperary pipeline to update annotation pickles.
this scripts is currently only suitable for single thread process
and takes 10 hours to process 7384 annotated fasta.

note: it's supposed to be executed at main directory. kept here for neatness
be wary of of relative import (from .hhblits_annotation ) and file names('data/hhblits_annotation.log')
if you want to reproduce it
'''
# %%
from .hhblits_annotation import (
    hhblits_annotation,group_hit,generate_neomodels,
    connect_hit,connect_hitregion,Path)
import pandas as pd
import numpy as np
from neomodel import config,db
import logging
import sys
import pickle as pkl
from time import sleep
logging.basicConfig(filename='data/hhblits_commit.log', 
                    filemode='a',           
                    level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(module)s - %(filename)s - %(message)s')
logger=logging.getLogger()
# %%
# annotations=pd.read_pickle('data/hhblits_ori/ISAV|RNA8|AF404340.pkl')
# grouped_annotations=group_hit(annotations)
config.DATABASE_URL = 'bolt://neo4j:WBrtpKCUW28e@52.4.3.233:7687'
to_fasta=lambda x:f"/home/hugheslab1/zfdeng/pangenome/CoreData/genome_fasta/{x}:genome.fasta"
# %%

# print(file)
def main(file:str,threshold=0.8):
    stem=Path(file).stem
    logger.info(f'{stem}: commit START!')
    annotations:pd.DataFrame=pd.read_pickle(file)
    ori_len=len(annotations)
    if len(annotations)<1:
        logger.warning(f'{stem}: BREAK! no hit ""')
        return
    else:
        logger.info(f'{stem}: original hit counts: {ori_len}')
    
    if len(annotations[annotations['probab']>threshold])<1:
        max_probab=annotations['probab'].max()
        logger.warning(f'{stem}: BREAK! no valid hit, max_probab={max_probab:.2f}')
        return
    try:
        grouped_annotations=group_hit(annotations,threshold)
        filter_len=len(grouped_annotations)
        logger.info(f'{stem}: filtered hit count: {filter_len}')
        group_len=len(grouped_annotations['group_id'].unique())
        logger.info(f'{stem}: group_number: {group_len}')
    except Exception as e:
        logger.error(f'{stem}: Fail! group_hit failed as: {e}')
    
    try:
        with db.transaction:
            generate_neomodels(to_fasta(stem),grouped_annotations)
            connect_hit(grouped_annotations)
            connect_hitregion(grouped_annotations)
        logger.info(f'{stem}: FINISHED!')
        return 
    except Exception as e:
        logger.info(f'{stem}: initial commit failed with: {e}; retry after sleep.')
        sleep(5)
    try:  
        with db.transaction:
            generate_neomodels(to_fasta(stem),grouped_annotations)
            connect_hit(grouped_annotations)
            connect_hitregion(grouped_annotations)
        logger.info(f'{stem}: commit FINISHED!')
    except Exception as e:
        logger.error(f'{stem}: Fail! neomodel commit failed as: {e}')
        
if __name__=='__main__':
    # i=int(sys.argv[1])
    from tqdm import tqdm
    for file in tqdm(open('data/oris_list','r').readlines()):
        file=file.strip()
        main(file)
# from workflow.commitseq.interproscan_annotation import get_valid_seg,hextranslate
# from Bio.SeqIO import read
# to_fasta=lambda x:f"/home/hugheslab1/zfdeng/pangenome/CoreData/genome_fasta/{x}:genome.fasta"
# for i in get_valid_seg(hextranslate(read(
#     to_fasta('ISAV|RNA8|AF404340'),'fasta').seq)):
#     print(len(i[1]))
# %%
