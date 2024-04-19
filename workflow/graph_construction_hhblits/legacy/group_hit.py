'''
legacy scratch for `group hit` procedure.
useful codes in this script has been integrated into `hhblits_annotation.py`
'''
# %%
import pickle as pkl
from typing import Union,Tuple,List
from Bio.Seq import Seq
from Bio import SeqIO
from pathlib import Path
import pandas as pd
import numpy as np
import logging
logger=logging.getLogger()
# %%
def proindice_to_genindice(genome:Union[str,Seq],prob:int,proe:int,transannot:str)->Tuple[int,int]:
    frame,direction=transannot
    frame=int(frame)
    b_,e_=frame+prob*3,frame+proe*3
    if direction=='s':
        return b_,e_
    else:
        return len(genome)-e_, len(genome)-b_

def genindice_to_proindice(genome:Union[str,Seq],genb:int,gene:int,transannot:str)->Tuple[int,int]:
    frame,direction=transannot
    frame=int(frame)
    
    if direction=='s':
        b_,e_=(genb-frame)//3,(gene-frame)//3
        return b_,e_
    else:
        return (len(genome)-frame-gene)//3,(len(genome)-frame-genb)//3

def genome_annot_indice(genome:Union[str,Seq],genb:int,gene:int,transannot:str)->Tuple[int,int]:
    'return: genseq segments '
    frame,direction=transannot
    # frame=int(frame)
    o=genome[genb:gene]
    o=o if direction=='s' else o[::-1]
    return o

def split_hit_line(hit_line:str)->List[str]:
    i=hit_line
    assert len(i)>94,f'invalid hit line: {i}'
    return [i[:4],i[4:34],i[34:40],i[40:48],i[48:56],i[56:63],
            i[63:69],i[69:74],i[74:85],i[85:94],i[94:]]
    
def parse_hitline(hit:List[str])->dict:
    #process domain accessions
    domain_accession,*_=hit[1].split(';')
    domain_accession=domain_accession.split('.')[0].strip()
    if len(_)>0:
        domain_annotation=_[0].strip()
    else:
        domain_annotation=domain_accession
    be = [int(i.strip())-1 for i in hit[8].split('-')]
    hmm_be = [int(i.strip())-1 for i in hit[9].split('-')]
    # process hit results
    hit_dict = dict(
    domain_accession=domain_accession,
    domain_annotation=domain_annotation,
    prob=float(hit[2].strip()),
    e=float(hit[3].strip()),
    p=float(hit[4].strip()),
    score=float(hit[5].strip()),
    ss=float(hit[6].strip()),
    begin=be[0],end=be[1],
    hmmbegin=hmm_be[0],hmmend=hmm_be[1],
    hmmlength=int(hit[10].strip().strip('(').strip(')'))
    )
    return hit_dict
    
def parse_segmeta(segmeta:Union[str,Path]):
    '''
    segmeta: the name of segfasta
        e.g. PosixPath('/tmp/tmpeazoz_z_/xx#0s#2069#3485.fasta')
    '''
    _,transannot,prob,proe=Path(segmeta).stem.split('#')
    meta_dict=dict(transannot=transannot,seg_b=int(prob),seg_e=int(proe))
    return meta_dict

from Bio.SeqRecord import SeqRecord
import numpy as np
from functools import partial

def parse_hhblits(blits_result:Tuple[Path,List[Tuple[Path,List[str]]]],
                  genome_dir:Path=Path('/home/hugheslab1/zfdeng/pangenome/CoreData/genome_fasta')):
    '''
    TODO refactor the genome fetching codes
    '''
    stem=blits_result[0].stem.replace(':segs','')
    genome_fasta_path=genome_dir/f'{stem}:genome.fasta'
    genome:SeqRecord=SeqIO.read(genome_fasta_path,'fasta')
    parses:List[dict]=[]
    for seg_res in blits_result[1]:
        if len(seg_res[1])>0:
            segmeta=parse_segmeta(seg_res[0])
            for h in seg_res[1]:
                hit=parse_hitline(h)
                parses.append(segmeta.copy())
                parses[-1].update(hit)
    hits_df=pd.DataFrame(parses)
    hits_df['gen_b'],hits_df['gen_e']=np.vectorize(partial(
        proindice_to_genindice,genome=genome))(
        prob=hits_df['begin']+hits_df['seg_b'],
        proe=hits_df['end']+hits_df['seg_b'],
        transannot=hits_df['transannot'])
    hits_df['gen_hmm_b'],hits_df['gen_hmm_e']=np.vectorize(partial(
        proindice_to_genindice,genome=genome))(
        prob=hits_df['begin']-hits_df['hmmbegin']+1+hits_df['seg_b'],
        proe=hits_df['end']+hits_df['hmmlength']-hits_df['hmmend']+hits_df['seg_b'],
        transannot=hits_df['transannot'])
    # hits_df['clan_accession'],hits_df['clan_annotation']=np.vectorize(query_clan)(
    #     hits_df['domain_accession'],hits_df['domain_annotation'])
    return genome,hits_df
                
annotation_file='../../data/blit_out.pkl'
annotations:List[Tuple[Path,List[Tuple[Path,List[str]]]]]=pkl.load(open(annotation_file,'rb'))
genome,hits_df=parse_hhblits(annotations[0])


def group_hit(hits_df:pd.DataFrame):
    o=[]
    for transannot,subg in hits_df.groupby('transannot'):
        # subg=hits_df[hits_df['transannot']=='1a'].copy(deep=True)
        subg.sort_values(by='gen_b',inplace=True)
        break_point=np.append(np.where(np.maximum.accumulate(subg.iloc[:-1]['gen_e'].to_numpy()
                    )-subg.iloc[1:]['gen_b'].to_numpy() < 0),len(subg)-1)
        differences = np.diff(break_point, prepend=-1).reshape(-1)
        subg['group_id'] = np.repeat(np.arange(differences.shape[-1]), differences)
        # subg['group_id'] = subg['transannot']+'-'+np.repeat(np.arange(differences.shape[-1]), differences).astype(str)
        subg['group_begin'] = subg.groupby('group_id')['gen_b'].transform('min')
        subg['group_end'] = subg.groupby('group_id')['gen_e'].transform('max')
        if ((subg.groupby('group_id')['group_end'].first().to_numpy()[1:
            ]-subg.groupby('group_id')['group_begin'].first().to_numpy()[:-1])>0).all():
            logger.error(f'abnormal grouping!: {hits_df}')
        o.append(subg)
    grouped_hits_df=pd.concat(o,axis=0)
    return grouped_hits_df
# %% 
if __name__=='__main__':
    genome,hits_df=parse_hhblits(annotations[1])
    grouped_hit=group_hit(hits_df)

# %%



