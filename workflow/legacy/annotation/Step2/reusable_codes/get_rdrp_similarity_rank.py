'''
archieved as a reminder,
Please update Everything at Step1 
'''

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from glob import glob
import json
from typing import Generator,Tuple
import pandas as pd
from Bio import Align

def iter_match(genome_dict:dict)->Generator[Tuple[dict,dict],None,None]:
    '''
    same itermath as `extract_domains.py`
    put here for convienience
    '''
    for orf in genome_dict['results'][0]['openReadingFrames']:
        for match in orf['protein']['matches']:
            yield (orf,match)

def get_domains():
    '''
    cleansing needed
    domains comes from ???
    '''
    domains=pd.read_csv('nido-domains.csv')
    return domains
    
def get_fasta(accession:str='PF00680',file_limit='nido_subset/*:genome.json'):
    '''
    file_limit: restrict the file to search for the accession;
    accession: in each file, the accession to fetch;
    '''
    o=[]
    o_dict={}
    b=0
    for gfile in glob(file_limit):
        # print('\n\n##'+gfile.split('/')[1].replace(':genome.json',''))
        genome_name=gfile.split('/')[1].replace(':genome.json','')
        genome_dict=json.load(open(gfile))
        for orf,match in iter_match(genome_dict):
            if match['signature']['accession']==accession:
                seq=orf['protein']['sequence']
                if len(match['locations'])==1:
                    _=match['locations'][0]
                    b,e=_['start'],_['end']
                    rdrp_seq=seq[b:e]
                    rdrp_name=genome_name
                    o.append(f'>{rdrp_name}\n{rdrp_seq}')
                    o_dict[rdrp_name]=rdrp_seq
                else:
                    for i,_ in enumerate(match['locations']):
                        b,e=_['start'],_['end']
                        rdrp_seq=seq[b:e]
                        rdrp_name=genome_name+'#'+str(i)
                        o.append(f'>{rdrp_name}\n{rdrp_seq}')
                        o_dict[rdrp_name]=rdrp_seq
    # print('\n'.join(o),file=open('nido-rdrp.fasta','w'))  
    return o_dict,o

def write_fasta():
    '''
    cleansing needed
    '''
    o={}
    domains=get_domains()
    used_domain=domains[domains['genome_name']=='SARS-CoV-2||MN908947']['domain_accession']
    for i in used_domain:
        o[i]=get_fasta(i,'nido_subset/SARS-CoV-2||MN908947:genome.json')[0]
    for k,v in o.items():
        with open(f'cov19-hits/{k}-cov2.fasta','w') as f:
            fastas='\n'.join([f'>{k1}\n{v1}' for k1,v1 in v.items()])
            f.write(fastas)

'''
diamond usage:

(entry name for cov19:SARS-CoV-2||MN908947)

./diamond blastp -d nido-rdrp-reference \
-q sars_rdrp_hmm.fasta \
-o matches.tsv \
--id 0 \
--max-target-seqs 300 \
--header verbose \
--min-score 0 \
--query-cover 0 \
--subject-cover 0 \
--evalue 1 
'''

def parse_diamond(matchfile:str='diamond-matches.tsv')->pd.DataFrame:
    '''
    matchfile:
        output of diamond blastp
        with '--header', 'verbose', 
        and 'tsv' suffix
    '''
    head_lines=open(matchfile,'r').readlines()[2].strip().split(': ')[1]
    diamond_aligns=pd.read_csv(matchfile,skiprows=3,names=head_lines.split(', '),delim_whitespace=True)
    return diamond_aligns

def get_no_aligns():
    '''
    cleansing needed
    '''
    o_dict,o=get_fasta()
    diamond_aligns=parse_diamond()
    fail_to_match=[i for i in o_dict.keys() if i not in diamond_aligns['Subject ID'].to_list()]
    return fail_to_match
    
def get_no_aligns_sort_list():
    '''
    cleansing needed
    '''
    o_dict,o=get_fasta()
    fail_to_match=get_no_aligns()
    s_dict={}
    ref_seq=o_dict['SARS-CoV-2||MN908947']
    for i in fail_to_match:    
        query_seq=o_dict[i]
        aligner = Align.PairwiseAligner()
        alignments = aligner.align(query_seq, ref_seq)
        s_dict[i]=alignments[0].score
    
    no_aligns_sort=[(k,v) for k,v in s_dict.items()]
    no_aligns_sort.sort(key=lambda x:x[1],reverse=True)
    
    return no_aligns_sort

def write_list():
    '''
    cleansing needed
    '''
    diamond_aligns=parse_diamond()
    no_aligns_sort=get_no_aligns_sort_list()
    with open('sort-rdrp.list','w') as f:
        f.write('id,match\n')
        for i,j in zip(diamond_aligns['Subject ID'],diamond_aligns['Percentage of identical matches']):
            f.write(i+','+f'{j:.2f}'+'\n')
        for i in no_aligns_sort:
            f.write(i[0]+','+str(i[1])+'\n')  