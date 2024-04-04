from pathlib import Path
import pandas as pd
from subprocess import run
from glob import glob
import json
from typing import Tuple,List,Generator

        
def iter_match(genome_dict:dict)->Generator[Tuple[dict,dict],None,None]:
    '''
    in a scan result json dict (of nt),
    iter through all orf-match pairs
    '''
    for orf in genome_dict['results'][0]['openReadingFrames']:
        for match in orf['protein']['matches']:
            yield (orf,match)
            
def get_domains(file_list:List[str],
                ostem='nido-domains'):
    '''
    get all PFAM hit from files matching the `file_limit` 
    and save to ${ostem}.csv
    
    entries of output scv: 
    genome_name,genome_length,domain_accession,strand,
    start,end,hmmStart,hmmEnd,evalue,domain_annotation
    
    accession: interproscan id
    '''
    o=[]
    for gfile in file_list:
        genome_name=Path(gfile).name.replace(":genome.json","")
        genome_dict=json.load(open(gfile))
        genome_length=len(genome_dict['results'][0]['sequence'])
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
                        loc['evalue'])
                    entry={}
                    entry['genome_name']=genome_name
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
                    entry['domain_annotation']=match_info[1]
                    o.append(entry)
    
    domains=pd.DataFrame(o)
    opath=Path(ostem).with_suffix('.csv')
    domains.to_csv(opath,index=False)
    
    acan_dict={}
    acan_dict['accession'],acan_dict['annotation']=[],[]
    for accession in domains['domain_accession'].unique():
        annot=domains[domains['domain_accession']==accession].iloc[0]['domain_annotation']
        acan_dict['accession'].append(accession)
        acan_dict['annotation'].append(annot)
        
    domain_annotations=pd.DataFrame(acan_dict)
    domain_annotations.to_csv(Path(ostem+'-acan').with_suffix('.csv'),index=False)

def test():
    get_domains(glob('../test/scan_result/*:genome.json'),'../test/domain')

