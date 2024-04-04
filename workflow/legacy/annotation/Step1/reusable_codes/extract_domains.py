from glob import glob
from typing import Generator,Tuple
import json
def iter_match(genome_dict:dict)->Generator[Tuple[dict,dict],None,None]:
    '''
    in a scan result json dict (of nt),
    iter through all orf-match pairs
    '''
    for orf in genome_dict['results'][0]['openReadingFrames']:
        for match in orf['protein']['matches']:
            yield (orf,match)

o=[]
b=0

### seems to be ambiguous, polish needed

for gfile in glob('nido_subset/*:genome.json'):
    # print('\n\n##'+gfile.split('/')[1].replace(':genome.json',''))
    genome_name=gfile.split('/')[1].replace(':genome.json','')
    genome_dict=json.load(open(gfile))
    genome_length=len(genome_dict['results'][0]['sequence'])
    for orf,match in iter_match(genome_dict):
        orf_info=orf['start'],orf['end'],orf['strand']
    #     if orf_info[2]!='SENSE':
    #         b=1
    #         break
    # if b==1:
    #     break
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
                


_match=[]
'''
entry: genome-id, genome-begin, genome-end, strand, accession,hmm-start,hmm-end,evalue
accession: name:description
'''

        
for gfile in glob('nido_subset/*:genome.json'):
    print('\n\n##'+gfile.split('/')[1].replace(':genome.json',''))
    genome_dict=json.load(open(gfile))
    for orf in genome_dict['results'][0]['openReadingFrames']:
        print(orf['start'],orf['end'],orf['strand'])
        for match in orf['protein']['matches']:
            infos=(match['signature']['accession'],
                f"{match['signature']['name']}:{match['signature']['description']}",
                match['signature']['signatureLibraryRelease']['library'])
            if infos[0]=='PF00680':
                _match.append(match)
            if infos[2]=='PFAM': #'PROSITE_PROFILES'
                for loc in match['locations']:
                    data=(loc['start'],
                        loc['end'],
                        loc['hmmStart'],
                        loc['hmmEnd'],
                        loc['evalue'],
                        )#match['evalue']
                    # o.append()
                    print(f'\t {infos}',end='\t')
                    print(f'\t {data}')