from typing import List,Tuple
import logging
import sys, errno, re, json, ssl
from urllib import request
from urllib.error import HTTPError
from time import sleep
import networkx as nx
import pickle as pkl
logger=logging.getLogger()

def fetch_pfam(base_url:str)->List[Tuple[str,dict]]:
    context = ssl._create_unverified_context()
    attempts = 0
    output=[]
    next:str=base_url
    while next:
        try:
            req = request.Request(next, headers={"Accept": "application/json"})
            res = request.urlopen(req, context=context)
            # If the API times out due a long running query
            if res.status == 408:
                # wait just over a minute
                logger.warning('408: waiting')
                sleep(31)
                # then continue this loop with the same URL
                continue
            elif res.status == 204:
                logger.warning('204: break')
                #no data so leave loop
                break
            payload:dict = json.loads(res.read().decode())
            attempts = 0
            # if not next:
            #     last_page = True
        except HTTPError as e:
            if e.code == 408:
                sleep(61)
                continue
            else:
                # If there is a different HTTP error, it wil re-try 3 times before failing
                logger.warn(f'HTTPError: attempt = {attempts}')
                if attempts < 3:
                    attempts += 1
                    sleep(61)
                    continue
                else:
                    logger.error(f"HTTPError: LAST URL = {next}")
                    raise e
        output.append((next,payload))
        next = payload.get('next',None)
        if next:
            logger.log(20,f'next url: {next}')
            sleep(0.5)
    return output

def dump_clan():
    clans_payload=fetch_pfam('https://www.ebi.ac.uk:443/interpro/api/set/pfam')
    clan_catalog={}
    for _,payload in clans_payload:
        for metadata_ in payload['results']:
            metadata=metadata_['metadata']
            clan_catalog[metadata['accession']]={'name':metadata['name']}
    for i in clan_catalog.keys():
        clan_info=fetch_pfam(f'https://www.ebi.ac.uk:443/interpro/api/set/pfam/{i}')
        clan_catalog[i]['description']=clan_info[0][1]['metadata']['description']
        clan_catalog[i]['relationships']=clan_info[0][1]['metadata']['relationships']
        

    G=nx.DiGraph()
    for k,v in clan_catalog.items():
        G.add_node(k,
                label='clan',
                name=v['name'],
                description=v['description'])
        for n in v['relationships']['nodes']:
            if n['accession'] not in G:
                G.add_node(n['accession'],
                    label='entry',
                    name=n['short_name'],
                    description=v['name'])
            G.add_edge(k,n['accession'],
                    label='hasMember',
                    score=n['score'])
        for l in v['relationships']['links']:
            G.add_edge(l['source'],
                    l['target'],
                    label='hasLink',
                    score=l['score'])

    pkl.dump(G,open('clan_graph.pkl','wb'))