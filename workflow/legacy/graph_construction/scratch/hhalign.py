from subprocess import run
from pathlib import Path
from io import StringIO
from neomodel import db,config
import tqdm
from multiprocessing import Pool
config.DATABASE_URL = 'bolt://neo4j:WBrtpKCUW28e@44.206.130.87:7687'
from neomodel.integration.pandas import to_dataframe

q='''MATCH (n:FuncDomainSet)
    RETURN n.accession as family
    '''
all_sets=to_dataframe(db.cypher_query(q,resolve_objects=True))

i:str
j:str

from collections import defaultdict
from typing import Tuple


garbage=Path('null')
garbage.mkdir(exist_ok=True)
def hhalign(intuple:Tuple[str,str],a3m_dir=Path('poc/a3m'),):
    i,j=intuple
    r=run(['/home/hugheslab1/zfdeng/pangengraph/hh-suite/build/src/hhalign',
                '-i',(a3m_dir/f'{i}.a3m'), '-t',(a3m_dir/f'{j}.a3m'),
                '-o',garbage/f'{i}-{j}.hhr',
                '-hide_cons','-hide_pred',
                ], #'-show_ssconf'
                capture_output=True,
                )
    return (intuple,r)

# def tmp(i):
#     print(i)
    
od=defaultdict(dict)
intuples=[(i,j) for i in all_sets['family'] for j in all_sets['family'] if i!=j]


pbar=tqdm.tqdm(total=len(intuples))
def update_progress(ret):
    # print('mycallback is called with {}\n'.format(ret))
    # od[ret[0][0]][ret[0][1]]=ret[1]
    pbar.update()
    
if __name__ == '__main__':
    pool=Pool(processes=54,maxtasksperchild=100)
    res=[]

    for i in intuples:
        r=pool.apply_async(hhalign,(i,),callback=update_progress)
        res.append(r)
    pool.close()
    pool.join()
    # for r in res:
    #     r.wait()
import pickle as pkl
pkl.dump([i.get() for i in res],open('hhalign_result.pkl','wb'))