import pickle as pkl
from typing import List
from Bio.Seq import Seq

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

def fetch_seq(f:str)->str:
    '''
    f: efetch result entry
    
    output: nt seq in the entry
    '''
    gb:dict=pkl.load(open(f,'rb'))
    seq:str=gb['GBSet']['GBSeq']['GBSeq_sequence']
    return seq


translate= lambda x:Seq(x).translate()._data.decode()

def hextranslate(g:str)->List[str]:
    '''
    genome: input nt seq
    output: list of sense trans*3 + antisense trans*3
    '''
    o=[]
    genome=Seq(g)
    genome_r:Seq=genome.reverse_complement()
    for i in [0,1,2]:
        o.append(genome[i:].translate()._data.decode())
        o.append(genome_r[i:].translate()._data.decode())
    return o

def get_hex(f:str)->str:
    return hextranslate(fetch_seq(f))


from typing import Tuple,List
from Bio.Seq import Seq
translate= lambda x:Seq(x).translate()._data.decode()
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
from typing import Tuple,Dict,Union
routine_dict=Dict[str,Union[str,int,float]]


transannot=['0s','0a','1s','1a','2s','2a']
transannot_indice={j:i for i,j in enumerate(transannot)}

def get_genome_length(translist:List[str])->int:
    return len(translist[0])

def get_valid_seg(translist:List[str],name:str='xx',thresh:int=100)->List[Tuple[routine_dict,str]]:
    '''
    split translist
    get head dict as :
    head={'name':name,
        'transannot':annot,
        'prob':b,
        'proe':e}
    and corresponding seqs
    '''
    output=[]
    for i,annot in zip(translist,transannot):
        b=e=0
        for seg in i.split('*'):
            e=b+len(seg)
            if len(seg)>thresh:
                head={'name':name,
                      'transannot':annot,
                      'prob':b,
                      'proe':e}
                output.append((head,seg))
            b=e+1
    return output


def scanindice_to_proindice(prob:int,scanb:int,scane:int):
    # TODO
    delta=0 
    return prob+delta+scanb,prob+delta+scane

def seg_to_fasta(seg:List[Tuple[routine_dict,str]])->str:
    '''
    seg: from `get_valid_seg`
    '''
    o=[]
    for s in seg:
       head='> ' + ','.join([f'{k}={v}' for k,v in s[0].items()])
       o.extend([head,s[1]])
    return '\n'.join(o)

def proindice_to_genindice(genome:str,prob:int,proe:int,transannot:str)->Tuple[int,int]:
    frame,direction=transannot
    frame=int(frame)
    b_,e_=frame+prob*3,frame+proe*3
    if direction=='s':
        return b_,e_
    else:
        return len(genome)-e_, len(genome)-b_

def genindice_to_proindice(genome:str,genb:int,gene:int,transannot:str)->Tuple[int,int]:
    frame,direction=transannot
    frame=int(frame)
    
    if direction=='s':
        b_,e_=(genb-frame)//3,(gene-frame)//3
        return b_,e_
    else:
        return (len(genome)-frame-gene)//3,(len(genome)-frame-genb)//3
        

def genome_annot_indice(genome:str,genb:int,gene:int,transannot:str)->Tuple[int,int]:
    'return: genseq segments '
    frame,direction=transannot
    # frame=int(frame)
    o=genome[genb:gene]
    o=o if direction=='s' else o[::-1]
    return o


# def hextranslate(genome:str)->List[str]:
#     o=[]
#     genome_r=genome[::-1]
#     for i in [0,1,2]:
#         o.append(translate(genome[i:]))
#         o.append(translate(genome_r[i:]))
#     return o

def test_seg():
    genome=fetch_seq('genbank_meta/AbV2||KY357507.pkl')
    t=hextranslate(genome)
    segs=get_valid_seg(t,thresh=100)
    for o,(head,seq) in enumerate(segs):
        prob,proe,annot=head['prob'],head['proe'],head['transannot']
        genb,gene=proindice_to_genindice(genome,prob,proe,annot)
        assert genindice_to_proindice(genome,genb,gene,annot)==(prob,proe),f'{o},{prob},{proe},{annot}'
        assert translate(genome_annot_indice(genome,genb,gene,annot))==seq,f'{o},{prob},{proe},{annot}'