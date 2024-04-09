from subprocess import run
from pathlib import Path
from tqdm import tqdm
from multiprocessing import Pool

def convert_sto(sto:Path):
    s=sto.stem
    o=Path('a3m')/f'{s}.a3m'
    r=run(['/home/hugheslab1/zfdeng/pangengraph/hh-suite/scripts/reformat.pl',
           "sto","a3m",sto,o],capture_output=True)
    return r

def split_sum_sto():
    '''
    DON'T RUN it directy! Note the file names
    sto file: from https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/
    '''
    with open('Pfam-A.seed','r',encoding='latin1') as f:
        i=0
        o=[]
        fn='holder'
        last_fn=fn
        for l in tqdm.tqdm(f.readlines()):
            if l.split()[-1].startswith('PF') and l.startswith('#=GF AC'):
                fn=l.split()[-1]
            if l=='//\n':
                if last_fn!=fn:
                    with open(f'sto/{fn}.sto','w') as fo:
                        fo.write(''.join(o))
                else:
                    with open(f'sto/warning-line{i}.sto','w') as fo:
                        fo.write(''.join(o))
                o=[]
                last_fn=fn
            else:
                o.append(l)
            i+=1
    stos=[i for i in Path('sto').iterdir()]
    pool=Pool(processes=32,maxtasksperchild=100)
    pool.map_async(convert_sto,stos)

    for i in Path('pfam_self_compile/a3m').iterdir():
        i.rename(i.with_stem(i.stem.split('.')[0]))
        
        