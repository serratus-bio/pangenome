# from glob import glob

import pandas as pd
from subprocess import run,PIPE
from pathlib import Path

domains=pd.read_csv('nido-domains.csv')
domains.columns
for i in Path('cov19-hits').iterdir():
    stem=i.stem
    if '-cov2' not in stem and i.suffix=='.fasta':
        _=run(['./diamond',
                'makedb',
                '--in',
                str(i),
                '-d',
                str(i.with_name(f'{stem}-reference'))],
                stdout=PIPE,
                stderr=PIPE)
        _=run(['./diamond','blastp',
               '-d',str(i.with_name(f'{stem}-reference')),
                '-q',str(i.with_name(f'{stem}-cov2.fasta')),
                '-o',str(i.with_name(f'{stem}-match.tsv')),
                '--id', '0' ,
                '--max-target-seqs', '300', 
                '--header', 'verbose', 
                '--min-score', '0', 
                '--query-cover', '0', 
                '--subject-cover', '0', 
                '--evalue','1'],
                stdout=PIPE,
                stderr=PIPE)