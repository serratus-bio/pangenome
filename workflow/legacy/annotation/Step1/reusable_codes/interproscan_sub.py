from pathlib import Path
head='''#!/usr/bin/env bash
source "/home/hugheslab1/zfdeng/miniconda3/etc/profile.d/conda.sh"
conda init
conda activate interproscan

workdir=`pwd`
'''
workdir='/home/hugheslab1/zfdeng/pangenome'
single_line='''source $interproscan -appl Pfam,AntiFam,ProSiteProfiles  -i "genome_fasta/{stem}.fasta" -f JSON -b "scan_result/{stem}" -t {type} -cpu 4; cd $workdir'''

odir=Path('cluster_scripts/continue_run')

tasknum=40
files=[]
for i in range(tasknum):
    files.append([head])
m=0
for meta in Path('genome_fasta/').iterdir():
    stem=meta.stem
    o=Path('scan_result')/(stem+'.json')
    if not o.is_file():
        t='n' if stem.endswith('genome') else 'p'
        line=single_line.format(stem=stem,type=t)
        _=(m%tasknum)
        files[_].append(line)
        m+=1
for i in range(tasknum):
    print('\n'.join(files[i]),file=open((odir/f'subtask{i}.sh').absolute(),'w'))