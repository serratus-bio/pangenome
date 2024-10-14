from workflow.graph_construction_hhblits import hhblits_annotation as ha
import pandas as pd
from glob import glob
from pathlib import Path
annot_dir=Path('data/hhblits_pcs')
for infasta in Path('data/CoreData/curated_genome_fasta').iterdir():
    fasta_id=infasta.stem.replace(':genome','')
    inannotations=annot_dir/f'{fasta_id}.pkl'
    if inannotations.exists():
        annotations:pd.DataFrame=pd.read_pickle(inannotations)
    else:
        annotations = None
        
    ha.generate_neomodels(infasta,annotations)
    if annotations is not None:
        ha.connect_hit(annotations)
        ha.connect_hitregion(annotations)