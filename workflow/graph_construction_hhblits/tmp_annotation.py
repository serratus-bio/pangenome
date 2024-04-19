'''
this is a temperary pipeline to run hhblits annotation and dump the output into pkl files.
this script is safe to run parallel (multiple processes in bash / multiple threads by Pool)
adjust `main`'s `cpu` to allocate resources to single "Fasta"'s annotations

note: it's supposed to be executed at main directory. kept here for neatness
be wary of of relative import (from .hhblits_annotation ) and file names('data/hhblits_annotation.log')
if you want to reproduce it
'''
# %%
# from workflow.commitseq.hhblits_annotation import (
#     hhblits_annotation,group_hit,generate_neomodels,
#     Path)
from .hhblits_annotation import (
    hhblits_annotation,group_hit,generate_neomodels,
    Path)
import logging
logger=logging.getLogger()
logging.basicConfig(filename='data/hhblits_annotation.log', 
                    filemode='a',           
                    level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(module)s - %(filename)s - %(message)s')
def main(infasta:str,cpu:int=4,
         ori_dir=Path('data/hhblits_ori'),
         pcs_dir=Path('data/hhblits_pcs')):
    stem=Path(infasta).stem.replace(':genome','')
    
    try:
        annotations=hhblits_annotation(infasta,cpu)
        annotations.to_pickle(ori_dir/f'{stem}.pkl')
    #     if len(annotations)<1:
    #         logger.warning(f'BREAK! no hit for "{infasta}"')
    #         return
    #     if len(annotations[annotations['probab']>0.8])<1:
    #         logger.warning(f'no valid hit for "{infasta}"')
    #         max_probab=annotations['probab'].max()
    #         logger.warning(f'max_probab: "{max_probab}"')
    #         return
    #     grouped_annotations=group_hit(annotations)
    #     grouped_annotations.to_pickle(pcs_dir/f'{stem}.pkl')
    except Exception as e:
        logger.error(f'{stem} failed with error: {e}')
    # generate_neomodels(infasta,grouped_annotations)
    
# %%
if __name__=='__main__':
    import sys
    # main(sys.argv[1])
    from multiprocessing import Pool
    from tqdm import tqdm
    from glob import glob
    fastas=glob("/home/hugheslab1/zfdeng/pangenome/CoreData/genome_fasta/*:genome.fasta")
    qbar=tqdm(total=len(fastas))
    def callback(i):
        qbar.update()
    pool=Pool(processes=5,maxtasksperchild=500)
    for i in fastas:
        pool.apply_async(main,(i,),callback=callback)
    pool.close()
    pool.join()
    
# grouped_annotations=main('/home/hugheslab1/zfdeng/pangengraph_2/_data/genome_fasta/ScVLA||J04692:genome.fasta')
    # grouped_annotations
# to_fasta=lambda x:f"/home/hugheslab1/zfdeng/pangenome/CoreData/genome_fasta/{x}:genome.fasta"
# main(to_fasta('SNV|M|L37903'))
# %%