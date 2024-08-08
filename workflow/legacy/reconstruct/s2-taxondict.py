import pandas as pd
import pickle as pkl
from workflow.fetchdata import vmr_parse as vp 

vmr=pd.read_csv('data/VMR_MSL38_v2.csv')
r_vmr=vmr[vmr['Genome composition'].apply(lambda x:'RNA' in x)].fillna('Null')
r_vmr['true_access']=r_vmr['Virus GENBANK accession'].apply(vp.get_genbank_id)
r_vmr['true_name']=r_vmr['Virus name abbreviation(s)'].apply(vp.get_correct_name)
taxo_dict={}
for v,k_ in zip(r_vmr[['Realm', 'Subrealm', 'Kingdom', 'Subkingdom',
              'Phylum', 'Subphylum', 'Class', 'Subclass', 'Order', 'Suborder',
              'Family', 'Subfamily', 'Genus', 'Subgenus', 'Species']].apply(lambda x:';'.join(x),axis=1),r_vmr['true_access']):
    for k in k_.values():
        if k!='Null':
            taxo_dict[k]=v
            
pkl.dump(taxo_dict,open('data/genid_taxonomy_dict.pkl','wb'))