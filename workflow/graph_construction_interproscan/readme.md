pipeline for the 1st version of neo4j graph.

## Overview:
### the schema:
the structure of neo4j graph.
Note this version is unannotated since it's on the brink of obsoletion

### the initial commit:
commit_genomes.py
    run interproscan for row fasta(or json pulled from genbank) and save the raw annotaion file(.csv).
    parse annotaion file into neomodel objects (nodes & connection) and commit them into neo4j instance. 

### the patches:
connect_domains.py
    update direct linkage between `FuncDomain` (`hasNextDomain`) and between `FuncDomainSet` (`hasDownstreamSet`).
    (in fact the information is redundant and would be better if calculate on-the-fly in downstream analysis)
    *note*: this version of submission only entailed unsegmented genomes.
update_clan.py
    update `FuncDomainClan` nodes and their connections with `FuncDomainSet`
update_hitscore.py
    something should be done with `commit_genomes`.
    update hitscore (generated during interproscan) between `Domain` and `DownstreamSet`
update_hmmlength.py
    something should be done with `commit_genomes`.
    update hmm model's length from `pfam36_profile_length.pkl` in data/
update_homology.py
    run mmseqs on every `FuncDomain` in each `FuncDomainSet` (on nt seqs instead of aa seqs)
    may have a bug (doesn't take strand into consideration)
update_pfam_type.py
    something should be done with `commit_genomes`.
    update hmm model's pfam_type (family, domain, motif,...) from `pfam_type_dict.pkl` in data/
update_taxonomy.py
    something should be done with `commit_genomes`.
    update taxomony from `pfam_type_dict.pkl` in data/
### the legacy:
legacy/annotate_genome.py: 
    useful codes are integrated into `genid_taxonomy_dict.py`.
    run large scale interproscan annotation on the cluster   