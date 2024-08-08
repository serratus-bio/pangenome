pipeline for the 2st version of neo4j graph.
## Overview:
### the schema:
    the structure of neo4j graph.
    check docstring and `help_text`for description of nodes, relationships and properties 
### hhblits_annotation.py
    main pipeline,check its head string for details 
    TODO: 
        Genome instance (this time all fractured genomes are taken into consideration as well)
        Connection of HitFamily & Calculation of Hit-Hit similarity (considered to run on-the fly of analysis, instead during genome commit.)
### tmp_*.py
    tmp pipeline to run blocks in hhblits_annotation separately,
    due to the fact that neo4j commit is not so compatible with parallel running.