'''
check docstring and `help_text`for description of nodes, relationships and properties 
'''
# %%
from __future__ import annotations
from neomodel import (StructuredNode, StructuredRel,
                    OneOrMore, ZeroOrOne,ZeroOrMore,One,
                    RelationshipTo,RelationshipFrom,
                    StringProperty, FloatProperty,
                    IntegerProperty,BooleanProperty)
from neomodel import config, db
from neomodel import clear_neo4j_database
import pandas as pd
from typing import Optional

# %% 
TAXO_ORDER=['Realm', 'Subrealm', 'Kingdom', 'Subkingdom',
    'Phylum', 'Subphylum', 'Class', 'Subclass', 'Order', 'Suborder',
    'Family', 'Subfamily', 'Genus', 'Subgenus', 'Species']

choices=lambda x:{i:'' for i in set(x)}

# %% 

class hasHit(StructuredRel):
    '''
    Between Fasta and Hit
    '''
    pass

class hasRegion(StructuredRel):
    '''
    Between Fasta and HitRegion
    '''
    regid=IntegerProperty(required=True)
    # pass

class hasDownstream(StructuredRel):
    '''
    Between HitRegion and HitRegion
    '''
    pass
     
class hasMember(StructuredRel):
    '''
    Between HitFamily and Hit
    '''
    template_neff=FloatProperty(required=True)
    probab=FloatProperty(required=True)
    identities=FloatProperty(required=True)
    similarity=FloatProperty(required=True)
    e_value=FloatProperty(required=True)
    p_value=FloatProperty(required=True)
    aligned_cols=FloatProperty(required=True)
    score=FloatProperty(required=True)
    sum_probs=FloatProperty(required=True)
    

class hasAffiliate(StructuredRel):
    '''
    Between HitFamily and Hit
    '''
    representation=BooleanProperty(required=True)
    sum_probs=FloatProperty(required=True)
    coverage=FloatProperty(required=True)

# class hasSynteny(StructuredRel):
#     pass

class hasAnalog(StructuredRel):
    '''
    Between Hit and Hit. 
    if two Hits in one HitFamily has alignment score beyond a certain threshold,
    there will be a hasAnalog rel between them.
    
    to be implemented.
    '''
    pass

class hasFasta(StructuredRel):
    '''
    Link Fasta of the same species to the species' genome.
    '''
    pass
# %%
class Fasta(StructuredNode):
    '''
    Root node of one nt sequence 
    '''
    name = StringProperty(required=True,unique_index=True)
    source = StringProperty(required=True)
    seq = StringProperty(required=True)
    accession = StringProperty(required=True)
    annotation = StringProperty()
    circular = BooleanProperty(default=False,help_text='if the fasta is a circular fasta')
    taxonomy = StringProperty(help_text=f'semicolon separated taxonomy, order:{TAXO_ORDER}')
    hits = RelationshipTo("Hit", "hasHit", cardinality=ZeroOrMore, model=hasHit)
    regions = RelationshipTo("HitRegion", "hasRegion", cardinality=OneOrMore, model=hasRegion)


class Genome(StructuredNode):
    '''
    To be implemented
    a set of fasta to define the genome of a species
    Not Impletemented
    '''
    __abstract_node__ = True
    #hasFasta


class Segment(StructuredNode):
    '''
    abstract, not shown in neo4j graph
    a certain part on the genome
    '''
    __abstract_node__ = True
    name = StringProperty(required=True,unique_index=True)
    begin = IntegerProperty(required=True)
    end = IntegerProperty(required=True)
    

class HitRegion(Segment):
    '''
    unoverlapped annotation region
    each HitRegion might have one or more Hit,
    ideally these hits have similar function, but not guaranteed
    '''
    downstream=RelationshipTo("HitRegion", "hasDownstream", cardinality=ZeroOrOne, model=hasDownstream)
    upstream=RelationshipFrom("HitRegion", "hasDownstream", cardinality=ZeroOrOne, model=hasDownstream)
    #temperary
    affiliates=RelationshipTo('Hit','hasAffiliate',cardinality=OneOrMore,model=hasAffiliate)


class Hit(Segment):
    '''
    original annotation entries.
    Filtered by hhblits' probab, but no further procession.
    could be overlapped with each other.
    '''
    strand = StringProperty(required=True,choices=choices(['s','a']))
    frame = StringProperty(required=True,choices=choices(['0','1','2']))
    
    hit_begin = IntegerProperty(required=True)
    hit_end = IntegerProperty(required=True)
    
    aligned_seq =  StringProperty(required=True,help_text='how Hit is aligned to its HitFamily')
    aligned_consensus = StringProperty(required=True,help_text='digit-wise hit score for Hit-HitFamily alignment')
    
    analog = RelationshipTo("Hit", "hasAnalog", cardinality=ZeroOrMore, model=hasAnalog)
    analog_of = RelationshipFrom("Hit", "hasAnalog", cardinality=ZeroOrMore, model=hasAnalog)

    family = RelationshipFrom("HitFamily", "hasAlignment", cardinality=ZeroOrMore, model=hasMember)


class Family(StructuredNode):
    '''
    abstract, not shown in neo4j graph
    a group of Segments
    '''
    __abstract_node__ = True
    name = StringProperty(required=True)


class HitFamily(Family):
    '''
    
    '''
    source = StringProperty(required=True)
    accession = StringProperty(unique_index=True)
    subtype = StringProperty(help_text=" hmm model's pfam_type (family, domain, motif,...)")
    annotation = StringProperty()
    std_length = IntegerProperty(help_text=" hmm model's length1")
    members =RelationshipTo("Hit", "hasMember", cardinality=OneOrMore, model=hasMember)
    
