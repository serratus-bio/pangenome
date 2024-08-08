
from __future__ import annotations
from neomodel import (StructuredNode, StructuredRel,
                    OneOrMore, ZeroOrOne,ZeroOrMore,One,
                    RelationshipTo,RelationshipFrom,
                    StringProperty, FloatProperty,
                    IntegerProperty,UniqueIdProperty)
from neomodel import config, db
from neomodel import clear_neo4j_database
import pandas as pd
from typing import Optional
# clear_neo4j_database(db,clear_constraints=True,clear_indexes=True)

# %%
    
class HasregRel(StructuredRel):
    '''
    
    '''
    regid = FloatProperty(required=True)
    iselemental = IntegerProperty(required=True)

class HasmemberRel(StructuredRel):
    evalue=FloatProperty()
    score=FloatProperty()

class HasantiRel(StructuredRel):
    pass

class HasdownstreamRel(StructuredRel):
    begin_at = IntegerProperty(required=True)
    nested = IntegerProperty(required=True)

class HasdownstreamsetRel(StructuredRel):
    degree = IntegerProperty(required=True)
        
class HasnextdomainsetRel(StructuredRel):
    degree = IntegerProperty(required=True)
    # @property
    # def degree(self):
    #     ''' might be too slow'''
    #     params={'s':self.start_node().name,'e':self.end_node().name}
    #     q=(' MATCH (n:FuncDomainSet {name:$s})-[:HasMember]->(x:FuncDomain)-[:HasDownstream]->()'
    #         '-[:HasDownstream]->(y:FuncDomain)<-[:HasMember]-(m:FuncDomainSet {name:$e}) \n'
    #         'RETURN count(x) as degree')
    #     # to_dataframe(db.cypher_query(q,resolve_objects=True))
    #     # to_dataframe(db.cypher_query(q,resolve_objects=True))
    #     return db.cypher_query(q, params, resolve_objects=False)[0][0][0]
      
class HasnextdomainRel(StructuredRel):
    begin_at = IntegerProperty(required=True)
    nested = IntegerProperty(required=True)

class HomologousToRel(StructuredRel):
    identity = FloatProperty(required=True)
    
class HasclanentryRel(StructuredRel):
    linkscore = FloatProperty(required=True) 
    # proportional to the number of sequences in the full alignment. 
    # see 2016 paper

class HasclanhomologyRel(StructuredRel):
    linkscore = FloatProperty(required=True) 
    # the evalue of hmm-hmm alignment
    # 
  
class Fasta(StructuredNode):
    '''
    Original Record of a full nt sequence 
    '''
    # __abstract_node__ = True
    # id=UniqueIdProperty()
    name = StringProperty(required=True)
    seq = StringProperty(required=True)
    
    hasreg = RelationshipTo("Region", "HasReg", cardinality=OneOrMore, model=HasregRel)
    hasanti = RelationshipTo("Fasta", "HasAnti", cardinality=ZeroOrOne, model=HasantiRel)
    isanti = RelationshipFrom("Fasta", "HasAnti", cardinality=ZeroOrOne, model=HasantiRel)
        
class Genome(Fasta):
    '''
    a complete genome sequence
    '''
    source = StringProperty(required=True)
    accession = StringProperty(unique_index=True)
    annotation = StringProperty()
    taxonomy = StringProperty()
    #circular


class Region(StructuredNode):
    '''
    a region in genome
    '''
    # __abstract_node__ = True
    # id=UniqueIdProperty()
    name = StringProperty(required=True) #considered to be unique
    b = IntegerProperty(required=True)
    e = IntegerProperty(required=True)
    
    regfrom = RelationshipFrom("Fasta", "HasReg", cardinality=One, model=HasregRel)
    memberof = RelationshipFrom("RegionSet", "HasMember", cardinality=One, model=HasmemberRel)
    
    hasdownstream = RelationshipTo("Region", "HasDownstream", cardinality=One, model=HasdownstreamRel)
    downstreamof = RelationshipFrom("Region", "HasDownstream", cardinality=One, model=HasdownstreamRel)
    
    @property
    def seq(self):
        mother_seq:str=self.regfrom[0].seq
        return mother_seq[self.b:self.e]

    def connect_regionset(self,regionprops:dict,rel_props:dict):
        raise NotImplementedError
    
    def connect_fasta(self,fastaprops:dict,rel_props:dict):
       genome=Genome.get_or_create(fastaprops)[0]
       self.regfrom.connect(genome,rel_props)
       
    def connect_fasta_shortcut(self,genome:Genome,regid:float,iselemental:int):
        self.regfrom.connect(genome,
        properties={"regid":regid,"iselemental":iselemental})
    
    def connect_last_region(self,last_region:Optional[Region],
            begin_at:int=0,nested:int=0):
        if isinstance(last_region,Region):
            # if not self.downstreamof.relationship(last_region):
            #     self.downstreamof.disconnect(last_region)
            self.downstreamof.connect(last_region,
                properties={"begin_at":begin_at,"nested":nested})
            
            last_regionset:RegionSet=last_region.memberof[0]
            this_regionset:RegionSet=self.memberof[0]
            rel:HasdownstreamsetRel=last_regionset.hasdownstreamset.relationship(this_regionset)
            if rel is None:
                last_regionset.hasdownstreamset.connect(this_regionset,
                properties={"degree":1})
            else:
                rel.degree=rel.degree+1
                rel.save()
        else:
            if last_region is not None:
                raise TypeError(f'last_region has wrong type: {type(last_region)}')
       
       
class FuncDomain(Region):
    memberof = RelationshipFrom("FuncDomainSet", "HasMember", cardinality=One, model=HasmemberRel)
    hasdownstream = RelationshipTo("DomainLinkage", "HasDownstream", cardinality=One, model=HasdownstreamRel)
    downstreamof = RelationshipFrom("DomainLinkage", "HasDownstream", cardinality=One, model=HasdownstreamRel)
    
    nextdomainset=RelationshipTo("FuncDomain", "hasNextDomain", cardinality=ZeroOrOne, model=HasnextdomainRel)
    lastdomainset=RelationshipFrom("FuncDomain", "hasNextDomain", cardinality=ZeroOrOne, model=HasnextdomainRel)
    
    homologous_to=RelationshipTo("FuncDomain", "homologousTo", cardinality=ZeroOrMore, model=HomologousToRel)
    homologous_ref_of=RelationshipFrom("FuncDomain", "homologousTo", cardinality=ZeroOrMore, model=HomologousToRel)
    
    def connect_regionset(self,regionprops:dict,rel_props:dict={}):
        domainset=FuncDomainSet.get_or_create(regionprops)[0]
        self.memberof.connect(domainset)

    
class HmmFuncDomain(FuncDomain):
    hmmstart = IntegerProperty(required=True)
    hmmend=IntegerProperty(required=True)
    #evalue!

class DomainLinkage(Region):
    memberof = RelationshipFrom("DomainLinkageSet", "HasMember", cardinality=One, model=HasmemberRel)
    hasdownstream = RelationshipTo("FuncDomain", "HasDownstream", cardinality=One, model=HasdownstreamRel)
    downstreamof = RelationshipFrom("FuncDomain", "HasDownstream", cardinality=One, model=HasdownstreamRel)
    def connect_regionset(self,regionprops:dict,rel_props:dict={}):
        domainset=DomainLinkageSet.get_or_create(regionprops)[0]
        self.memberof.connect(domainset)
    
    
class RegionSet(StructuredNode):
    # __abstract_node__ = True
    # id=UniqueIdProperty()
    name = StringProperty(required=True)
    
    hasdownstreamset=RelationshipTo("RegionSet", "hasDownstreamSet", cardinality=ZeroOrMore, model=HasdownstreamsetRel)
    downstreamsetof=RelationshipFrom("RegionSet", "hasDownstreamSet", cardinality=ZeroOrMore, model=HasdownstreamsetRel)
    hasmember=RelationshipTo("Region", "HasMember", cardinality=OneOrMore, model=HasmemberRel)
    
    
class FuncDomainSet(RegionSet):
    source = StringProperty(required=True)
    accession = StringProperty(unique_index=True)
    annotation = StringProperty()
    std_length = IntegerProperty()
    pfam_subtype = StringProperty()
    hasdownstreamset=RelationshipTo("DomainLinkageSet", "hasDownstreamSet", cardinality=ZeroOrMore, model=HasdownstreamsetRel)
    downstreamsetof=RelationshipFrom("DomainLinkageSet", "hasDownstreamSet", cardinality=ZeroOrMore, model=HasdownstreamsetRel)
    
    nextdomainset=RelationshipTo("FuncDomainSet", "hasNextDomainSet", cardinality=ZeroOrMore, model=HasnextdomainsetRel)
    lastdomainset=RelationshipFrom("FuncDomainSet", "hasNextDomainSet", cardinality=ZeroOrMore, model=HasnextdomainsetRel)
    
    clanentryof = RelationshipTo("FuncDomainClan", "hasClanEntry", cardinality=ZeroOrOne, model=HasclanentryRel) # clan ->self
    hasclanhomology = RelationshipTo("FuncDomainSet", "hasClanHomology", cardinality=ZeroOrMore, model=HasclanhomologyRel) # self -> peer
    clanhomologyof = RelationshipTo("FuncDomainSet", "hasClanHomology", cardinality=ZeroOrMore, model=HasclanhomologyRel) # peer -> self
    
    
class DomainLinkageSet(RegionSet):
    hasdownstreamset=RelationshipTo("FuncDomainSet", "hasDownstreamSet", cardinality=ZeroOrOne, model=HasdownstreamsetRel)
    downstreamsetof=RelationshipFrom("FuncDomainSet", "hasDownstreamSet", cardinality=ZeroOrOne, model=HasdownstreamsetRel)

    
class RegionSetClan(StructuredNode):
    name = StringProperty(required=True)
    hasclanentry = RelationshipTo("FuncDomainSet", "hasClanEntry", cardinality=ZeroOrMore, model=HasclanentryRel)


class FuncDomainClan(RegionSetClan):
    source = StringProperty(required=True)
    accession = StringProperty(unique_index=True)
    annotation = StringProperty()
    
    