#!/usr/bin/env python

"""
amnonscript
dbstruct.py

the structure for bactdb (to get info from the bacterial database (SRBactDB.py))
separate file so won't be erased when reloading
"""

__version__ = "0.1"

class dbstruct:
	Ontology={}
	OntologyNames={}
	OntologyIDs={}
	# hold the built graph structures for fast drawing (from InitOntoGraph())
	OntoGraph={}
	dbfile=''
	recPrecBalance=10

