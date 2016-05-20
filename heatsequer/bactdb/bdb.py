#!/usr/bin/env python

"""
amnonscript
heatsequer
bactdb.py

get info from the bacterial database
"""

__version__ = "1.0"

from ..utils.amnonutils import Debug,dictupper,listupper,delete
from ..utils.graphutils import simplifygraph

import numpy as np
import matplotlib.pyplot as plt
import biom
from matplotlib.pyplot import *
from networkx.drawing.nx_agraph import graphviz_layout
import csv
#from scipy import cluster
#from scipy import spatial
#from sklearn.preprocessing import scale
#import copy
import sqlite3
from collections import defaultdict
# for debugging - use XXX()
from pdb import set_trace as XXX
from scipy import stats
import sys


#def initdb(dbname="/Users/amnon/Python/SRBactDB/SRBactDB.db"):
def initdb(dbname="db/SRBactDB.db"):
	db=dbstart(dbname)
	db=initontologies(db)
	return db


def getontonames(db):
	"""
	returns a list of all the ontology fields to load and their ontology types
	input:
	db : dbstruct
	output:
	ontonames : list of strings
		the names of all fields to load for ontology plots
	ontotypes : list of strings
		the ontology name for each field
	"""
	ontonames=["HOST_TAXID","ENV_MATTER","ENV_FEATURE","ENV_BIOM","COUNTRY","BODY_SITE"]
	ontotypes=["NCBITAX",   "ENVO",      "ENVO",       "ENVO",    "GAZ",    "UBERON"]
	return ontonames,ontotypes


def initontologies(db):
	"""
	load and preprocess the ontologies for the database
	input:
	db : dbstruct
		the database (from dbstart())
	output:
	db : dbstruct
		after loading and initializing ontologies
	"""

	db=LoadOntologies(db)
	ontofields,ontonames = getontonames(db)
	for contofield,contotype in zip(ontofields,ontonames):
		db=InitOntologyGraph(db,[0],contofield,contotype)

	# db=InitOntologyGraph(db,[0],'HOST_TAXID','NCBITAX')
	# db=InitOntologyGraph(db,[0],'ENV_MATTER','ENVO')
	# db=InitOntologyGraph(db,[0],'ENV_BIOM','ENVO')
	# db=InitOntologyGraph(db,[0],'ENV_FEATURE','ENVO')
	# db=InitOntologyGraph(db,[0],'COUNTRY','GAZ')
	# db=InitOntologyGraph(db,[0],'BODY_SITE','UBERON')
	db.ontologyinit=True
	return db


#def dbconnect(db,dbname="/Users/amnon/Python/SRBactDB/SRBactDB.db"):
def dbconnect(db,dbname="db/SRBactDB.db"):
	"""
	connect to the database
	input:
	db - from dbstart()
	dbname - the name of the database to connect to
	"""

	db.dbfile=dbname
	# the database connection
	Debug(1,"Connecting to database ",dbname)
	db.con=sqlite3.connect(db.dbfile)
	Debug(1,"Connected")
	# and the cursor
	db.cur=db.con.cursor()
	# test if the database file exists:
	try:
		db.cur.execute("SELECT count(*) FROM sqlite_master WHERE type='table' AND name='Reads'")
		istable=db.cur.fetchone()
		if istable[0]==0:
			Debug(9,"Can't find Reads table in database file %s" % dbname)
			db=False
	except:
		Debug(9,"Can't test database file %s" % dbname)
		db=False

	return db


class dbstruct:
	def __init__(self):
		# True if ontology is initialized, False if not ready yet (need to call InitOntologies)
		self.ontologyinit=False

		self.Ontology={}
		self.OntologyNames={}
		self.OntologyIDs={}
		# hold the built graph structures for fast drawing (from InitOntoGraph())
		self.OntoGraph={}
		self.dbfile=''
		self.recPrecBalance=10


#def dbstart(dbname="/Users/amnon/Python/SRBactDB/SRBactDB.db"):
def dbstart(dbname="db/SRBactDB.db"):
	'''
	start the database structure and connect to database
	input:
	dbname : string
		the name of the database to connect to
	output:
	db : dbstruct
		the database variable
	'''

	db=dbstruct()

#	db=dbstruct.dbstruct()
	db=dbconnect(db,dbname)

	# the minimal frequency in order for a sequence to be present
	db.MINFREQ=0.0

	# minimal number of reads in sample order to include in database
	db.MINEXPREADS=500

	# the minimal number of saples in order to show a leaf
	db.MINLEAFSIZE=0
	db.ontograph={}
	return db


def LoadOntologyOBO(db,filename,ontologyname):
	"""
	load an ontology and store it in an ontology hash
	each ontology is a hash of ids containing lists of ids of parents
	and an additional hash of name to ids (including synonyms????)
	input:
	filename - of the .obo ontology
	ontologyname - for the hash of ontologies
	"""
	ofile=open(filename,"rU")
	ohash={}
	namehash={}
	idhash={}
	numids=0
	numisa=0
	numconsider=0
	numsyn=0
	numlocin=0
	for cline in ofile:
		cline=cline.strip()
		cline=cline.lower()
		terms=cline.split(': ')
		if terms[0]=='id':
			cid=terms[1]
			numids+=1
		if terms[0]=='name':
			cname=terms[1]
			cname=cname
			namehash[cname]=cid
			idhash[cid]=cname
		if terms[0]=='is_a':
			idandname=terms[1].split(' ! ')
			# to remove the " {XXX}" which sometimes appears
			simid=idandname[0].split(" ")[0]
#                simname=idandname[1].lower()
			if cid not in ohash:
				ohash[cid]=[]
			ohash[cid].append(simid)
			numisa+=1
		if terms[0]=='consider':
			simid=terms[1]
			if cid not in ohash:
				ohash[cid]=[]
			ohash[cid].append(simid)
			numconsider+=1
		if terms[0]=='synonym':
			nameandrem=terms[1]
			simdat=nameandrem.split('"')
			simname=simdat[1]
			namehash[simname]=cid
			numsyn+=1
		if terms[0]=='relationship':
			# get the type of reationship
			cont=terms[1].split(" ",1)
			if cont[0]=="part_of":
				idandname=cont[1].split(' ! ')
				# to remove the " {XXX}" which sometimes appears
				simid=idandname[0].split(" ")[0]
#                    simname=idandname[1].lower()
				if cid not in ohash:
					ohash[cid]=[]
					ohash[cid].append(simid)
					numisa+=1
			if cont[0]=="located_in":
				otherpart=cont[1].split(" ")
				simid=otherpart[0]
				if cid not in ohash:
					ohash[cid]=[]
					ohash[cid].append(simid)
					numlocin+=1
	Debug(0,"loaded ontology",ontologyname,"from file",filename)
	Debug(0,"ids",numids,"isa",numisa,"consider",numconsider,"synonym",numsyn,"located_in",numlocin)
	db.Ontology[ontologyname]=ohash
	db.OntologyNames[ontologyname]=namehash
	db.OntologyIDs[ontologyname]=idhash
	# add the term for sampls that dont contain the field
	db.Ontology[ontologyname]['not found']=[]
	db.OntologyNames[ontologyname]['not found']='not found'
	db.OntologyIDs[ontologyname]['not found']='not found'
	return db


def LoadOntologies(db):
	"""
	load the ontologies
	"""
	db=LoadOntologyOBO(db,"/Users/amnon/Databases/ontologies/envo.obo","ENVO")
	db=LoadOntologyOBO(db,"/Users/amnon/Databases/ontologies/gaz.obo","GAZ")
	db=LoadOntologyOBO(db,"/Users/amnon/Databases/ontologies/ncbitax.obo","NCBITAX")
	db=LoadOntologyOBO(db,"/Users/amnon/Databases/ontologies/uberon.obo","UBERON")
	Debug(1,"Loaded ontologies")
	return db


def InitOntologyGraph(db,pv,field,ontohashname):
	"""
	just init the ontology graph for the field
	keep it in a dictionary of field names in self.OntoGraph[fieldname]
	including the graph layout for fast drawing
	"""

	import networkx as nx

	if not ontohashname in db.Ontology:
		Debug(10,"Ontology not loaded",ontohashname)
		return()

	G=nx.DiGraph()
	nodes=[]
	nodelabs=[]
	edges=[]
	sizes=[]

	Debug(0,"Ontology test field:",field)
	db.cur.execute("SELECT Value,SampleID FROM Maps WHERE Field=?",[field])
	allvals=db.cur.fetchall()
	chash={}
	studies={}
	sids=[]
	nodestudies={}

	# add all samples to the ontology tree
	for cv in allvals:
		sids.append(cv[1])
		ontoname=cv[0].lower()
		# if we're in taxonomy, need to get name instead of id
		if ontohashname=="NCBITAX":
			ontoid="ncbitaxon:"+ontoname
			if ontoid in db.OntologyIDs[ontohashname]:
				ontoname=db.OntologyIDs[ontohashname][ontoid]
			else:
				Debug(2,"NCBITAX ontoname for ontoid",ontoname,"not found!!! in sample",cv[1])
				ontoname="not found"
				ontoid="not found"
		# otherwise just remove the header "GAZ:" etc.
		else:
			ontoname=ontoname[len(ontohashname)+1:]
			if ontoname in db.OntologyNames[ontohashname]:
				ontoid=db.OntologyNames[ontohashname][ontoname]
			else:
				Debug(2,"ontology",ontohashname," id for ontoname",ontoname,"not found!!!")
				ontoid="not found"
				ontoname="not found"

		if ontoid!="":
			chash=AddOntologySample(db.Ontology[ontohashname],chash,cv[1]-1,ontoid)
			# get the study for the sample and add to dict
			db.cur.execute("SELECT StudyID FROM Samples WHERE SampleID=?",[cv[1]])
			studyid=db.cur.fetchone()
			studies[int(cv[1]-1)]=studyid

	# add all unaccounted samples (i.e. don't have this field)
	for cpos in range(len(pv)):
		if cpos+1 not in sids:
			Debug(1,"in sample",cpos+1,"field was not found")
			chash=AddOntologySample(db.Ontology[ontohashname],chash,cpos,"not found")
			# get the study for the sample and add to dict
			db.cur.execute("SELECT StudyID FROM Samples WHERE SampleID=?",[cpos+1])
			studyid=db.cur.fetchone()
			studies[cpos]=studyid

	# and add nodes to graph
	# and create study list of each node
	for (conto,csamp) in chash.items():
		if conto in db.OntologyIDs[ontohashname]:
			contoname=db.OntologyIDs[ontohashname][conto]
		else:
			contoname="not found"
		Debug(1,"Value",conto,"which is",contoname)
		Debug(1,"present in",len(csamp),"samples")

		# add node and parents to the graph
		if contoname not in nodes:
			# init the number of studies per node
			alls=[]
			for cs in csamp:
				if cs in studies:
					alls.extend(studies[cs])
				else:
					Debug(1,"Sample",cs,"Does not have a study")
			nodestudies[contoname]=alls

			nname=contoname+"-"+str(len(csamp))
			if len(csamp)>db.MINLEAFSIZE:
				nodes.append(contoname)
				nodelabs.append(nname)
				sizes.append(len(csamp))

			if conto in db.Ontology[ontohashname]:
				for cpar in db.Ontology[ontohashname][conto]:
					if cpar in db.OntologyIDs[ontohashname]:
						parname=db.OntologyIDs[ontohashname][cpar]
					else:
						parname="not found"
					edges.append((contoname,parname))
			else:
				Debug(5,"didnt find onto",conto)
				if contoname in nodes:
					sizes[nodes.index(contoname)]+=len(csamp)
				else:
					debug(5,"and not in nodes")

	Debug(6,"Created node data")
	Debug(0,"Got nodes:",len(nodes),"sizes:",len(sizes))

	# create the node data to the graph
	fraccol=np.zeros(len(nodes))
	ndict={}
	for a in range(len(nodes)):
		Debug(0,a,nodes[a],sizes[a])
		G.add_node(a)
		ndict[a]=nodelabs[a]

	# Create the edge data for the graph
	edges2=[]
	for pairs in edges:
		if pairs[0] in nodes and pairs[1] in nodes:
			pos1=nodes.index(pairs[0])
			pos2=nodes.index(pairs[1])
			edges2.append((pos1,pos2))
	G.add_edges_from(edges2)
	Debug(6,"Created edge data")

	# add the node for each experiment and connect to leaves
	lastpos=len(nodes)
	oldlen=len(nodes)
	for a in range(oldlen):
		cnode=nodes[a]
		neighborset=[]
		for cneighbor in G.predecessors_iter(a):
			neighborset.extend(nodestudies[nodes[cneighbor]])
		newstudies=list(set(nodestudies[cnode]).difference(set(neighborset)))
		Debug(0,"node",cnode,"Has new studies",newstudies)

		for cstudy in newstudies:
			# find all the samples that contribute to the study
			ccpos=[ti for ti, tx in enumerate(nodestudies[cnode]) if tx == cstudy]
			Debug(0,"For study",cstudy,"we have",ccpos,"samples")
			# store the samples in this study node which contribute to it
			samplelist=[]
			for (ti,tx) in enumerate(nodestudies[cnode]):
				if tx==cstudy:
					try:
						samplelist.append(chash[db.OntologyNames[ontohashname][cnode]][ti])
					except:
						Debug(6,cnode,'Not found in ontology! - did not add samples!')

			# check each sample of the study if the sequence is present
			studytot=len(ccpos)
			# add all new experiment nodes to the node processed
			ndict[lastpos]="study-"+str(cstudy)
			fraccol=np.append(fraccol,0)
			sizes.append(studytot)
			nodes.append(lastpos)
			chash[lastpos]=samplelist
			G.add_node(lastpos)
			G.add_edge(lastpos,a)
			lastpos+=1

	# # remove non informative nodes (1 in 1 out)
	# delpos=simplifygraph(G)
	# nodes=delete(nodes,delpos)
	# sizes=delete(sizes,delpos)
	# for cdel in delpos:
	# 	del ndict[cdel]

	# and save the taxonomy tree data
	db.OntoGraph[field]={}
	db.OntoGraph[field]['G']=G
	db.OntoGraph[field]['chash']=chash
	db.OntoGraph[field]['ontohashname']=ontohashname
	db.OntoGraph[field]['nodes']=nodes
	db.OntoGraph[field]['sizes']=sizes
	db.OntoGraph[field]['ndict']=ndict
#	db.OntoGraph[field]['nodepositions']=nx.graphviz_layout(G)
	if sys.version_info[0]<3:
		db.OntoGraph[field]['nodepositions']=graphviz_layout(G)
	else:
		db.OntoGraph[field]['nodepositions']=[]
		Debug(9,"Can't init the graph layout since using python 3")
	return db


def PlotOntologyGraph(db,seq,field,tofig=False,toax=False):
	"""
	Plot the ontology graph for sequence seq using ontology field field
	must init the graphs first using InitOntologyGraph()
	tofig - the figure to plot to or False to plot new fig
	"""

	import networkx as nx

	if not db.ontologyinit:
		Debug(8,"Ontologies not initialized - initializing now (may take a minute)")
		db=initontologies(db)
		Debug(8,"Ontologies initialized")

	pv=GetSeqVec(db,seq)
	if len(pv)==0:
		Debug(8,'No reads for bacteria in BactDB!')
		return
	nz=np.where(pv>db.MINFREQ)
	totalappear=len(nz[0])
	Debug(7,"--------Seq",seq)
	Debug(7,"Found in",totalappear,"Samples")
	Debug(7,"Fraction=",float(totalappear)/len(pv))
	Debug(1,nz[0])

	chash=db.OntoGraph[field]['chash']
	cnodepositions=db.OntoGraph[field]['nodepositions']
	ontohashname=db.OntoGraph[field]['ontohashname']
	nodes=db.OntoGraph[field]['nodes']
	sizes=db.OntoGraph[field]['sizes']
	ndict=db.OntoGraph[field]['ndict']
	G=db.OntoGraph[field]['G']

	# check if we don't care for the useme field
	useme=np.ones(len(pv))

	db.selection_max_enrich=0
	db.selection_field="none"
	db.selection_value="none"
	db.selection_samples=[]
	db.selection_fraction=[]

	fracs={}
	for (conto,csamp) in chash.items():
		if conto in db.OntologyIDs[ontohashname]:
			contoname=db.OntologyIDs[ontohashname][conto]
		else:
			try:
				contoname=int(conto)
				contoname=conto
				Debug(1,'study item',conto)
			except:
				contoname="not found"
				Debug(1,'not found',conto)
		Debug(1,"Value",conto,"which is",contoname)
		Debug(1,"present in",len(csamp),"samples")
		# get the fraction of samples which have the sequence for this node
		(cpval,cfracenrich,cnzincval,cfraccval)=GetStatsUse(db,pv,csamp,useme,field,contoname)
		fracs[contoname]=cfraccval

	# # create the node data to the graph
	fraccol=np.zeros(len(nodes))
	for a in range(len(nodes)):
		Debug(0,a,nodes[a],sizes[a])
		if nodes[a] in fracs:
			fraccol[a]=fracs[nodes[a]]
			Debug(1,"color for",nodes[a],"is",fraccol[a])
		else:
			Debug(6,"fraction not found for",nodes[a])

	# save the graph
	tg=G.copy()
	sizedict={}
	colordict={}
	labsizedict={}
	labshapedict={}
	for cnode in tg.nodes():
		sizedict[cnode]=sizes[cnode]
		colordict[cnode]=int(1000*fraccol[cnode]/(max(fraccol)+0.00001))
		if ndict[cnode][:5]=='study':
			labsizedict[cnode]=1
			labshapedict[cnode]=1
		else:
			labsizedict[cnode]=int(1000*fraccol[cnode]/(max(fraccol)+0.00001))
			labshapedict[cnode]=0
	nx.set_node_attributes(tg,'nodenames',ndict)
	nx.set_node_attributes(tg,'nodecolors',sizedict)
	nx.set_node_attributes(tg,'nodesizes',colordict)
	nx.set_node_attributes(tg,'nodelabelsizes',labsizedict)
	nx.set_node_attributes(tg,'nodeshapes',labshapedict)
	simplifygraph(tg,'nodeshapes')
	tg=tg.to_undirected()
#	nx.set_node_attributes(tg,'origsizes',sizes)
	nx.write_graphml(tg,'/Users/amnon/test.graphml')
	Debug(6,'saved graphml ~/test.graphml')

	# and draw the graph!!!
	if toax:
		nx.draw_networkx(G,nodes=range(len(nodes)),node_size=1000*fraccol/(max(fraccol)+0.00001),labels=ndict,node_color=sizes,vmin=-1,vmax=50,with_labels=True,font_color='g',pos=cnodepositions,ax=toax)
		return

	if tofig:
		figure(tofig.number)
		clf()
	else:
		plt.figure()


	nx.draw_networkx(G,nodes=range(len(nodes)),node_size=1000*fraccol/(max(fraccol)+0.00001),labels=ndict,node_color=sizes,vmin=-1,vmax=50,with_labels=True,font_color='g',pos=cnodepositions)



def GetStatsUse(db,pv,csamp,usepos,field,value):
	"""
	# Get the statistics for a given division using part of the samples
	# input:
	# pv - the total reads vector (from GetGGIDVec)
	# csamp - an array of sample positions whereto test enrichment
	# usepos - the positions to use - array of size pv (positions with 0 are ignored in the calculation)
	# field,value - the field and value name to be stored
	# output:
	# pval - the p-value to egt such enrichment (binomial)
	# franenrich - the fold enrichment compared to other samples
	# nzincval - number of elements in the catergory
	# fraccval - fraction of present in the category
	"""
	uupos=np.where(usepos>0)[0]
	nz=np.where(pv[uupos]>db.MINFREQ)
	totalappear=len(nz[0])
	Debug(0,"Total appear is ",totalappear)
	totalsamples=len(np.where(usepos>0)[0])
	Debug(0,"Total samples",totalsamples)
	if totalsamples==0:
		return(1,0,0,0)
	fractotal=np.float32(totalappear)/totalsamples
	Debug(0,"Fraction total samples",fractotal)

	ucsamp=np.zeros(len(usepos))
	ucsamp[csamp]=1
	ucsamp=ucsamp*usepos
	ucsamp=np.where(ucsamp>0)[0]
	Debug(0,"Length of ucsamp",len(ucsamp))
	Debug(0,ucsamp)
	tmp=np.where(pv[ucsamp]>db.MINFREQ)
	nzincval=len(tmp[0])
	Debug(1,"*Present in",nzincval,"Of the focus samples")
	totincval=len(ucsamp)
	if totincval==0:
		Debug(3,"No samples in the intersect",field,value)
		return(1,0,0,0)
	fraccval=np.float32(nzincval)/totincval
	pval=1-stats.binom.cdf(nzincval-1,totincval,fractotal)
	Debug(1,"present in other",totalappear-nzincval,"samples in other",totalsamples-totincval)

	# enrichemtn metric
	fracenrich=fraccval*(totalsamples-totincval)/(db.recPrecBalance+np.float32(totalappear-nzincval))

	# F-metric
	recall=float(nzincval)/(totalappear+0.000001)
	precision=float(nzincval)/(totincval)
	fracenrich=2*recall*precision/(recall+precision+0.0000001)

	# dani distance metric
#        fracenrich=self.DistDani(nzincval,totincval,totalappear,totalsamples)

	Debug(1,"nzincval",nzincval,"totincval",totincval,"fraccval",fraccval,"pval",pval,"enrichment",fracenrich)
	if pval<0.001 and nzincval>2 and fracenrich>0.3:
		Debug(5,field,value,"Fraction",fraccval,"Enrichment",fracenrich,"present in",nzincval)
	# store the new values if it is the best performance yet
	if pval<0.1 and nzincval>2 and fracenrich>db.selection_max_enrich:
		Debug(0,"selected",field,value,"fracenrich=",fracenrich,"old val was",db.selection_max_enrich)
		db.selection_max_enrich=fracenrich
		db.selection_samples=csamp
		db.selection_pval=pval
		db.selection_fraction=fraccval
		db.selection_field=field
		db.selection_value=value
	return(pval,fracenrich,nzincval,fraccval)


def AddOntologySample(ohash,chash,sampleid,ontid):
	"""
	add sample sampleid to ontid and all it's parents
	input:
	ohash - the ontology hash (from LoadOntologyOBO)
	chash - the current hash structure to fill with sampleids
	sampleid - the sampleid to add
	ontid - the ontologyid (and all its parents) to fill with the sample id
	output:
	chash - the new hash
	"""

	if not ontid in chash:
		chash[ontid]=[]
	chash[ontid].append(sampleid)
	if not ontid in ohash:
		Debug(0,"parent not found for",ontid)
		return(chash)
	# add sample to hash
	Debug(0,"found parent for",ontid)
	parlist=list(ohash[ontid])
	for parents in parlist:
		chash=AddOntologySample(ohash,chash,sampleid,parents)
	return(chash)


def GetSeqID(db,osequence,insert=False):
	'''
	Get sequence id from database
	input:
	db - from dbconnect()
	osequence - the sequence to look for
	insert - True to insert sequence to database if not found, 0 to return 0 instead
	'''

	sequence=osequence.upper()
	db.cur.execute("SELECT SeqID FROM Sequences WHERE Sequence = ?",[sequence])
	res=db.cur.fetchone()
	if res:
		seqid=res[0]
		Debug(0,"Found sequence",sequence,"SeqID=",seqid)
	else:
		if insert:
			Debug(0,"Sequence",sequence,"Not found - creating in table")
			db.cur.execute("INSERT INTO Sequences (Sequence) VALUES (?)",[sequence])
			seqid=db.cur.lastrowid
			Debug(0,"SequenceID is",seqid)
		else:
			Debug(3,"Sequence not found in database",sequence)
			seqid=0
	return(seqid)


def GetSeqIDVec(db,seqid):
	'''
	Get the read vector from the database for a given sequence id
	input:
	seqid - the sequence id to look for in the reads database
	output:
	pv - the row of reads for the seqid in the database (column is defined by SampleID-1)
	'''

	# Get the number of samples in the database (for comaprison array)
	db.cur.execute("SELECT MAX(SampleID) FROM Samples")
	res=db.cur.fetchone()
	numsamples=res[0]
	Debug(0,"Maximal SampleID is",numsamples)
	# take all reads for the seqID
	Debug(0,"Getting reads")
	pv=np.zeros(numsamples)
	db.cur.execute("SELECT SampleID,Reads FROM Reads WHERE SeqID=?",[seqid])
	Debug(0,"Putting in vector")
	found=0
	for crv in db.cur:
		pv[crv[0]-1]+=crv[1]
		found+=1
	if found==0:
		Debug(5,"SeqID",seqid,"Has no reads in database")
	Debug(2,"Found in",found,"Read entries")
	return(pv)


def GetSeqVec(db,seq):
	try:
		seqid=int(seq)
	except:
		if len(seq)!=89:
			seq=seq[0:89]
			Debug(1,"reducing to 89bp")
		Debug(0,"Getting vector for sequence",seq)
		db.cur.execute("SELECT SeqID from Sequences WHERE Sequence=?",[seq])
		res=db.cur.fetchone()
		if not(res):
			Debug(3,"Sequence",seq,"Not found!")
			return([])
		seqid=res[0]
	Debug(0,"found sequence, seqid=",seqid)
	pv=GetSeqIDVec(db,seqid)
	return(pv)


def StudyNameFromID(db,studyid):
	db.cur.execute('SELECT ExpName FROM Experiments WHERE StudyID=?',[studyid])
	res=db.cur.fetchall()
	nres=[]
	for a in res:
		nres.append(a[0])
	Debug(3,nres)
	return(nres[0])


def GetSeqInfo(db,seq):
	"""
	Get info from the database about the sequence seq
	input:
	db - from dbconnect()
	seq - the sequence to search for (min 89 bp)
	output:
	totalappear - total number of samples in which it appears
	numstudies - number of different studies where it appears
	allstudies - ids of all the studies where it appears
	studysamples - a dict (by studyid) of samples in the studies
	totdbsamples - the total number of samples in the database
	"""

	# Get the read vector for the samples
	pv=GetSeqVec(db,seq)
	if len(pv)==0:
		return(0,"not in db",0,"not in db",0)
	nz=np.where(pv>db.MINFREQ)
	totalappear=len(nz[0])
	Debug(2,"--------Seq",seq)
	Debug(3,"Found in",totalappear,"Samples")
	Debug(3,"Fraction=",float(totalappear)/len(pv))
	totdbsamples=len(pv)
	if totalappear==1:
		Debug(5,'Only in one sample',nz[0])

	allstudies=[]
	studysamples=defaultdict(list)
	for cpos in nz[0]:
	# get the study for the sample and add to dict
		db.cur.execute("SELECT StudyID FROM Samples WHERE SampleID=?",[int(cpos)+1])
		studyid=db.cur.fetchone()
		allstudies.append(studyid)
		studysamples[studyid[0]].append(int(cpos)+1)
	for (ck,cv) in studysamples.items():
		Debug(1,'study',ck,StudyNameFromID(db,ck))
		Debug(1,'samples:',len(cv))
	numstudies=len(set(allstudies))
	return([totalappear,numstudies,list(set(allstudies)),studysamples,totdbsamples])

	# get the most frequent samples
	nzp=pv[nz]
	si=np.argsort(nzp)
	sv=nzp[si]
	siorig=nz[si]
	for cid in siorig:
		db.cur.execute("SELECT StudyID FROM Samples WHERE SampleID=?",[int(cpos)+1])
		studyid=db.cur.fetchone()


def GetSeqInfoSummary(bactdb,seq,field='HOST_TAXID',sortednodes=[]):
	import networkx as nx

	pv=GetSeqVec(bactdb,seq)
	if len(pv)==0:
		return(0,"not in db",0,"not in db",0)
	nz=np.where(pv>bactdb.MINFREQ)
	totalappear=len(nz[0])
	fractotal=float(totalappear)/len(pv)
	Debug(7,"--------Seq",seq)
	Debug(7,"Found in",totalappear,"Samples")
	Debug(7,"Fraction=",fractotal)
	Debug(0,nz[0])
	if fractotal==0:
		return(0,"not present",0,"not present",0)

	chash=bactdb.OntoGraph[field]['chash']
	ontohashname=bactdb.OntoGraph[field]['ontohashname']
	nodes=bactdb.OntoGraph[field]['nodes']
	sizes=bactdb.OntoGraph[field]['sizes']
	ndict=bactdbf.OntoGraph[field]['ndict']
	G=bactdb.OntoGraph[field]['G']

	# go over all ontology values (in the tree) starting from the leaves
	# we're guaranteed to get all children before the parent
	maxgroupsize=0
	maxgroupname="not found"
	maxfrac=0
	maxnz=0
	if len(sortednodes)==0:
		sortednodes=nx.topological_sort(G,reverse=False)
	sfreqs=zeros(len(sortednodes))
	totfreqs=zeros(len(sortednodes))
	subtreesize=zeros(len(sortednodes))
	leafnum=zeros(len(sortednodes))
	totbactfreq=0
	for cnode in sortednodes:
		children=G.predecessors(cnode)
		if len(children)>0:
			sfreqs[cnode]=mean(sfreqs[children])
			totfreqs[cnode]=sum(totfreqs[children])
			leafnum[cnode]=sum(leafnum[children])
#                if len(children)>1:
#                    subtreesize[cnode]=max(subtreesize[children])+1
#                else:
#                    subtreesize[cnode]=max(subtreesize[children])
			for cchild in children:
				if sfreqs[cchild]>0:
					subtreesize[cnode]+=subtreesize[cchild]
			if sfreqs[cnode]>=0.25:
				if subtreesize[cnode]>maxgroupsize:
					# count to see if more than 1 child has the bacteria - otherwise it doesn't count
					numnz=0
					for cchild in children:
						if sfreqs[cchild]>0:
							numnz+=1
					if numnz>1:
						maxgroupsize=subtreesize[cnode]
						maxgroupname=nodes[cnode]
						maxfrac=sfreqs[cnode]
						maxnz=numnz
		else:
			if nodes[cnode] in chash:
				csamp=chash[cnode]
				tmp=where(pv[csamp]>self.MINFREQ)
				nzcsamp=len(tmp[0])
				sfreqs[cnode]=float(nzcsamp)/len(csamp)
				# add leaf node to total count (for min subtree finding)
				totfreqs[cnode]=sfreqs[cnode]
				totbactfreq+=sfreqs[cnode]
				subtreesize[cnode]=1
				leafnum[cnode]=1
				if sfreqs[cnode]>0.25:
					if maxgroupsize<=1:
						if maxfrac<sfreqs[cnode]:
							maxgroupsize=subtreesize[cnode]
							maxgroupname=str(nodes[G.neighbors(cnode)[0]])
							maxfrac=sfreqs[cnode]
			else:
				self.Debug(7,"Did not find node in tree",nodes[cnode])
	self.Debug(3,"max group size is "+str(maxgroupsize)+" and name is "+maxgroupname)

	# now look for the smallest tree with > threshold samples containing the bacteria
	mingroupsize=len(sortednodes)
	minnode=sortednodes[-1]
	mingroupname='error. not found'
	self.Debug(3,"total freq is ",totbactfreq)
	for cnode in sortednodes:
		if totfreqs[cnode]>0.75*totbactfreq:
			if mingroupsize>leafnum[cnode]:
				minnode=cnode
				mingroupsize=leafnum[cnode]
				mingroupname=str(nodes[cnode])
	self.Debug(3,"min group num leaves is "+str(mingroupsize)+" and name is "+mingroupname)

	return(fractotal,maxgroupname,maxgroupsize,mingroupname,mingroupsize)



def GetDBSource(db,seqs):
	"""
	Get the biggest sets covering bacteria from the list seqs from the database
	input:
	db - the database
	seqs - a list of sequences to test
	output:
	dat - an array (sequeces * db samples) of frequency of sequence in sample
	"""

	# get the number of samples in the database
	for cseq in seqs:
		pv=GetSeqVec(db,cseq)
		if len(pv)>0:
			break

	# fill the frequency array for all sequences
	dat=np.zeros((len(seqs),len(pv)))
	for idx,cseq in enumerate(seqs):
		pv=GetSeqVec(db,cseq)
		if len(pv)==0:
			continue
		dat[idx,:]=pv

	return dat


def GetSampleStudy(db,samp):
	"""
	Get the study id and name for sampleID pos
	input:
	db
	samp - the sampleID

	output:
	sid - study id
	sname - studyname
	"""

	db.cur.execute("SELECT StudyID FROM Samples WHERE SampleID=?",[samp])
	studyid=db.cur.fetchone()
	sid = studyid[0]
	sname=StudyNameFromID(db,sid)

	return sid,sname

def SamplesInStudy(db,studyid):
	"""
	Get the list of samples in a study
	input:
	db
	studyid - the id of the study to examine

	output:
	samples - a list of sampleids in the study
	"""
	db.cur.execute('SELECT SampleID FROM Samples WHERE StudyID=?',[studyid])
	res=db.cur.fetchall()
	samples=[]
	for a in res:
		samples.append(a[0])
	return samples


def GetSampleMapField(db,sampleid):
	"""
	Get the value of the mapping file field for the given sample
	input:
	db
	sampleid - the sampleid to look for
	field - field name (i.e. "ENV_MATTER")

	output:
	val - the value of the field
	"""

	vals={}
	db.cur.execute('SELECT Field,Value FROM Maps WHERE SampleID=?',[sampleid])
	res=db.cur.fetchall()
	for a in res:
		vals[a[0]]=a[1]
	return(vals)


def GetSeqListInfo(db,seqs,info='samples'):
	"""
	Get information about a sequence list.
	input:
	db
	seqs - a list of sequences (ACGT)
	info - the info type to collect:
		'samples' - the samples for each sequence
		'studies' - the studies where each seq appears
		'types' - env_matter+host_taxid

	output:
	res - a dict with info type as keys, each containing an array (1 entry per sequence) with the freq of the seq in the entry
	"""

	res={}
	studyreads={}
	for idx,cseq in enumerate(seqs):
		if info=='samples':
			# Get the read vector for the samples
			pv=GetSeqVec(db,cseq)
			if len(pv)==0:
				continue
			nz=np.where(pv>db.MINFREQ)
			totalappear=len(nz[0])
			Debug(2,"--------Seq",cseq)
			Debug(3,"Found in",totalappear,"Samples")
			Debug(3,"Fraction=",float(totalappear)/len(pv))
			for csamp in nz[0]:
				res.setdefault(csamp,np.zeros(len(seqs)))[idx]=pv[csamp]
		elif info=='studies':
			totappear,numstudies,allstudies,studysamples,totdbsamples=GetSeqInfo(db,cseq)
			if totappear>0:
				sres=studysamples.items()
				vlens=[]
				for cv in sres:
					cstudy=cv[0]
					cnumreads=cv[1]
					if not cstudy in studyreads:
						studyreads[cstudy]=SamplesInStudy(db,cstudy)
					totsamps=studyreads[cstudy]
					vlens.append(float(len(cnumreads))/len(totsamps))
				sv,si=isort(vlens,reverse=True)
				for cind in si:
					studyid=sres[cind][0]
					if vlens[cind]>0.25:
						res.setdefault(studyid,np.zeros(len(seqs)))[idx]+=1
		elif info=='types':
			pv=GetSeqVec(db,cseq)
			if len(pv)==0:
				continue
			nz=np.where(pv>db.MINFREQ)
			for csamp in nz[0]:
				db.cur.execute("SELECT Field,Value FROM Maps WHERE SampleID=?",[int(csamp)+1])
				allvals=db.cur.fetchall()
				matter='NA'
				host='NA'
				for cv in allvals:
					if cv[0].lower()=='env_matter':
						matter=cv[1]
					if cv[0].lower()=='host_taxid':
						host=cv[1]
				res.setdefault(matter+'-'+host,np.zeros(len(seqs)))[idx]+=1

		else:
			Debug(9,"info type not supported",info)
			return False
	return res



#####################################
# database creation/addition methods
#####################################
def CreateTables(dbfilename,areyousure='no'):
	"""
	Create the database tables
	NOTE: will delete old database!!!!
	input:
	dbfilename - name of the database file to create
	areyousure - must be 'yes' in order to create the database

	output:
	db - the database structure
	"""

	assert areyousure=='yes'
	db=dbstart(dbfilename)
	db.cur.execute("DROP TABLE IF EXISTS Samples")
	db.cur.execute("CREATE TABLE Samples(SampleID INTEGER PRIMARY KEY AUTOINCREMENT,StudyID INT,SampleName TEXT,FileName TEXT,TotalReads INTEGER)")
	db.cur.execute("DROP TABLE IF EXISTS Reads")
	db.cur.execute("CREATE TABLE Reads(SeqID INT NOT NULL,SampleID INT NOT NULL,Reads DOUBLE NOT NULL)")
	db.cur.execute("DROP TABLE IF EXISTS Maps")
	db.cur.execute("CREATE TABLE Maps(UID INTEGER PRIMARY KEY AUTOINCREMENT,SampleID INTEGER,Field TEXT,Value TEXT)")
	db.cur.execute("DROP TABLE IF EXISTS Experiments")
	db.cur.execute("CREATE TABLE Experiments(StudyID INTEGER,ExpName TEXT,MapFileName TEXT)")
	db.cur.execute("DROP TABLE IF EXISTS Sequences")
	db.cur.execute("CREATE TABLE Sequences(SeqID INTEGER PRIMARY KEY AUTOINCREMENT,Sequence TEXT)")

	# create the indexes
	# note we do not create the read table indexes since it slows insertion down
	# can create later by using:
	# CreateReadInds()
	db.cur.execute("CREATE INDEX SeqIDInd ON Sequences (SeqID)")
	db.cur.execute("CREATE INDEX SeqInd ON Sequences (Sequence)")
	db.cur.execute("CREATE INDEX MapSampleInd ON Maps (SampleID)")
	db.cur.execute("CREATE INDEX SampleStudyInd ON Samples (StudyID)")
	db.cur.execute("CREATE INDEX SampleNameInd ON Samples (SampleName)")

	db.con.commit()
	return db


def CreateReadInds(db):
	"""
	create the read indexes (after db is ready - to make insertions faster)
	"""
	Debug(6,"Creating read sequence index")
	db.cur.execute("CREATE INDEX ReadSeqInd ON Reads (SeqID)")
	Debug(6,"Creating read sample index")
	db.cur.execute("CREATE INDEX ReadSampleInd ON Reads (SampleID)")
	db.con.commit()
	Debug(6,"Done creating indexes")



def GetSuidFromSampleID(db,sampleid,expid,addnew=False):
	"""
	get the suid (sample unique identifier) from the Samples table in the database
	if it does not exist, create it
	input:
	sampleid - the name of the sample (from mapping file)
	expid - the studyid from where the sample originated can also do without expid(=0)
	addnew - True to add the new sample to the table if not found, false to not add and return False
	output:
	suid - the sample id or False if not added
	isnew - True if a new sample, False if already existed
	"""

	isnew=True
	if expid>0:
		db.cur.execute("SELECT SampleID FROM Samples WHERE SampleName = ? AND StudyID = ?",(sampleid,expid))
	else:
		db.cur.execute("SELECT SampleID FROM Samples WHERE SampleName = ?",sampleid)
	res=db.cur.fetchone()
	if res:
		isnew=False
		suid=res[0]
		Debug(0,"Found SampleID",sampleid,"SUID=",suid)
		# see we don't have a double entry for same sampleid
		res=db.cur.fetchone()
		if res:
			self.Debug(10,"SampleID found twice!!!")
	else:
		if addnew:
			db.cur.execute("INSERT INTO Samples (SampleName,StudyID,TotalReads) VALUES (?,?,?)",(sampleid,expid,0))
			suid=db.cur.lastrowid
			Debug(0,"SampleID",sampleid,"Not found so added, new SUID=",suid)
		else:
			Debug(6,"Sample not found and not added",sampleid)
			suid=False
	return suid,isnew


def AddSample(db,sampleid,seqid,count):
	"""
	add count reads to position sampleid,seqid
	important: if it already exists, it is overwritten!!!!
	input:
	db
	sampleid - from the Samples table
	seqid - int (from the Sequences table)
	count - the fraction it was observed
	"""
	db.cur.execute("INSERT INTO Reads (SeqID,SampleID,Reads) VALUES(?,?,?)",(seqid,sampleid,count))
	Debug(0,"Added seqid ",seqid,"to sample",sampleid,"counts",count)



def AddFileToDB(db,expid,filename):
	"""
	load a deblurred sample file into the database
	File needs to be .ref.fa format
	(fasta, XXX;size=SIZE;\nSequence (one line))
	and one line per sequence (after the > line - not split)
	input:
	db
	expid - from where the sample is taken - integer
	sampleid - the name of the sample
	filename - file containing the .ref.fa of the sample
	"""
	# get the total number of reads
	totalreads=0
	Debug(2,"* scanning total reads from file",filename)
	afile=open(filename,'rU')
	for cline in afile:
		# get the size from the mark duplicaed in usearch by looking for ;size=NNN;
		sizematch=re.search('(?<=;size=)\d+',cline)
		numReads=int(sizematch.group(0))
		totalreads += numReads
		seqline=afile.next().strip()
	afile.close()
	Debug(2,"total reads",totalreads)
	if totalreads==0:
		Debug(3,"*** no reads in file",filename)
		return()
	if totalreads<self.MINEXPREADS:
		Debug(2,"*** not enough reads in file",filename,' it was',totalreads)
		return()

	# get sampleid from file name
	# just remove the '.fasta.ref.fa'
	sampleid=filename[:-13]
	sampleid=basename(sampleid)
	# get the sampleUID and upate total reads
	suid,isnew=GetSuidFromSampleID(db,sampleid,expid,addnew=True)
	db.cur.execute("UPDATE Samples SET TotalReads=? WHERE SampleName=? AND StudyID=?",(totalreads,suid,expid))
	if not isnew:
		Debug(6,"Sample already exists - deleting read file entries")
		db.cur.execute("DELETE FROM Reads WHERE SampleID=?",[suid])

	# and now read the sequences
	afile=open(filename,'rU')
	Debug(0,"inserting reads into table in file",filename)
	for cline in afile:
		sizematch=re.search('(?<=;size=)\d+',cline)
		numReads=int(sizematch.group(0))
		seqline=afile.next().strip()
		# get the greengenes ID and number of reads
		seqID=GetSeqID(db,seqline,insert=True)
		AddSample(db,suid,seqID,float32(numReads)/totalreads)

	afile.close()
	Debug(2,"finished reading",filename)
	db.con.commit()
	db.Debug(2,"Commited")


def GetStudyIDFromMap(mapfilename):
	"""
	get the study_id value from the mapping file if it exists
	otherwise, return False
	"""

	mf=open(mapfilename,'rU')
	reader=csv.DictReader(mf,delimiter='\t')
	cline=reader.next()
	cline=dictupper(cline)
	if 'STUDY_ID' in cline:
		studyid=cline['STUDY_ID']
	else:
		Debug(6,"STUDY_ID not found in map file",mapfilename)
		studyid=False
	mf.close()
	return studyid


def GetAllSampleInfo(db,fields):
	"""
	get the values of fields for every sample in database
	input:
	db - dbstruct
	fields : list of strings
		the list of fields to get the values for
	output:
	valsdict : dict of list of ints, keyed by concatanation of field values
		the positions for each value (concatanation of fields)
	"""
	db.cur.execute("SELECT MAX(SampleID) FROM Samples")
	res=db.cur.fetchone()
	numsamples=res[0]
	Debug(0,"Maximal SampleID is",numsamples)

	# get the values for each sample
	fields=listupper(fields)
	valsdict={}
	for csample in range(numsamples):
		db.cur.execute("SELECT Field,Value FROM Maps WHERE SampleID=?",[csample+1])
		allvals=db.cur.fetchall()
		cresdict={}
		for cres in allvals:
			cfield=cres[0].upper()
			cval=cres[1]
			cresdict[cfield]=cval
		cans=""
		for cfield in fields:
			if cfield in cresdict:
				cans+=cresdict[cfield]
			else:
				cans+="NF"
			cans+=';'
		if cans not in valsdict:
			valsdict[cans]=[]
		valsdict[cans].append(csample)
	return valsdict


def GetSamplesFromValues(db,fvdict):
	"""
	Get a list of samples that have values in their fields matching fv dict (all of the fields, any of the values in each field)
	input:
	db - dbstruct
	fvdict - dict
		keyd by field, value is a list of possible values for the field
	output:
	oksamps : list of int
		positions (in samples vector) that for all fields have one of the values matching. sampleids are the int+1
	"""
	# get # of reads
	db.cur.execute("SELECT MAX(SampleID) FROM Samples")
	res=db.cur.fetchone()
	numsamples=res[0]
	Debug(0,"Maximal SampleID is",numsamples)

	# convert keys and values to upper case
	fvdict=dictupper(fvdict)
	for k,v in fvdict.items():
		fvdict[k]=listupper(fvdict[k])

	# take all reads for the seqID
	Debug(0,"Getting reads")
	oksamps=[]
	totfound=0
	totnotfound=0
	for csample in range(numsamples):
		db.cur.execute("SELECT Field,Value FROM Maps WHERE SampleID=?",[csample+1])
		allvals=db.cur.fetchall()
		foundit={}
		for k in fvdict.keys():
			foundit[k]=False
		for cres in allvals:
			if not cres[0]:
				continue
			if not cres[1]:
				continue
			cfield=cres[0].upper()
			cval=cres[1].upper()
			if cfield in fvdict:
				if cval in fvdict[cfield]:
					foundit[cfield]=True
		foundall=True
		for v in foundit.values():
			foundall=foundall and v
		if foundall:
			oksamps.append(csample)
			totfound+=1
		else:
			totnotfound+=1
	Debug(6,"Found %d matching samples, %d non matching" % (totfound,totnotfound))
	return oksamps


def AddMap(db,experimentname,mapfilename,studyid=False,deleteifpresent=False,samplesadded=False):
	"""
	Add a tsv mapping file to the database db mapping file and experiment tables

	input:
	db - the database
	experimentname - name of the experiment (for the experiment table)
	mapfilename - name of the mapping file to add
	studyid - the id of the study or false to get from mapping file
	"""

	Debug(1,"Adding experiment to mapping database")
	# get the studyid
	if not studyid:
		studyid=GetStudyIDFromMap(mapfilename)
		if not studyid:
			raise ValueError("study_id not supplied and not in mapping file")
			return()

	# test if studyid already in database:
	db.cur.execute("SELECT * FROM Experiments WHERE StudyID=?",[studyid])
	res=db.cur.fetchone()
	if res:
		if not deleteifpresent:
			Debug(10,"Experiment",experimentname,"from mapfile",mapfilename,"already in database. id=",studyid)
			return()
		# if present and deleteifpresent - delete it!
		else:
			db.cur.execute("DELETE FROM Experiments WHERE StudyID=?",[studyid])
			Debug(6,"Study already in database - deleted it",studyid)

	# add study to Experiments table
	db.cur.execute("INSERT INTO Experiments (StudyID,ExpName,MapFileName) VALUES (?,?,?)",(studyid,experimentname,mapfilename))
	Debug(1,"Added experiment to experiments table")

	mf=open(mapfilename,'rU')
	reader=csv.DictReader(mf,delimiter='\t')
	for cline in reader:
		cline=dictupper(cline)
		# get the sampleids and the studyid
		try:
			sampleid=cline['#SAMPLEID']
			Debug(0,"SampleID",sampleid,"found in map",mapfilename)
		except:
			Debug(10,"#SampleID not found in map file",mapfilename)
			return()
		# if we have a list of samples with enough reads (from addbiomtable) and this sample didn't have enough reads - don't add it
		if samplesadded:
			if sampleid not in samplesadded:
				continue
		suid,isnew=GetSuidFromSampleID(db,sampleid=sampleid,expid=studyid,addnew=True)
		if not isnew:
			Debug(1,"Sample already exists - deleting mapping file entries")
			db.cur.execute("DELETE FROM Maps WHERE SampleID=?",[suid])
		for (field,val) in cline.items():
			db.cur.execute("INSERT INTO Maps (SampleID,Field,Value) VALUES (?,?,?)",(suid,field,val))
	mf.close()
	db.con.commit()
	Debug(1,"Added experiment to mapping database")


def AddBiomToDb(db,studyid,tablename,seqlength=False):
	"""
	load a deblurred biom table into the database
	input:
	db
	studyid - experiment for the biom table (studyid or other int)
	tablename - the experiment biom table file name (needs to be deblurred!)
	seqlength - False to use the actual lengths, >0 to clip all sequences to that length

	output:
	samplesadded - dict keyed with the samplenames of samples added to the database (i.e. >MINEXPREADS reads) value is the sampleid
	"""

	# minimal total number of reads per experiment in order to add to database
	READTHRESHABS=4
	# minimal mean number of reads per sample in order to add to database
	READTHRESHPERC=4.0/(10000*100)
	# minimal number of reads per sample in order to add sample to database
	MINEXPREADS=1000

	# load the biom table
	table = biom.load_table(tablename)
	samples=table.ids(axis='sample')
	seqs=table.ids(axis='observation')
	data=table.matrix_data.todense().A
	totreads=np.sum(data,axis=0)
	tottotreads=np.sum(totreads)
	samplesadded={}
	Debug(2,"Loaded table %s, %d sequences, %d samples" % (tablename,len(seqs),len(samples)))
	for idx,csamp in enumerate(samples):
		if totreads[idx]<MINEXPREADS:
			Debug(1,"not enough reads in sample %s (%d)" % (csamp,totreads[idx]))
			continue
		suid,isnew=GetSuidFromSampleID(db,csamp,studyid,addnew=True)
		samplesadded[csamp]=suid
		# and delete reads if already in the database
		if not isnew:
			db.cur.execute("SELECT TotalReads FROM Samples WHERE SampleID = ?",[suid])
			res=db.cur.fetchone()
			if res>0:
				Debug(6,"Sample already exists - deleting read file entries")
				db.cur.execute("DELETE FROM Reads WHERE SampleID=?",[suid])
		# and update the total reads
		db.cur.execute("UPDATE Samples SET TotalReads=? WHERE SampleID=?",(totreads[idx],suid))

	seqdict={}
	for idx,cseq in enumerate(seqs):
		if seqlength:
			# trim the sequence
			cseq=cseq[:seqlength]
			# and make sure it is in upper case
			cseq=cseq.upper()
		if cseq not in seqdict:
			seqdict[cseq]=[]
		seqdict[cseq].append(idx)

	Debug(2,"After trimming to len %s, %d sequences remaining" % (str(seqlength),len(seqdict)))
	for k,v in seqdict.items():
		if len(v)>1:
			allfreq=np.sum(data[v,:],axis=0)
		else:
			allfreq=data[v,:].flatten()
		Debug(1,"V is %s (%d)" % (v,len(v)))
		Debug(1,"length of allfreq is %d" % np.size(allfreq,0))
		# test if we have not enough reads don't save this sequence
		ctot=np.sum(allfreq)
		if ctot<READTHRESHABS or ctot<tottotreads*READTHRESHPERC:
			Debug(1,"Not enough reads for sequence %s (%f)" % (k,ctot))
			continue

		seqID=GetSeqID(db,k,insert=True)
		for sampidx,csamp in enumerate(samples):
			Debug(0,"Sample %s" % csamp)
			if csamp not in samplesadded:
				continue
			if allfreq[sampidx]==0:
				continue
			AddSample(db,samplesadded[csamp],seqID,np.float32(allfreq[sampidx])/totreads[sampidx])
			Debug(1,"added sequence %s in sample (%s)" % (k,csamp))

	Debug(2,"added biom table",tablename)
	Debug(1,"deleting samples with not enough reads")
	db.con.commit()
	Debug(2,"Commited")
	return samplesadded
