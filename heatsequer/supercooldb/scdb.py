#!/usr/bin/env python

"""
amnonscript
heatsequer
supercooldb.py

the sql cooldb implementation
"""

__version__ = "0.1"


from ..utils.amnonutils import Debug,dictupper,listupper,delete
from ..utils.oboparse import Parser
from ..utils.ontologygraph import ontologysubtreeids,ontologytotree,getnodeparents,ontotreetonames
from ..experiment.expclass import getheatsequerdir

import numpy as np
import matplotlib.pyplot as plt
import csv
import sqlite3
from collections import defaultdict
# for debugging - use XXX()
from pdb import set_trace as XXX
import datetime
import pickle
import requests
import os


class scdbstruct:
	def __init__(self):
		# the filename used for the database
		self.dbfile=''
		# the ontology dict (key is name and value is id)
		self.ontology={}
		# the ontology from id dict (key is id, and value is first name with this id) - used for save
		self.ontologyfromid={}
		# the dict of ontology graphs (load with loadontotree)
		self.ontodict={}
		# the names of the ontology files used:
		self.ontologyfiles=['/Users/amnon/Databases/ontologies/doid.obo','/Users/amnon/Databases/ontologies/envo.obo','/Users/amnon/Databases/ontologies/uberon.obo','/Users/amnon/Databases/ontologies/efo.obo','/Users/amnon/Databases/ontologies/po.obo','/Users/amnon/Databases/ontologies/gaz.obo']
		# the database server url
		# self.dburl='http://localhost:5000'
#		self.dburl='http://amnonim.webfactional.com/scdb:29708'
		self.dburl='http://amnonim.webfactional.com/scdb_main'
#		self.dburl='http://amnonim.webfactional.com/scdb_develop'


def addontology(scdb,ontology,ontoprefix='',namelist={}):
	"""
	add an obo ontology file to scdb
	input:
	scdb : scdbstruct
		from scdbstart()
	ontology : str
		name of the obo ontology file to add
	ontoprefix : str
		the ontology prefix (i.e. ENVO) to show at end of each string, or '' for autodetect (default)
	namelist : dict of ids
		if non-empty, keep only items with ids from namelist (get namelist from ontologysubtreeids() )

	output:
	scdb :scdbstruct
		with the new ontology items added to scdb.ontology dict
	"""
	parser=Parser(open(ontology))
	for citem in parser:
		cid=citem.tags["id"][0]
		if len(ontoprefix)==0:
			tt=cid.split(':')
			if len(tt)>1:
				ontoprefix=tt[0]
				Debug(2,'Found ontology prefix %s' % ontoprefix)
		if namelist:
			if cid not in namelist:
				continue
		names=[]
		if "name" in citem.tags:
			names.extend(citem.tags["name"])
			origname=citem.tags["name"][0]
			scdb.ontologyfromid[cid]=origname
		else:
			Debug(6,"ontology item id %s does not have a name" % cid)
			origname="NA"
		if "synonym" in citem.tags:
			names.extend(citem.tags["synonym"])

		for cname in names:
			Debug(1,"%s %s" % (cname,cid))
			oname=cname+' :'+ontoprefix
			if cname!=origname:
				oname+='('+origname+')'
			if oname in scdb.ontology:
				Debug(1,"name %s id %s already in ontology list for id %s" % (oname,cid,scdb.ontology[oname]))
			scdb.ontology[oname]=cid

	return scdb


def loaddbonto(db,ontofile=None,ontofromidfile=None):
	"""
	load the pickled ontologies to the scdb structure
	input:
	db : from scdbstart
	ontofile : str
		name of the ontologies term file (from saveontologies) or None for default location
	ontofromidfile : str
		name of the ontologies reverse dict file (from saveontologies) or None for default location

	output:
	db : scdbstruct
		with the loaded ontologies fields
	"""
	if ontofile is None:
		ontofile=os.path.join(getheatsequerdir(),'db/ontology.pickle')
	if ontofromidfile is None:
		ontofromidfile=os.path.join(getheatsequerdir(),'db/ontologyfromid.pickle')
	Debug(6,'loading ontology pickles')
	Debug(6,'Files %s and %s' % (ontofile,ontofromidfile))
	db.ontology=pickle.load(open(ontofile,'rb'))
	db.ontologyfromid=pickle.load(open(ontofromidfile,'rb'))
	Debug(6,'ontologies loaded')
	return db


def loadontologies(scdb,pickleit=True,ontologies=[]):
	"""
	load the ontologies into the scdb class
	input:
	scdb : scdbstruct
		from scdbstart()
	pickleit : bool
		True (default) to pickle the loaded ontologies to default location, False to not pickle
	ontologylist : list of str
		names of the obo ontology files to load or empty to use the default files (from scdb.ontologyfiles)
	"""
	if not ontologies:
		ontologies=scdb.ontologyfiles

	scdb.ontology={}
	scdb.ontologyfromid={}
	for contology in ontologies:
		addontology(scdb,contology)
	if pickleit:
		saveontologies(scdb)


def saveontologies(scdb,ontofile=None,ontofromidfile=None):
	"""
	save the ontologies to pickle files.
	use after loadontologies()
	"""
	if ontofile is None:
		ontofile=os.path.join(getheatsequerdir(),'db/ontology.pickle')
	if ontofromidfile is None:
		ontofromidfile=os.path.join(getheatsequerdir(),'db/ontologyfromid.pickle')

	fl=open(ontofile,'wb')
	pickle.dump(scdb.ontology,fl,protocol=2)
	fl.close()
	fl=open(ontofromidfile,'wb')
	pickle.dump(scdb.ontologyfromid,fl,protocol=2)
	fl.close()


def dbstart(dbname="db/supercooldb.db"):
	'''
	start the database structure and connect to database
	input:
	dbname : string
		the name of the database to connect to
	output:
	db : dbstruct
		the database variable
	'''

	scdb=scdbstruct()
#	scdb=dbconnect(scdb,dbname)
	return scdb


def addexpdata(db,data,studyid=None):
	"""
	add new data entries (for a new study)
	input:
	db : scdbstruct
		from dbstart()
	data : list of tuples (Type:Value)
		a list of tuples of (Type,Value) to add to Data table (i.e. ("PUBMEDID","322455") etc)
	studyid : list of int
		the ids in which this study appears (from finddataid)

	output:
	suid : int
		the value of DataID for the new study (from Data table)
	"""
	# we need to get a new identifier for all entries in the study
	# there should be a more elegant way to do it
	Debug(2,"addexpdata for %d enteries" % len(data))
	if studyid is None:
		# add new study
		Debug(2,"addexpdata for a new study")
	else:
		Debug(2,'addexpdata for existing study %d' % studyid)
	rdata={}
	rdata['expId']=studyid
	rdata['details']=data
	res=requests.post(db.dburl+'/experiments/add_details',json=rdata)
	if res.status_code==200:
		newid=res.json()['expId']
		Debug(2,'experiment added. id is %d' % newid)
		return newid
	else:
		Debug(8,'error adding experiment. msg: %s' % res.content)
		return None


def addannotations(db,expid,sequences,annotationtype,annotations,submittername='NA',description='',method='',primerid=0,agenttype='HeatSequer',private='n'):
	"""
	Add a new manual curation to the database
	input:
	db : scdbstruct
		from dbstart()
	expid : int or dict of (Type:Value)
		if int - the value of DataID from Data table, otherwise a list of (Type,Value) tuples to add to Data table
	sequences : list of ACGT
		the sequences to curate
	annotationtype : str
		the curation type (COMMON,DIFFEXP,CONTAM,HIGHFREQ,PATHOGEN)
	annotations : list of Type,Value
		The curations to add to the CurationList table (Type,Value)
	submittername : str
		Name of the submitter (first,last) or NA
	description : str
		text description of the curation entry (i.e. "lower in whole wheat pita bread")
	method : str
		text description of how the curation was detected - only if needed
	primerid : int
		the PrimerID from Primers table of the sequences (usually 1 - the V4 515F,806R)
	agenttype : str
		the program submitting the curation
	private : str (optional)
		'n' (default) or 'y'

	output:
	curationid : int
		the CurationID (in Curations table) of the new curation, or 0 if not added
	data : int
		the DataID from Data table
	"""
	Debug(2,"addannotation - %d sequences" % len(sequences))
	if len(sequences)==0:
		Debug(6,"No sequences to annotate!")
		return 0,0
	if len(annotations)==0:
		Debug(6,"No annotations to add. still adding...")
	if not type(expid) is int:
		Debug(6,"looking for studyid %s in data" % expid)
		expid=addexpdata(db,expid)
		if expid is None:
			Debug(8,'problem adding new experiment data')
			return 0,0

	# add the curation
	rdata={}
	rdata['expId']=expid
	rdata['sequences']=sequences
	rdata['region']=primerid
	rdata['annotationType']=annotationtype
	rdata['method']=method
	rdata['agentType']=agenttype
	rdata['description']=description
	rdata['private']=private
	rdata['annotationList']=annotations

	res=requests.post(db.dburl+'/annotations/add',json=rdata)
	if res.status_code==200:
		newid=res.json()['annotationId']
		Debug(1,"Finished adding experiment id %d annotationid %d" % (expid,newid))
		return res,newid
	Debug(8,'problem adding annotations for experiment id %d' % expid)
	Debug(8,res.content)
	return 0,0


def finddataid(db,datamd5='',mapmd5='',getall=False):
	"""
	find the data id for the data/map md5 (which are calculated on load)
	note the md5s don't change following filtering/normalization/etc... - only the original data
	input:
	scdb : from startdb()
	datamd5 : str
		from Experiment.datamd5
	mapmd5 : str
		from Experiment.mapmd5
	getall : bool (optional)
		False (default) to get only 1st id, True to get a list of all

	output:
	expids: int (if getall=False - default) or list of int (if getall=True)
		an id or a list of ids of matching dataID indices (or None if no match)
	"""
	Debug(1,'findexpid for datamd5 %s mapmd5 %s' % (datamd5,mapmd5))
	details=[]
	if datamd5:
		details.append(['DataMD5',datamd5])
	if mapmd5:
		details.append(['MapMD5',mapmd5])
	if len(details)==0:
		Debug(6,'Error. MapMD5 and DataMD5 both missing from finddataid')
		return None

	rdata={}
	rdata['details']=details
	res=requests.get(db.dburl+'/experiments/get_id',json=rdata)
	if res.status_code==200:
		expids=res.json()['expId']
		if not getall:
			if len(expids)>1:
				Debug(6,'Problem. Found %d matches for data' % len(expids))
			Debug(2,'Found study id %d' % expids[0])
			return expids[0]
		Debug(2,"Found %d matches to data" % len(expids))
		return expids
	Debug(8,'Error getting expid from details')
	return None


def getexperimentinfo(db,expid):
	"""
	get the information about a given study dataid
	input:
	db : from dbstart()
	dataid : int
		The dataid on the study (DataID field)

	output:
	info : list of (str,str,str)
		list of tuples for each entry in the study:
		type,value,descstring about dataid
		empty if dataid not found
	"""
	Debug(1,'get experiment details for expid %d' % expid)
	rdata={}
	rdata['expId']=expid
	res=requests.get(db.dburl+'/experiments/get_details',json=rdata)
	if res.status_code==200:
		details=res.json()['details']
		Debug(2,'Found %d details for experiment %d' % (len(details),expid))
		return details
	return []


def getexpannotations(db,expid):
	"""
	get the list of annotations for study studyid

	input:
	db : from dbstart()
	expid : int
		The dataid of the study

	output:
	info: list of str
		the list of curations for this study (1 item per curation)
	"""
	Debug(1,'get experiment annotations for expid %d' % expid)
	rdata={}
	rdata['expId']=expid
	res=requests.get(db.dburl+'/experiments/get_annotations',json=rdata)
	if res.status_code!=200:
		Debug(6,'error getting annotations for experiment %d' % expid)
		return []
	annotations=res.json()['annotations']
	Debug(2,'Found %d annotations for experiment %d' % (len(annotations),expid))
	# make it into a nice list of str
	info=[]
	for cann in annotations:
		cstr='date:%s description:%s user:%s private:%s' % (cann['date'],cann['description'],cann['userid'],cann['private'])
		info.append(cstr)
	return info


def getannotationseqs(db,annotationid):
	"""
	get the list of sequences to which an annotation relates

	input:
	db : from dbstart()
	annotationid : int
		the unqiue id of the annotation (annotationid in the annotation details)

	output:
	seqids: list of int
		list of sequences to which the annotation is attached
	"""
	Debug(1,'get annotationseqs for annotationid %d' % annotationid)
	rdata={}
	rdata['annotationid']=annotationid
	res=requests.get(db.dburl+'/annotations/get_sequences',json=rdata)
	if res.status_code!=200:
		Debug(6,'error getting sequences for annotation %d' % annotationid)
		Debug(6,res.content)
		return []
	seqids=res.json()['seqids']
	Debug(2,'Found %d sequences for annotationid %d' % (len(seqids),annotationid))
	return seqids


def getseqannotations(db,sequence):
	"""
	Get the manual curations for a sequence

	input:
	db : from scdbstart()
	sequence : str (ACGT)

	output:
		curs : list of list of (curation dict,list of [Type,Value] of curation details)
	"""
	Debug(1,'get sequence annotations for sequence %s' % sequence)
	rdata={}
	rdata['sequence']=sequence
	print('***'+db.dburl+'/sequences/get_annotations')
	res=requests.get(db.dburl+'/sequences/get_annotations',json=rdata)
	if res.status_code!=200:
		Debug(6,'error getting annotations for sequence %s' % sequence)
		return []
	print(res.json())
	annotations=res.json()['annotations']
	Debug(2,'Found %d annotations for sequence %s' % (len(annotations),sequence))
	return annotations


def getannotationstrings(db,sequence):
	"""
	get a nice string summary of a curation

	input:
	db : from scdbstart()
	sequence : str (ACGT)

	output:
	shortdesc : list of (dict,str) (annotationdetails,annotationsummary)
		a list of:
			annotationdetails : dict
				'annotationid' : int, the annotation id in the database
				'annotationtype : str
				...
			annotationsummary : str
				a short summary of the annotation
	"""
	shortdesc=[]
	annotations=getseqannotations(db,sequence)
	for cann in annotations:
		annotationdetails=cann
#		annotationdetails['annotationid']=cann['annotationid']
#		for k,v in cann.items():
#			annotationdetails[k]=v
#		annotationdetails['annotationtype']=cann['annotationtype']
		cdesc=''
		if cann['description']:
			cdesc+=cann['description']+' ('
		if cann['annotationtype']=='diffexp':
			chigh=[]
			clow=[]
			call=[]
			for cdet in cann['details']:
				if cdet[0]=='all':
					call.append(cdet[1])
					continue
				if cdet[0]=='low':
					clow.append(cdet[1])
					continue
				if cdet[0]=='high':
					chigh.append(cdet[1])
					continue
			cdesc+=' high in '
			for cval in chigh:
				cdesc+=cval+' '
			cdesc+=' compared to '
			for cval in clow:
				cdesc+=cval+' '
			cdesc+=' in '
			for cval in call:
				cdesc+=cval+' '
		elif cann['annotationtype']=='isa':
			cdesc+=' is a '
			for cdet in cann['details']:
				cdesc+='cdet,'
		elif cann['annotationtype']=='contamination':
			cdesc+='contamination'
		else:
			cdesc+=cann['annotationtype']+' '
			for cdet in cann['details']:
				cdesc=cdesc+' '+cdet[1]+','
		shortdesc.append( (annotationdetails,cdesc) )
	return shortdesc

###################################################


def getseqcurations(db,sequence):
	"""
	Get the manual curations for a sequence

	input:
	db : from scdbstart()
	sequence : str (ACGT)

	output:
		curs : list of list of (curation dict,list of [Type,Value] of curation details)
	"""
	curs=[]
	seqid=getseq(db,sequence,insert=False)
	if seqid==0:
		Debug(2,'Sequence not found')
		return curs
	curids=getseqcurationids(db,seqid)
	for cid in curids:
		cdat=select_column_and_value(db.cur,"SELECT * FROM Curations WHERE CurationID = ?",[cid])
		if cdat=={}:
			Debug(8,'no curation found for curationid %d' % cid)
			continue
		ccur=cdat
		Debug(2,cdat)
		ccurdetails=[]
		curation=db.cur.execute('SELECT Type,Value from CurationList WHERE CurationID = ?',[cid])
		for ccuration in curation.fetchall():
			Debug(2,ccuration)
			ccurdetails.append([ccuration[0],ccuration[1]])
		curs.append([ccur,ccurdetails])
	return curs


def getcurationstrings(db,sequence):
	"""
	get a nice string summary of a curation

	input:
	db : from scdbstart()
	sequence : str (ACGT)

	output:
	shortdesc : list of str
		a short summary of the curations (1 item per curation)
	"""
	shortdesc=[]
	curs=getseqcurations(db,sequence)
	for ccur in curs:
		cdesc=''
		curinfo=ccur[0]
		curdetails=ccur[1]
		if curinfo['Description']:
			cdesc+=curinfo['Description']+' ('
		if curinfo['CurType']=='COMMON':
			cdesc+='common in '
			for cdet in curdetails:
				cdesc+=cdet[1]+' '
		elif curinfo['CurType']=='HIGHFREQ':
			cdesc+='high freq in '
			for cdet in curdetails:
				cdesc+=cdet[1]+' '
		elif curinfo['CurType']=='DIFFEXP':
			chigh=[]
			clow=[]
			call=[]
			for cdet in curdetails:
				if cdet[0]=='ALL':
					call.append(cdet[1])
					continue
				if cdet[0]=='LOW':
					clow.append(cdet[1])
					continue
				if cdet[0]=='HIGH':
					chigh.append(cdet[1])
					continue
			cdesc+=' high in '
			for cval in chigh:
				cdesc+=cval+' '
			cdesc+=' compared to '
			for cval in clow:
				cdesc+=cval+' '
			cdesc+=' in '
			for cval in call:
				cdesc+=cval+' '
		else:
			Debug(2,'unknown curation %s.' % curinfo['CurType'])
		shortdesc.append(cdesc)
	return shortdesc


def select_column_and_value(cur, sql, parameters=()):
	"""
	get a dict with field as key, value as values for an sql query
	"""
	execute = cur.execute(sql, parameters)
	fetch = execute.fetchone()

	if fetch is None:
		return {}

	return {k[0]: v for k, v in list(zip(execute.description, fetch))}


def createontologytree(db,ontologies=[],outname=None):
	"""
	load the ontology tree graphs into the scdb and store them in a pickle dict

	input:
	db : from scdbstart()
	ontologies : list of str
		list of obo ontologies to use or empty to use the default files (db.ontologyfiles)
	outname : str
		name of the output pickle file or None for default location

	output:
	db - scdbstruct
		with the ontodict field added
	"""
	if outname is None:
		outname=os.path.join(getheatsequerdir(),'db/ontologygraph.pickle')

	if not ontologies:
		ontologies=db.ontologyfiles

	if not db.ontologyfromid:
		db=loaddbonto(db)

	ontodict={}
	for conto in ontologies:
		Debug(6,'Processing ontology %s' % conto)
		g=ontologytotree(conto)
		g=ontotreetonames(g,db.ontologyfromid)
		ontodict[conto]=g
	Debug(6,'ontologies loaded. saving to pickel %s' % outname)
	fl=open(outname,'wb')
	pickle.dump(ontodict,fl,protocol=2)
	Debug(6,'ontologies pickled')
	db.ontodict=ontodict
	return db



def loadontotrees(db,ontopickle=None):
	"""
	load the ontology dict pickle file

	input:
	db : from scdbstart()
	ontopickle : str
		the pickled ontology dict filename (from createontologytree()) or None for default location

	output:
	db : scdbstruct
		with the ontology dict in ontodict
	"""
	if ontopickle is None:
		ontopickle=os.path.join(getheatsequerdir(),'db/ontologygraph.pickle')
	Debug(6,'loadding ontology trees')
	fl=open(ontopickle,'rb')
	db.ontodict=pickle.load(fl)
	Debug(6,'loaded %d trees' % len(db.ontodict))
	return db


def getontoparents(db,term):
	"""
	get all the parent terms (including original term) for a given ontology term.
	look in all ontology trees in db.ontodict

	input:
	db : from scdbstart()
	term : str
		the ontology term (lower case)

	output:
	terms : list of str
		all the parent terms of this item
	"""
	terms=[]
	for conto in db.ontodict.values():
		parents=getnodeparents(conto,term)
		terms+=parents
	terms=list(set(terms))
	return terms


def getcurationontologies(db,sequence):
	"""
	get the curation for items and all upstream (ontology-wise) items for a given sequence

	input:
	db : from scdbstart()
	sequence : str (ACGT)

	output:
	onto : dict of key=str and value dict of key=str and value=int
		a dictionary of ['up','down','contaminant'] of dict of ontology item and number of times total we see it
	"""
	if not db.ontodict:
		loadontotrees(db)
	onto={}
	curs=getseqcurations(db,sequence)
	for ccur in curs:
		curinfo=ccur[0]
		curdetails=ccur[1]
		ctype='other'
		lookseparate=False
		if curinfo['CurType']=='COMMON':
			ctype='up'
		elif curinfo['CurType']=='HIGHFREQ':
			ctype='up'
		else:
			lookseparate=True
		for cdet in curdetails:
			if lookseparate:
				if cdet[0]=='ALL':
					ctype='up'
				elif cdet[0]=='LOW':
					ctype='down'
				elif cdet[0]=='HIGH':
					ctype='up'
				else:
					ctype='other'
			if ctype not in onto:
				onto[ctype]={}
			ontoparents=getontoparents(db,cdet[1])
			for conto in ontoparents:
				if conto not in onto[ctype]:
					onto[ctype][conto]=1
				else:
					onto[ctype][conto]+=1
	return onto


def delete_annotation(db,annotationid):
	'''
	delete an annotation from the database.

	input:
	db :
	annotationid : int
		the annotationid to delete

	output:
	'''
	Debug(1,'delete annotation for annotatioid %d' % annotationid)
	rdata={}
	rdata['annotationid']=annotationid
	res=requests.post(db.dburl+'/annotations/delete',json=rdata)
	if res.status_code!=200:
		Debug(6,'error deleting annotationid %d' % annotationid)
		Debug(6,res.content)
		return []
	Debug(2,'Annotation %d deleted' % annotationid)




def get_experiment_annotations(db,exp):
	'''
	get annotations on all sequences in the experiment

	input:
	db:
	exp : Experiment

	output:
	'''
	Debug(1,'get experiment sequence annotations for %d sequences' % len(exp.seqs))
	rdata={}
	rdata['sequences']=exp.seqs
	res=requests.get(db.dburl+'/sequences/get_list_annotations',json=rdata)
	if res.status_code!=200:
		Debug(6,'error getting list annotations')
		Debug(6,res.content)
		return []
	return res.json()['seqannotations']


def convert_olddb_to_server(db,olddbfilename='db/supercooldb.db'):
	'''
	load the old local sqlite3 database to the server

	input:
	db
	olddbfilename : str
		the sqlite3 database filename
	'''
	con=sqlite3.connect(olddbfilename)
	Debug(7,"Connected to database")
	# and the cursor
	cur=con.cursor()

	# add all studies
	# get studyids
	Debug(7,'getting study ids')
	cur.execute("SELECT DataID FROM Data")
	allvals=cur.fetchall()
	Debug(7,'found %d details' % len(allvals))
	studyids=set()
	for cres in allvals:
		studyids.add(cres[0])
	Debug(7,'found %d unique ids' % len(studyids))
	Debug(7,'processing per study data')
	studyoldnew={}
	for cid in studyids:
		cur.execute("SELECT Type,Value FROM Data WHERE DataID = ?",[cid])
		allvals=cur.fetchall()
		Debug(7,'found %d details for studyid %d' % (len(allvals),cid))
		Debug(7,'adding study data')
		Debug(7,'looking for md5 in previous studies')
		for cres in allvals:
			if cres[0]=='DataMD5':
				datamd5=cres[1]
			if cres[0]=='MapMD5':
				mapmd5=cres[1]
		newid=finddataid(db,datamd5=datamd5,mapmd5=mapmd5,getall=False)
		if newid is not None:
			Debug(7,'already exists. id=%d' % newid[0])
		# newid = addexpdata(db,list(allvals))
		newid=0
		Debug(7,'added data new studyid %d' % newid)
		studyoldnew[cid]=newid
