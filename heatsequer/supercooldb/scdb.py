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

import numpy as np
import matplotlib.pyplot as plt
import csv
import sqlite3
from collections import defaultdict
# for debugging - use XXX()
from pdb import set_trace as XXX
import datetime
import pickle


class scdbstruct:
	def __init__(self):
		# the filename used for the database
		self.dbfile=''
		# the ontology dict (key is name and value is id)
		self.ontology={}
		# the ontology from id dict (key is id, and value is first name with this id) - used for save
		self.ontologyfromid={}



def addontology(scdb,ontology,ontoprefix=''):
	"""
	add an obo ontology file to scdb
	input:
	scdb : scdbstruct
		from scdbstart()
	ontology : str
		name of the obo ontology file to add
	ontoname : str
		the ontology prefix (i.e. ENVO) to show at end of each string, or '' for autodetect (default)

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


def loadontologies(scdb,pickleit=True,ontologies=['/Users/amnon/Databases/ontologies/doid.obo','/Users/amnon/Databases/ontologies/envo.obo','/Users/amnon/Databases/ontologies/uberon.obo','/Users/amnon/Databases/ontologies/efo.obo','/Users/amnon/Databases/ontologies/po.obo','/Users/amnon/Databases/ontologies/gaz.obo']):
	"""
	load the ontologies into the scfb class
	input:
	scdb : scdbstruct
		from scdbstart()
	pickleit : bool
		True (default) to pickle the loaded ontologies to default location, False to not pickle
	ontologylist : list of str
		names of the obo ontology files to load
	"""
	scdb.ontology={}
	scdb.ontologyfromid={}
	for contology in ontologies:
		addontology(scdb,contology)
	if pickleit:
		saveontologies(scdb)


def saveontologies(scdb,ontofile='/Users/amnon/Python/git/heatsequer/db/ontology.pickle',ontofromidfile='/Users/amnon/Python/git/heatsequer/db/ontologyfromid.pickle'):
	"""
	save the ontologies to pickle files.
	use after loadontologies()
	"""
	fl=open(ontofile,'wb')
	pickle.dump(scdb.ontology,fl)
	fl.close()
	fl=open(ontofromidfile,'wb')
	pickle.dump(scdb.ontologyfromid,fl)
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
	scdb=dbconnect(scdb,dbname)
	return scdb


def dbconnect(scdb,dbname="db/supercooldb.db"):
	"""
	connect to the database
	input:
	scdb : scdbstruct
		from scdbstart()
	dbname - the name of the database to connect to
	"""

	scdb.dbfile=dbname
	# the database connection
	Debug(1,"Connecting to database ",dbname)
	scdb.con=sqlite3.connect(scdb.dbfile)
	Debug(1,"Connected")
	# and the cursor
	scdb.cur=scdb.con.cursor()
	# test if the database file exists:
	try:
		scdb.cur.execute("SELECT count(*) FROM sqlite_master WHERE type='table' AND name='Curations'")
		istable=scdb.cur.fetchone()
		if istable[0]==0:
			Debug(9,"Can't find Reads table in database file %s" % dbname)
			scdb=False
	except:
		Debug(9,"Can't test database file %s" % dbname)
		scdb=False
	return scdb


def CreateTables(dbfilename='db/supercooldb.db',areyousure='no'):
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

	# Primers table - info about each primer set used
	# PrimerID - uniqueID - the key for other table
	# FPrimer - forward primer sequence ACGT (upper case)
	# RPrimer - reverse primer sequence ACGT (upper case)
	# Region - name of the region (i.e. V4, etc)
	db.cur.execute("DROP TABLE IF EXISTS Primers")
	db.cur.execute("CREATE TABLE Primers(PrimerID INTEGER PRIMARY KEY AUTOINCREMENT,FPrimer TEXT,RPrimer TEXT,Region TEXT)")
	db.cur.execute("CREATE INDEX PrimerIDInd ON Primers (PrimerID)")
	# add the V4 EMP primers
	db.cur.execute("INSERT INTO Primers (PrimerID,FPrimer,RPrimer,Region) VALUES (?,?,?,?)",(1,"GTGCCAGC[AC]GCCGCGGTAA","ATTAGA[AT]ACCC[CGT][AGT]GTAGTCC","V4"))

	# Data table - info about the experiment sequence data used
	# DataUniqueID - uniqueID - the key for other table
	# DataID - same for all data info for the same experiment (so can have multiple entries for the same experiment - paper,qiita,sra etc.)
	# also contains "DataMD5" and "MapMD5" with appropriate values
	# Type - the type of reference (SRA,DRYAD,QIITA,PUBMEDID,HTML)
	# Value - the value for the Type field selected (i.e. 10192,www.pita.com etc.)
	db.cur.execute("DROP TABLE IF EXISTS Data")
	db.cur.execute("CREATE TABLE Data(DataUniqueID INTEGER PRIMARY KEY AUTOINCREMENT,DataID INTEGER,Type TEXT,Value TEXT)")
	db.cur.execute("CREATE INDEX DataIDInd ON Data (DataID)")
	db.cur.execute("CREATE INDEX ValueInd ON Data (Value)")

	# Sequences table - info about each sequence with manual curation
	# SeqID - uniqueID - the key for other table
	# Sequence - ACGT (upper case)
	# Primers - link to PrimerID of Primers table for forward and reverse primer sequences
	# Length - lenght of the sequence
	# Taxonomy - the taxonomy string for the sequence
	db.cur.execute("DROP TABLE IF EXISTS Sequences")
	db.cur.execute("CREATE TABLE Sequences(SeqID INTEGER PRIMARY KEY AUTOINCREMENT,Sequence TEXT,PrimerID INTEGER, Length INTEGER, Taxonomy TEXT)")
	db.cur.execute("CREATE INDEX SeqIDInd ON Sequences (SeqID)")
	db.cur.execute("CREATE INDEX SeqInd ON Sequences (Sequence)")

	# SeqCurations table - manual curation ids for each sequenceid
	# SeqID - the sequence id (from Sequences table)
	# CurationID - the curation id (from Curations table)
	db.cur.execute("DROP TABLE IF EXISTS SeqCurations")
	db.cur.execute("CREATE TABLE SeqCurations(SeqID INTEGER,CurationID INTEGER)")

	# Curations table - manual curation data
	# CurationID - the curation id (from Curations). also used a identifier for all CurationList items for this curation
	# Date - the date when this curation was added (ISO8601 strings ("YYYY-MM-DD HH:MM:SS.SSS"))
	# SubmitterName - name of the user who added this curation (first,last) or NA
	# DataID - the identifier (for Data table) of all entries related to the dataset used for this analysis
	# Description - a text description of the observation (i.e. higher in pigs with IBD)
	# Method - how the annotation was achieved
	# SubmitterAgent - name of the program submitting the curation (e.g. HeatSequer)
	# CurType - the curation type ('DIFFEXP,COMMON,CONTAM,HIGHFREQ,PATHOGEN')
	db.cur.execute("DROP TABLE IF EXISTS Curations")
	db.cur.execute("CREATE TABLE Curations(CurationID INTEGER PRIMARY KEY AUTOINCREMENT,Date TEXT,SubmitterName TEXT,DataID INTEGER,Description TEXT,Method TEXT,SubmitterAgent TEXT,CurType TEXT)")
	db.cur.execute("CREATE INDEX CurationIDInd ON Curations (CurationID)")
	db.cur.execute("CREATE INDEX CurDataIDInd ON Curations (DataID)")

	# CurationList table - the list of curations (i.e. contamination, higher in, all, etc)
	# CurationID - for linking to other table, all entries with the same value are for the same curation
	# Type - type of curation (ALL - all bacteria in the experiment, HIGHERIN, LOWERIN, ISA)
	# Value - for the type (for ISA - COMMON, HIGHFREQ, CONTAMINATION, PATHOGEN, etc.)
	db.cur.execute("DROP TABLE IF EXISTS CurationList")
	db.cur.execute("CREATE TABLE CurationList(CurationID INTEGER, Type TEXT,Value TEXT)")
	db.cur.execute("CREATE INDEX CurationListIDInd ON CurationList (CurationID)")
	db.cur.execute("CREATE INDEX CurationListValue ON CurationList (Value)")

	db.con.commit()
	return db


def getseq(db,seq,primerid=1,insert=False):
	"""
	get the sequence id for the sequence seq in the database.
	if it doesn't exist - create it and return the new id

	input:
	db : scdbstruct
		from dbstart()
	seq : sequence (ACGT)
		the sequence to find/add
	primerid : int
		the PrimerID from Primers table of the sequenced region
	insert : bool
		True to insert the sequence if not in table, False to return 0 and not insert

	output:
	seqid : int
		the SeqID from the Sequences table for the sequence
	"""
	seq=seq.upper()
	db.cur.execute("SELECT SeqID FROM Sequences WHERE Sequence = ?",[seq])
	res=db.cur.fetchone()
	if res:
		seqid=res[0]
		Debug(0,"Found sequence %s,SeqID=%s" % (seq,seqid))
	else:
		if insert:
			Debug(0,"Sequence %s not found - creating in table" % seq)
			db.cur.execute("INSERT INTO Sequences (Sequence,PrimerID,Length) VALUES (?,?,?)",[seq,primerid,len(seq)])
			seqid=db.cur.lastrowid
			Debug(0,"SequenceID is",seqid)
		else:
			Debug(3,"Sequence %s not found in database" % seq)
			seqid=0
	return(seqid)



def adddata(db,data,studyid=[]):
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
	Debug(2,"adddata for %d enteries" % len(data))
	if len(studyid)==0:
		db.cur.execute("INSERT INTO Data (DataID,Type,Value) VALUES (?,?,?)",(0,"tmp","tmp"))
		suid=db.cur.lastrowid
		db.cur.execute("DELETE FROM Data WHERE DataUniqueID=?",[suid])
		Debug(2,"New adddata DataID: %d" % suid)
	else:
		suid=studyid[0]
		Debug(2,'adding new data to existing study id %d' % suid)
	for (k,v) in data:
		db.cur.execute("INSERT INTO Data (DataID,Type,Value) VALUES (?,?,?)",(suid,k,v))
	db.con.commit()
	Debug(2,"Data added")
	return suid


def addcuration(db,data,sequences,curtype,curations,submittername='NA',description='',method='',primerid=1,submitteragent='HeatSequer'):
	"""
	Add a new manual curation to the database
	input:
	db : scdbstruct
		from dbstart()
	data : int or dict of (Type:Value)
		if int - the value of DataID from Data table, otherwise a list of (Type,Value) tuples to add to Data table
	sequences : list of ACGT
		the sequences to curate
	curtype : str
		the curation type (COMMON,DIFFEXP,CONTAM,HIGHFREQ,PATHOGEN)
	curations : list of Type,Value
		The curations to add to the CurationList table (Type,Value)
	submittername : str
		Name of the submitter (first,last) or NA
	description : str
		text description of the curation entry (i.e. "lower in whole wheat pita bread")
	method : str
		text description of how the curation was detected - only if needed
	primerid : int
		the PrimerID from Primers table of the sequences (usually 1 - the V4 515F,806R)
	submitteragent : str
		the program submitting the curation

	output:
	curationid : int
		the CurationID (in Curations table) of the new curation, or 0 if not added
	data : int
		the DataID from Data table
	"""
	Debug(2,"addcuration - %d sequences" % len(sequences))
	if len(sequences)==0:
		Debug(6,"No sequences to annotate!")
		return 0,0
	if len(curations)==0:
		Debug(6,"No currations to add. still adding...")
	if not type(data) is int:
		Debug(6,"looking for studyid %s in data" % data)
		data=adddata(db,data)

	# add the curation
	cdatetime=datetime.datetime.now().replace(microsecond=0).isoformat()
	db.cur.execute("INSERT INTO Curations (Date,submittername,DataID,Description,Method,submitteragent,curtype) VALUES (?,?,?,?,?,?,?)",(cdatetime,submittername,data,description,method,submitteragent,curtype))
	curationid=db.cur.lastrowid
	# and add the curationlist entries
	Debug(1,"Adding %d curations" % len(curations))
	for (k,v) in curations:
		db.cur.execute("INSERT INTO CurationList (CurationID,Type,Value) VALUES (?,?,?)",(curationid,k,v))

	# add seqcuration entry for each sequence (create the sequence if doesn't exist)
	Debug(1,"Adding %d sequences to SeqCuration table" % len(sequences))
	for cseq in sequences:
		cseqid=getseq(db,cseq,primerid=primerid,insert=True)
		db.cur.execute("INSERT INTO SeqCurations (SeqID,CurationID) VALUES (?,?)",(cseqid,curationid))

	db.con.commit()
	Debug(1,"Finished adding")
	return curationid,data


def finddataid(db,datamd5='',mapmd5=''):
	"""
	find the data id for the data/map md5 (which are calculated on load)
	note the md5s don't change following filtering/normalization/etc... - only the original data
	input:
	scdb : from startdb()
	datamd5 : str
		from Experiment.datamd5
	mapmd5 : str
		from Experiment.mapmd5

	output:
	outlist:
		a list of ints of matching dataID indices (or empty if no match)
	"""
	Debug(1,'finddataid for data %s map %s' % (datamd5,mapmd5))
	outlist=[]
	if datamd5:
		db.cur.execute("SELECT DataID FROM Data WHERE Value = ?",[datamd5])
		allvals=db.cur.fetchall()
		for cres in allvals:
			outlist.append(cres[0])
	if mapmd5:
		db.cur.execute("SELECT DataID FROM Data WHERE Value = ?",[mapmd5])
		allvals=db.cur.fetchall()
		for cres in allvals:
			outlist.append(cres[0])
	outlist=list(set(outlist))
	Debug(2,"found %d matches to data" % len(outlist))
	return outlist


def getdatainfo(db,dataid):
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
	info=[]
	db.cur.execute("SELECT Type,Value FROM Data WHERE DataID = ?",[dataid])
	allvals=db.cur.fetchall()
	for cres in allvals:
		info.append((cres[0],cres[1],'%s:%s' % (cres[0],cres[1])))
	return info



def getstudyannotations(db,studyid):
	"""
	get the list of annotations for study studyid

	input:
	db : from dbstart()
	studyid : int
		The dataid of the study

	output:
	info: list of str
		the list of curations for this study (1 item per curation)
	"""
	info=[]
	db.cur.execute("SELECT Date,SubmitterName,Description FROM Curations WHERE DataID = ?",[studyid])

	allvals=db.cur.fetchall()
	for cres in allvals:
		info.append('%s:%s:%s' % (cres[0],cres[1],cres[2]))
	return info


def getseqcurationids(db,seqid):
	"""
	get the curation ids for the sequence with the id seqid

	input:
	db : from sbdtart()
	seqid : int
		the sequence identifier (from getseq)

	output:
	curids : list of int
		ids of all the curations for this sequence
	"""
	ids=[]
	db.cur.execute("SELECT CurationID FROM SeqCurations WHERE SeqID = ?",[seqid])
	allvals=db.cur.fetchall()
	for cres in allvals:
		ids.append(cres[0])
	Debug(2,'found %d curations' % len(ids))
	return ids


def getseqcurations(db,sequence):
	"""
	Get the manual curations for a sequence

	input:
	db : from sbdtart()
	sequence : str (ACGT)

	output:

	"""
	seqid=getseq(db,sequence,insert=False)
	if seqid==0:
		Debug(2,'Sequence not found')
		return
	curids=getseqcurationids(db,seqid)
	for cid in curids:
		cdat=select_column_and_value(db.cur,"SELECT * FROM Curations WHERE CurationID = ?",[cid])
		if cdat=={}:
			Debug(8,'no curation found for curationid %d' % cid)
			continue
		print(cdat)
		curation=db.cur.execute('SELECT * from CurationList WHERE CurationID = ?',[cid])
		for ccuration in curation.fetchall():
			print(ccuration)


def select_column_and_value(cur, sql, parameters=()):
	"""
	get a dict with field as key, value as values for an sql query
	"""
	execute = cur.execute(sql, parameters)
	fetch = execute.fetchone()

	if fetch is None:
		return {}

	return {k[0]: v for k, v in list(zip(execute.description, fetch))}
