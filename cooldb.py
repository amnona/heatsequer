#!/usr/bin/env python

"""
amnonscript
cooldb.py

analyzing the cool sequences database (manual curation)
located at ~/Databases/coolseqs/coolseqs.txt
"""

__version__ = "0.1"

import amnonutils as au

import sys
import numpy as np
import matplotlib.pyplot as plt
import biom
from matplotlib.pyplot import *
import csv
from scipy import cluster
from scipy import spatial
from sklearn.preprocessing import scale
import datetime
import copy
import sqlite3
from collections import defaultdict
# for debugging - use XXX()
from pdb import set_trace as XXX

class cooldb:
	def __init__(self):
		# name of the database file
		self.dbfile=''

		# the data (one line per sequence)
		self.dat=[]

		# 89 bp hash of sequences pointing to location in dat
		self.shortseqdict={}


def loaddb(dbname='./db/coolseqs.txt'):
	db=cooldb()
	db.dbfile=dbname
	au.Debug(0,'Loadding coolseq database',db.dbfile)
	dbf = open(db.dbfile, 'rU')
	reader = csv.DictReader(dbf, delimiter='\t')
	idx=0
	for cline in reader:
		db.dat.append(cline)
		cseq=cline['sequence']
		cseq=cseq[:89]
		if cseq in db.shortseqdict:
			db.shortseqdict[cseq].append(idx)
		else:
			db.shortseqdict[cseq]=[idx]
		idx+=1
	dbf.close()
	au.Debug(1,'Loaded coolseq database','num sequences',idx)
	return db


def getseqinfo(db,seq):
	'''
	get all known information for a given sequence (minimum length 89bp)
	input:
	db - from loaddb()
	seq - the sequence to search for
	output:
	info - a list of known info about the sequence
	'''

	info=[]
	sseq=seq[0:89]
	if not sseq in db.shortseqdict:
		au.Debug(0,'sequence not found in short hash')
		return info
	for cpos in db.shortseqdict[sseq]:
		fulldbseq=db.dat[cpos]['sequence']
		mlen=min(len(fulldbseq),len(seq))
		seq=seq[0:mlen]
		fulldbseq=fulldbseq[0:mlen]
		if not seq==fulldbseq:
			au.Debug(0,'sequence matches short hash but not identical',cpos)
			continue
		info.append(db.dat[cpos]['bacteria_description'])
		au.Debug(0,'found in position',cpos)
	return info

def saveseq(db,seq,taxonomy,filename,description,ggid=False,expdescription=False):
	""""
	Save a sequence into the database
	input:
	db - the initiated database file
	seq - the sequence (ACGT)
	ggid - the hashed id or False to calculate here
	taxonomy - the taxonomy for the sequence
	filename - the name of the biom table file
	expdescription - description of the experiment (similar to filename if False)
	description - text about the bacteria

	output:
	db - the modified database list (for immediate use - no need to reload)!
	"""

	if not ggid:
		ggid=au.mlhash(seq)
	if not expdescription:
		expdescription=description
	seq89=seq[:89]
	ggid89=au.mlhash(seq89)
	ctime=datetime.date.today()
	fl=open(db.dbfile,'a')
	fl.write('%s\t' % str(ggid))
	fl.write('%d\t' % len(seq))
	fl.write('%s\t' % seq)
	fl.write('%s\t' % str(ggid89))
	fl.write('%s\t' % seq89)
	fl.write('%s\t' % taxonomy)
	fl.write('%s\t' % filename)
	fl.write('%s\t' % expdescription)
	fl.write('%s\t' % ctime.isoformat())
	fl.write('%s\n' % description)
	fl.close()

	cdat={}
	cdat['GGID']=str(ggid)
	cdat['length']=str(len(seq))
	cdat['sequence']=seq
	cdat['ggid89']=str(ggid89)
	cdat['seq89']=seq89
	cdat['taxonomy']=taxonomy
	cdat['filename']=filename
	cdat['exp_description']=expdescription
	cdat['add_date']=ctime.isoformat()
	cdat['bacteria_description']=description

	idx=len(db.dat)
	db.dat.append(cdat)
	if seq89 in db.shortseqdict:
		db.shortseqdict[seq89].append(idx)
	else:
		db.shortseqdict[seq89]=[idx]

	return db
