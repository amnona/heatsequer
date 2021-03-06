#!/usr/bin/env python


"""
heatsequer cooldb manual curation database functions
"""

# amnonscript

__version__ = "0.9"

from ..utils.amnonutils import mlhash,Debug,listtodict,fdr
from ..utils.sequence import seqdist

import numpy as np
import biom
import csv
from scipy import stats
import datetime
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

		# True if we also loaded the greengenes ids from getggids
		self.biomggloaded=False

		# number of times each description appeats
		self.numtimes={}

def loaddb(dbname='./db/coolseqs.txt'):
#def loaddb(dbname='./db/coolseqs.gg97-135.txt'):
	db=cooldb()
	db.dbfile=dbname
	Debug(0,'Loadding coolseq database',db.dbfile)
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
	Debug(1,'Loaded coolseq database','num sequences',idx)
	return db


def savecoolseqs(expdat,cdb,seqs,description):
	"""
	save sequences to the cooldb database
	input:
	expdat : Experiment
	seqs : list of sequence strings
		a list of sequence strings (ACGT) to save
	description : string
		the description for these sequences to be put in coolseqs database
	"""

	numwritten=0
	for cseq in seqs:
		if cseq in expdat.seqdict:
			seqpos=expdat.seqdict[cseq]
			saveseq(db=cdb,seq=cseq,ggid=expdat.sids[seqpos],taxonomy=expdat.tax[seqpos],filename=expdat.tablefilename,expdescription=expdat.studyname,description=description)
			numwritten+=1
		else:
			Debug(9,'sequence:\n%s\nnot in experiment!!!' % cseq)
	Debug(6,'Wrote %d sequences to cooldb' % numwritten)


def getseqinfo(db,seq):
	'''
	get all known information for a given sequence (minimum length 89bp)
	input:
	db - from loaddb()
	seq - the sequence to search for
	output:
	info - a list of known info about the sequence
	'''

	# if seq is numeric, it is a greengenesid - need the database after getggids
	if seq.isdigit():
		Debug(0,'Looking for info for ggid %s' % seq)
		if not db.biomggloaded:
			Debug(9,'GreenGenes IDs were not loaded for cooldb!')
		ggid=int(seq)
		info=[]
		for cidx in range(len(db.dat)):
			if 'biomggid' in db.dat[cidx]:
				if int(db.dat[cidx]['biomggid'])==ggid:
					info.append(db.dat[cidx]['bacteria_description'])
		return info
	# not a ggid so look for the sequence
	Debug(0,'Looking for info for sequence',seq)
	info=[]
	sseq=seq[0:89]
	if not sseq in db.shortseqdict:
		Debug(0,'sequence not found in short hash')
		return info
	for cpos in db.shortseqdict[sseq]:
		fulldbseq=db.dat[cpos]['sequence']
		mlen=min(len(fulldbseq),len(seq))
		seq=seq[0:mlen]
		fulldbseq=fulldbseq[0:mlen]
		if not seq==fulldbseq:
			Debug(0,'sequence matches short hash but not identical',cpos)
			continue
		info.append(db.dat[cpos]['bacteria_description'])
		Debug(0,'found in position',cpos)
	return info


def getcloseseqinfo(db,oseq,maxdist):
	"""
	get database info from all sequences in database close up to maxdist mismatches from oseq
	input:
	db - the initiated database file
	oseq : sequence (ACGT string)
		the sequence to find neoghbors of in the database
	maxdist : integer
		find only sequences with hamming<=maxdist
	output:
	info : list of strings
		database entries of all close sequences
	"""
	info=[]
	for cdbent in db.dat:
		cseq=cdbent['sequence']
		cdist=seqdist(oseq,cseq)
		if cdist<=maxdist:
			info.append('(%d) %s' % (cdist,cdbent['bacteria_description']))
			Debug(1,'Found similar sequence dist %d' % cdist)
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
		ggid=mlhash(seq)
	if not expdescription:
		expdescription=description
	seq89=seq[:89]
	ggid89=mlhash(seq89)
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


def exportfasta(db,filename):
	"""
	export the cooldb to fasta format
	used for conversion to greengenes
	the header is the actual sequence!
	input:
	db - the cooldb
	filename - name of the output fasta file
	"""

	fl=open(filename,'w')
	for cent in db.dat:
		fl.write('>%s\n' % cent['sequence'])
		fl.write('%s\n' % cent['sequence'])
	fl.close()


def getggids(db,biomname):
	"""
	load the closed reference biom table for the sequences in exportfasta
	and convert the database to greengenesIDs.
	input:
	db - the cooldb used for the exportfasta
	biomname - the biom table file obtained by running pick_closed_reference on the exportfasta file

	output:
	db - the new cooldb with the 'biomggid' field added
	"""

	# load the biom table
	table = biom.load_table(biomname)

	# each sample is a sequence from cooldb
	tseqs = table.ids(axis='sample')
	tggids= table.ids(axis='observation')

	for cseq in tseqs:
		# find the position of the sequence in cooldb
		spos=-1
		sseq=cseq[0:89]
		if not sseq in db.shortseqdict:
			Debug(8,'sequence %s not found in database short hash' % cseq)
			return db
		for cpos in db.shortseqdict[sseq]:
			fulldbseq=db.dat[cpos]['sequence']
			if cseq==fulldbseq:
				spos=cpos
				break
		if spos==-1:
			Debug(8,'sequence %s not in database' % cseq)
			return db
		reads=table.data(cseq,axis='sample')
		pnz=np.where(reads>0)
		if not len(pnz[0]) == 1:
			Debug(9,'same sequence assigned to multiple ggids!',cseq)
			return db

		db.dat[spos]['biomggid']=tggids[pnz[0]][0]

	db.biomggloaded=True
	return db



def savedb(db,filename):
	"""
	save the coolseqs database
	for use after getggids to add the biomggid field to the database
	input:
	db - the coolseq database
	filename - name of the output file to write
	"""

	fl=open(filename,'w')
	fields=db.dat[0].keys()
	# write the header (tsv)
	for cfield in fields:
		fl.write('%s\t' % cfield)
	fl.write('\n')
	for cdat in db.dat:
		for cfield in fields:
			if cfield in cdat:
				fl.write('%s\t' % cdat[cfield])
			else:
				fl.write('%s\t' % '0')
		fl.write('\n')
	fl.close()
	Debug(4,'Saved cooldb',filename)


def preparealldbvals(db):
	"""
	calculate the amount of times each description appears
	inputL
	db

	output:
	db - with the added db.numtimes dict for each description
	"""
	db.numtimes={}
	for cdat in db.dat:
		if cdat['bacteria_description'] in db.numtimes:
			db.numtimes+=1
		else:
			db.numtimes=1
	return db


def testenrichment(db,allseqs=[],group=[],maxfval=0.1,freqs=[]):
	"""
	test for database entries enriched in sequences from group compared to allseqs (all the sequences) or all the database
	input:
	db - the coolseq database
	allseqs - a list of all the sequences or [] to compare to all the database
	group - a list of subgroup of sequences to test vs allseqs
	maxfval - the maximal fdr value
	freqs - the mean frequency of each sequence, or [] to not use frequency weighing

	output:
	newplist - a sorted list of dict for annotaions which are below fdr ('description','pval','observed','expected')
	"""

	if len(allseqs)>0:
		alllen=len(allseqs)
		Debug(6,'comparing list of %d seqs to %d sequences' % (alllen,len(group)))
	else:
		alllen=len(db.dat)
		Debug(6,'comparing list of %d seqs to all %d sequences' % (alllen,len(group)))
	glen=len(group)

	if len(allseqs)>0:
		allshort=[]
		for cseq in allseqs:
			allshort.append(cseq[:89])
		asdict=listtodict(allshort)

	groupshort=[]
	for cseq in group:
		groupshort.append(cseq[:89])
	gsdict=listtodict(groupshort)

	dbdesc=[]
	for cdat in db.dat:
		dbdesc.append(cdat['bacteria_description'])

	descdict=listtodict(dbdesc)

	allp=[]
	pv=[]
	for k,v in descdict.items():
		allmatch=0
		groupmatch=0
		usedseq={}
		usedseq2={}
		usedseq3={}
		for cdatidx in v:
			cdat=db.dat[cdatidx]
			cseq=cdat['sequence']
			if cseq in usedseq:
				continue
			usedseq[cseq]=True
			cseqs=cseq[:89]
			if len(allseqs)>0:
				if cseqs not in asdict:
					continue
				for opos in asdict[cseqs]:
					oseq=allseqs[opos]
					mlen=min(len(oseq),len(cseq))
					if not oseq[:mlen]==cseq[:mlen]:
						continue
					if oseq in usedseq2:
						continue
					allmatch+=1
					usedseq2[oseq]=True
			else:
				allmatch+=1
			if cseqs not in gsdict:
				continue
			for opos in gsdict[cseqs]:
				oseq=group[opos]
				mlen=min(len(oseq),len(cseq))
				if not oseq[:mlen]==cseq[:mlen]:
					continue
				if oseq in usedseq3:
					continue
				usedseq3[oseq]=True
				if len(freqs)>0:
					groupmatch+=freqs[opos]
				else:
					groupmatch+=1
		if allmatch==0:
			continue
		pnull=float(allmatch)/alllen
		p1=1-stats.binom.cdf(groupmatch-1,glen,pnull)
		if len(allseqs)>0:
			p2=stats.binom.cdf(glen-groupmatch,glen,1-pnull)
		else:
			p2=stats.binom.cdf(glen-groupmatch,glen,1-pnull)
#		p=min(p1,p2)
#		p=p2
		p=p1
		allp.append(p)
		cpv={}
		cpv['pval']=p
		cpv['observed']=groupmatch
		cpv['expected']=pnull*glen
		cpv['description']=k
		pv.append(cpv)

	fval=fdr(allp)
	keep=np.where(np.array(fval)<=maxfval)
	plist=[]
	rat=[]
	for cidx in keep[0]:
		plist.append(pv[cidx])
		rat.append(np.abs(float(pv[cidx]['observed']-pv[cidx]['expected']))/np.mean([pv[cidx]['observed'],pv[cidx]['expected']]))
	si=np.argsort(rat)
	si=si[::-1]
	newplist=[]
	for idx,crat in enumerate(rat):
		print(plist[si[idx]])
		newplist.append(plist[si[idx]])
	return(newplist)


