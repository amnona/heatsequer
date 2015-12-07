#!/usr/bin/env python


"""
heatsequer experiment filtering module
"""

# amnonscript

__version__ = "0.9"


import heatsequer as hs

import numpy as np


def filterminreads(exp,minreads,logit=True):
	"""
	filter away all bacteria that contain less than minreads in all samples together (out of 10k/samples)
	input:
	exp : Experiment
	minreads : float
		the minimum number of reads total for all samples (and out of 10k/sample) for a bacteria to be kept
	logit : bool
		True to add to command log, False to not (if called from another heatsequer function)
	output:
	newexp - the filtered experiment
	"""
	params=locals()

	numreads=np.sum(exp.data,axis=1)
	keep=np.where(numreads>=minreads)
	newexp=hs.reorderbacteria(exp,keep[0])
	if logit:
		newexp.filters.append('filter min reads %d' % minreads)
		hs.addcommand(newexp,"filterminreads",params=params,replaceparams={'exp':exp})
	hs.Debug(6,'%d Bacteria left' % len(newexp.sids))
	return newexp

def filterpresence(expdat,frac):
	"""
	filter away bacteria present in less than frac of the samples
	input:
	expdat : Experiment
	frac : float
		the minimal fraction of samples to appear in for a beacteria to be kept

	output:
	newexp : Experiment
		the filtered experiment
	"""
	params=locals()

	fracreads=(np.sum(expdat.data>0,axis=1)+0.0)/len(expdat.samples)
	keep=np.where(fracreads>=frac)
	newexp=hs.reorderbacteria(expdat,keep[0])
	newexp.filters.append('filter presence %f' % frac)
	hs.addcommand(newexp,"filterpresence",params=params,replaceparams={'expdat':expdat})
	hs.Debug(6,'%d Bacteria left' % len(newexp.sids))
	return newexp


def filtermean(expdat,meanval):
	"""
	filter keeping bacteria with a mean frequency >= meanval
	input:
	expdat : Experiment
		the experiment
	meanval : float
		the minimum mean reads of per sample (and out of 10k/sample) for a bacteria to be kept
	output:
	newexp : Experiment
		the filtered experiment
	"""
	params=locals()

	meanreads=np.mean(expdat.data,axis=1)
	keep=np.where(meanreads>=meanval)
	newexp=hs.reorderbacteria(expdat,keep[0])
	newexp.filters.append('filter mean reads %f' % meanval)
	hs.addcommand(newexp,"filtermean",params=params,replaceparams={'expdat':expdat})
	hs.Debug(6,'%d Bacteria left' % len(newexp.sids))
	return newexp


def filterorigreads(expdat,minreads,inplace=False):
	"""
	filter away all samples that contained originally (before normalization) less than minreads
	input:
	expdat : Experiment
	minreads : float
		the minimum number of reads of the sample in the biom table to filter if less (usually int)
	inplace : Bool
		True to replace current experiment, False (default) to create a new one
	output:
	newexp : Experiment
		the filtered experiment
	"""
	params=locals()

	numreads=np.array(expdat.origreads)
	keep=np.where(numreads>=minreads)
	newexp=hs.reordersamples(expdat,keep[0])
	hs.Debug(6,'%d Samples left' % len(newexp.samples))
	hs.addcommand(newexp,"filterorigreads",params=params,replaceparams={'expdat':expdat})
	return newexp


def filtersamples(expdat,field,filtval,exact=True,exclude=False,numexpression=False):
	"""
	filter samples in experiment according to value in field
	input:
	exp : Experiment
	field : string
		name of the field to filter by
	filtval : string or list of strings
		the string to filter (if a list of strings, filter if any in the list)
	exact : bool
		True for exact match, False for substring
	exclude : bool
		False to keep only matching samples, True to exclude matching samples
	numexpression : bool
		True if val is a python expression, False if just a value. For an expression assume value is the beggining of the line (i.e. '<=5')
	"""
	params=locals()
	if not isinstance(filtval,list):
		filtval=[filtval]

	keep=[]
	for cidx,csamp in enumerate(expdat.samples):
		keepit=False
		for filt in filtval:
			if numexpression:
				cval=expdat.smap[csamp][field]
				if eval(cval+filt):
					keepit=True
			elif exact:
				if expdat.smap[csamp][field]==filt:
					keepit=True
			else:
				if filt in expdat.smap[csamp][field]:
					keepit=True
			# if exclude reverse the decision
		if exclude:
			keepit=not keepit
		if keepit:
			keep.append(cidx)
	newexp=hs.reordersamples(expdat,keep)
	fstr="filter data %s in %s " % (filt,field)
	if exact:
		fstr=fstr+"(exact)"
	else:
		fstr=fstr+"(substr)"
	if exclude:
		fstr+=" (exclude)"
	newexp.filters.append(fstr)
	hs.addcommand(newexp,"filtersamples",params=params,replaceparams={'expdat':expdat})
	hs.Debug(6,'%d Samples left' % len(newexp.samples))
	return newexp


def filterid(expdat,sids,exclude=False):
	"""
	filter bacteria keeping only ones in sids
	input:
	expdat : Experiment
	sids : list of integers
		the list of (hashed) sequence ids
	exclude : bool
		False to keep these bacteria, True to filter away
	output:
	newexp : Experiment
		the filtered experiment
	"""
	params=locals()

	if not type(sids) is list:
		sids=[sids]
	keep=[]
	hs.Debug(1,'filter ids',sids)
	for cid in sids:
		for idx,tid in enumerate(expdat.sids):
			if tid==cid:
				keep.append(idx)
	if exclude:
		keep=set(range(len(expdat.sids))).difference(keep)
	keep=list(set(keep))
	hs.Debug(1,'keep pos',keep)
	newexp=hs.reorderbacteria(expdat,keep)
	if exclude:
		newexp.filters.append('Filter %d ids (exclude)' % len(sids))
	else:
		newexp.filters.append('Filter %d ids' % len(sids))
	hs.addcommand(newexp,"filterid",params=params,replaceparams={'expdat':expdat})
	return newexp

def filtertaxonomy(exp,tax,exact=False,exclude=False):
	"""
	filter bacteria matching a given taxonomy name
	input:
	exp : Experiment
	tax : string
		the taxonomy name to filter by
	exact : bool
		True for exact matches to tax string, false for substring
	exclude : bool
		True to throw away matching taxonomy, False to keep matching
	"""
	params=locals()

	match=[]
	for cidx,ctax in enumerate(exp.tax):
		keep=False
		if exact:
			if ctax==tax:
				keep=True
		else:
			if tax in ctax:
				keep=True
		if exclude:
			keep=not keep
		if keep:
			match.append(cidx)
	newexp=hs.reorderbacteria(exp,match)
	filt='filter taxonomy '
	if exact:
		filt+='exact match '
	if exclude:
		filt+='exclude '
	filt+=tax
	newexp.filters.append(filt)
	hs.Debug(6,'%d bacteria left' % len(newexp.sids))
	hs.addcommand(newexp,"filtertaxonomy",params=params,replaceparams={'exp':exp})
	return newexp


def clearexp(expdat):
	"""
	clear experiment from missing samples and bacteria
	remove samples with <1 reads and bacteria with total <1 reads
	input:
	expdat : Experiment
	output:
	newexp : Experiment
		the new experiment without <1 reads samples or bacteria
	"""
	params=locals()

	newexp=filterorigreads(expdat,1)
	newexp=filterminreads(expdat,1)
	newexp.filters.append('clear nonpresent bacteria and samples')
	hs.Debug(6,'%d bacteria left' % len(newexp.sids))
	hs.addcommand(newexp,"clearexp",params=params,replaceparams={'expdat':expdat})
	return newexp


def filterseqs(expdat,seqs,exclude=False,subseq=False,logit=True):
	"""
	filter sequences from the list seqs (keeping sequences appearing in the list)
	input:
	expdat : Experiment
	seqs : string
		a list of (ACGT) sequences to keep
	exclude : bool
		True to filter away instead of keep, False to keep
	subseq : bool
		if true, the sequences can be subsequence (slower). False - look only for exact match.
	logit : bool
		True to add to command log, false to not log it (if called from other heatsequer function)

	output:
	newexp : Experiment
		the filtered experiment
	"""
	params=locals()

	keeplist=[]
	for cseq in seqs:
		if subseq:
			for idx,coseq in enumerate(expdat.seqs):
				if len(cseq)>len(coseq):
					if coseq in cseq:
						keeplist.append(idx)
						break
				else:
					if cseq in coseq:
						keeplist.append(idx)
						break
		else:
			if cseq in expdat.seqdict:
				keeplist.append(expdat.seqdict[cseq])
			else:
				hs.Debug(7,'sequence not in experiment',cseq)
	if exclude:
		keeplist=list(set(range(len(expdat.seqs))).difference(keeplist))
	newexp=hs.reorderbacteria(expdat,keeplist)
	if logit:
		newexp.filters.append('filter sequences')
		hs.addcommand(newexp,"filterseqs",params=params,replaceparams={'expdat':expdat})
	return newexp


def filterseqsfromexp(expdat,filtexp,exclude=False,subseq=False,logit=True):
	"""
	filter sequences from expdat using the sequences in filtexp experiment (keeping sequences appearing in the other experiment)
	input:
	expdat : Experiment
	filtexp : string
		an experiment from which to decide which sequences to filter
	exclude : bool
		True to filter away instead of keep, False to keep
	subseq : bool
		if true, the sequences can be subsequence (slower). False - look only for exact match.
	logit : bool
		True to add to command log, false to not log it (if called from other heatsequer function)

	output:
	newexp : Experiment
		the filtered experiment
	"""
	params=locals()
	newexp=hs.filterseqs(expdat,filtexp.seqs,exclude=exclude,subseq=subseq,logit=False)
	if logit:
		newexp.filters.append('filter sequences')
		hs.addcommand(newexp,"filterseqsfromexp",params=params,replaceparams={'expdat':expdat,'filtexp':filtexp})
	return newexp


def filterknownbact(expdat,cdb,exclude=False):
	"""
	filter keeping only bacteria which we know about in cooldb
	input:
	expdat : Experiment
	cdb : cooldb
		the manual annotation database (fromn cooldb.loaddb)
	exclude : bool
		True to throw away known bacteria, False to keep only them
	output:
	newexp : Experiment
		the filtered experiment
	"""
	params=locals()

	known=[]
	for idx,cseq in enumerate(expdat.seqs):
		if len(hs.cooldb.getseqinfo(cdb,cseq))>0:
			known.append(idx)
	hs.Debug(2,'Found %d sequences known in cooldb' % len(known))
	if exclude:
		known=set(range(len(expdat.seqs))).difference(known)
	newexp=hs.reorderbacteria(expdat,known)
	if not exclude:
		newexp.filters.append('filter cooldb known bacteria')
	else:
		newexp.filters.append('filter exclude cooldb known bacteria')
	hs.Debug(6,'%d bacteria left' % len(newexp.sids))
	newexp.filters.append('keep only sequences from cooldb')
	hs.addcommand(newexp,"filterknownbact",params=params,replaceparams={'expdat':expdat})
	return newexp


def filtersamplesfromfile(expdat,filename,field='#SampleID',exclude=False):
	"""
	filter samples based on a list in a file (one line per sample)
	input:
	expdat
	filename - name of the text file (1 line per sample)
	field - field that the file contains
	exclude - true to throw away instead of keep
	output:
	newexp - the filtered experiment
	"""
	params=locals()

	fl=open(filename,'r')
	vals=[]
	for cline in fl:
		vals.append(cline)
	fl.close()

	keep=[]
	for cidx,csamp in enumerate(expdat.samples):
		keepit=False
		if expdat.smap[csamp][field] in vals:
			keepit=True
		# if exclude reverse the decision
		if exclude:
			keepit=not keepit
		if keepit:
			keep.append(cidx)
	newexp=hs.reordersamples(expdat,keep)

	fstr="filter data from file %s in %s " % (filename,field)
	if exclude:
		fstr+=" (exclude)"
	newexp.filters.append(fstr)
	hs.addcommand(newexp,"filtersamplesfromfile",params=params,replaceparams={'expdat':expdat})
	hs.Debug(6,'%d Samples left' % len(newexp.samples))
	return newexp


def filterfasta(expdat,filename,exclude=False,subseq=False):
	"""
	Filter bacteria based on a fasta file (keeping only sequences in the fasta file)
	input:
	expdat
	filename - the fasta file name
	exclude - False to keep only seqs from the file, True to throw away seqs from the file
	subseq - match subsequences (slower)

	output:
	newexp - the filtered experiment
	"""
	params=locals()

	seqs,headers=hs.readfastaseqs(filename)
	newexp=hs.filterseqs(expdat,seqs,exclude=exclude,subseq=subseq)
	filt='Filter sequences from file '+filename
	if exclude:
		filt+=' (Exclude)'
	if subseq:
		filt+=' (subseq)'
	newexp.filters.append(filt)
	hs.addcommand(newexp,"filterfasta",params=params,replaceparams={'expdat':expdat})
	return newexp


def filterbacteriafromfile(expdat,filename,exclude=False,subseq=False):
	"""
	filter bacteria from an experiment based on a file with sequences (one per line)
	input:
	expdat
	filename - name of the sequence file (1 per line)
	exclude - remove bacteria from the file instead of keeping them
	subseq - the sequences in the file can be subsequences of the experiment sequences (different lengths). but slower.

	output:
	newexp - the filtered experiment
	"""
	params=locals()

	fl=open(filename,'rU')
	seqs=[]
	for cline in fl:
		seqs.append(cline.strip())
	newexp=hs.filterseqs(expdat,seqs,exclude=exclude,subseq=False)
	filt='Filter sequences from file '+filename
	if exclude:
		filt+=' (Exclude)'
	if subseq:
		filt+=' (subseq)'
	newexp.filters.append(filt)
	hs.addcommand(newexp,"filterbacteriafromfile",params=params,replaceparams={'expdat':expdat})
	return newexp



def filterandnormalize(expdat,seqs,exclude=False,subseq=False,numreads=10000):
	"""
	filter away sequences in seqs
	and then renormalize the data and recalculate origreads per sample
	input:
	expdat
	seqs - a list of sequences (ACGT) to remove
	exclude - False to remove seqs, True to only keep seqs
	subseq - False to look for exact match only, True to look for subsequence match (slower)

	output:
	newexp - the experiment without seqs and renormalied to 10k reads/sample
	"""
	params=locals()

	newexp=hs.filterseqs(expdat,seqs,exclude=not(exclude),subseq=subseq)
	newexp=hs.normalizereads(newexp,fixorig=True,numreads=numreads)
	newexp.filters.append("filter sequences and normalize to numreads %d" % numreads)
	hs.addcommand(newexp,"filterandnormalize",params=params,replaceparams={'expdat':expdat})
	return newexp


def filterannotations(expdat,annotation,cdb,exclude=False):
	"""
	filter keeping only samples which have annotation in their cooldb description
	input:
	expdat
	annotation - substring of the annotation (case insensitive)
	cdb - the database of cool sequences (from cooldb.load())
	exclude - False to keep matching bacteria, True to remove matching bacteria

	output:
	newexp - the filtered experiment
	"""
	params=locals()

	keeplist=[]
	for idx,cseq in enumerate(expdat.seqs):
		keep=False
		info=hs.cooldb.getseqinfo(cdb,cseq)
		for cinfo in info:
			if annotation.lower() in str(cinfo).lower():
				keep=True
		if exclude:
			keep = not keep
		if keep:
			keeplist.append(idx)
	newexp=hs.reorderbacteria(expdat,keeplist)
	newexp.filters.append('Filter annotations %s' % annotation)
	hs.addcommand(newexp,"filterannotations",params=params,replaceparams={'expdat':expdat})
	hs.Debug(6,'%d bacteria found' % len(keeplist))
	return newexp
