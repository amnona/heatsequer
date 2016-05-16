#!/usr/bin/env python


"""
heatsequer experiment filtering module
"""

# amnonscript

__version__ = "0.9"


import heatsequer as hs

from scipy import stats
import numpy as np
from pdb import set_trace as XXX


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
		the minimum mean fraction of reads per sample (NOT out of 10k/sample) for a bacteria to be kept
	output:
	newexp : Experiment
		the filtered experiment
	"""
	params=locals()

	meanreads=np.mean(expdat.data,axis=1)
	meantotreads=np.mean(np.sum(expdat.data,0))
	keep=np.where(meanreads>=meanval*meantotreads)
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


def filtersamples(expdat,field,filtval,exact=True,exclude=False,numexpression=False,shownumoutput=True):
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
	shownumoutput : bool
		True (default) to show number of samples remaining, False to not show
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
				if len(cval)==0:
					continue
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
	if shownumoutput:
		hs.Debug(6,'%d Samples left' % len(newexp.samples))
	else:
		hs.Debug(1,'%d Samples left' % len(newexp.samples))
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


def filterknownbact(expdat,cdb=None,exclude=False):
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

	if cdb is None:
		cdb=hs.cdb
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
	field - field for the experiment that the file contains
	exclude - true to throw away instead of keep
	output:
	newexp - the filtered experiment
	"""
	params=locals()

	fl=open(filename,'rU')
	vals=[]
	for cline in fl:
		cline=cline.strip()
		vals.append(cline)
	fl.close()

	keep=[]
	for cidx,csamp in enumerate(expdat.samples):
		keepit=False
		if expdat.smap[csamp][field] in vals:
			hs.Debug(3,'found file value %s in sample %s' % (expdat.smap[csamp][field],csamp))
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


def filterannotations(expdat,annotation,cdb=None,exclude=False):
	"""
	filter keeping only samples which have annotation in their cooldb description
	input:
	expdat
	annotation - substring of the annotation (case insensitive)
	cdb - the database of cool sequences (from cooldb.load()) or None (default) to use the heatsequer loaded cdb
	exclude - False to keep matching bacteria, True to remove matching bacteria

	output:
	newexp - the filtered experiment
	"""
	params=locals()

	if cdb is None:
		cdb=hs.cdb
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


def filtersimilarsamples(expdat,field,method='mean'):
	"""
	join similar samples into one sample (i.e. to remove samples of same individual)
	input:
	expdat : Experiment
	field : string
		Name of the field containing the values (for which similar values will be joined)
	method : string
		What to do with samples with similar value. options:
		'mean' - replace with a sample containing the mean of the samples
		'median'- replace with a sample containing the median of the samples
		'random' - replace with a single random sample out of these samples
		'sum' - replace with sum of original reads in all samples, renormalized after to 10k and orignumreads updated
	output:
	newexp : Experiment
		like the input experiment but only one sample per unique value in field
	"""
	params=locals()

	newexp=hs.copyexp(expdat)
	if method=='sum':
		newexp=hs.toorigreads(newexp)
	uvals=hs.getfieldvals(expdat,field,ounique=True)
	keep=[]
	for cval in uvals:
		cpos=hs.findsamples(expdat,field,cval)
		if len(cpos)==1:
			keep.append(cpos[0])
			continue
		if method=='random':
			keep.append(cpos[np.random.randint(len(cpos))])
			continue
		# set the mapping file values
		cmap=expdat.smap[expdat.samples[cpos[0]]]
		for ccpos in cpos[1:]:
			for cfield in cmap.keys():
				if cmap[cfield]!=expdat.smap[expdat.samples[ccpos]][cfield]:
					cmap[cfield]='NA'
		if method=='mean':
			cval=np.mean(expdat.data[:,cpos],axis=1)
			newexp.data[:,cpos[0]]=cval
			keep.append(cpos[0])
		elif method=='median':
			cval=np.median(expdat.data[:,cpos],axis=1)
			newexp.data[:,cpos[0]]=cval
			keep.append(cpos[0])
		elif method=='sum':
			cval=np.sum(newexp.data[:,cpos],axis=1)
			newexp.data[:,cpos[0]]=cval
			newexp.origreads[cpos[0]]=np.sum(hs.reorder(expdat.origreads,cpos))
			keep.append(cpos[0])
		else:
			hs.Debug(9,'method %s not supported' % method)
			return False
		newexp.smap[expdat.samples[cpos[0]]]=cmap
	newexp=hs.reordersamples(newexp,keep)
	if method=='sum':
		newexp=hs.normalizereads(newexp)
	newexp.filters.append('Filter similar samples field %s method %s' % (field,method))
	hs.addcommand(newexp,"filtersimilarsamples",params=params,replaceparams={'expdat':expdat})
	hs.Debug(6,'%d samples before filtering, %d after' % (len(expdat.samples),len(newexp.samples)))
	return newexp


def filterwave(expdat,field=False,numeric=True,minfold=2,minlen=3,step=1,direction='up',posloc='start'):
	"""
	filter bacteria, keeping only ones that show a consecutive region of samples with higher/lower mean than other samples
	Done by scanning all windowlen/startpos options for each bacteria
	input:
	expdat : Experiment
	field : string
		The field to sort by or False to skip sorting
	numeric : bool
		For the sorting according to field (does not matter if field is False)
	minfold : float
		The minimal fold change for the window compared to the rest in order to keep
	step : int
		The skip between tested windows (to make it faster use a larger skip)
	minlen : int
		The minimal window len for over/under expression testing
	direction : string
		'both' - test both over and under expression in the window
		'up' - only overexpressed
		'down' - only underexpressed
	posloc : string
		The position to measure the beginning ('maxstart') or middle ('maxmid') of maximal wave
		or 'gstart' to use beginning of first window with >=minfold change

	output:
	newexp : Experiment
		The filtered experiment, sorted according to window start samples position
	"""
	params=locals()

	# sort if needed
	if field:
		newexp=hs.sortsamples(expdat,field,numeric=numeric)
	else:
		newexp=hs.copyexp(expdat)

	dat=newexp.data
	dat[dat<1]=1
	dat=np.log2(dat)
	numsamples=len(newexp.samples)
	numbact=len(newexp.seqs)
	maxdiff=np.zeros([numbact])
	maxpos=np.zeros([numbact])-1
	maxlen=np.zeros([numbact])
	for startpos in range(numsamples-minlen):
		for cwin in np.arange(minlen,numsamples-startpos,step):
			meanin=np.mean(dat[:,startpos:startpos+cwin],axis=1)
			nowin=[]
			if startpos>0:
				nonwin=np.arange(startpos-1)
			if startpos<numsamples:
				nowin=np.hstack([nowin,np.arange(startpos,numsamples-1)])
			nowin=nowin.astype(int)
			meanout=np.mean(dat[:,nowin],axis=1)
			cdiff=meanin-meanout
			if direction=='both':
				cdiff=np.abs(cdiff)
			elif direction=='down':
				cdiff=-cdiff
			if posloc=='gstart':
				usepos=np.logical_and(cdiff>=minfold,maxpos==-1)
				maxpos[usepos]=startpos
			elif posloc=='start':
				maxpos[cdiff>maxdiff]=startpos
			elif posloc=='mid':
				maxpos[cdiff>maxdiff]=startpos+int(cwin/2)
			else:
				hs.Debug('posloc nut supported %s' % posloc)
				return False
			maxlen[cdiff>maxdiff]=cwin
			maxdiff=np.maximum(maxdiff,cdiff)

	keep=np.where(maxdiff>=minfold)[0]
	keeppos=maxpos[keep]
	si=np.argsort(keeppos)
	keep=keep[si]
	for ci in keep:
		hs.Debug(6,'bacteria %s startpos %d len %d diff %f' % (newexp.tax[ci],maxpos[ci],maxlen[ci],maxdiff[ci]))
	newexp=hs.reorderbacteria(newexp,keep)
	newexp.filters.append('Filter wave field=%s minlen=%d' % (field,minlen))
	hs.addcommand(newexp,"filterwave",params=params,replaceparams={'expdat':expdat})
	hs.Debug(6,'%d samples before filtering, %d after' % (len(expdat.samples),len(newexp.samples)))
	return newexp


def filtern(expdat):
	"""
	delete sequences containing "N" from experiment and renormalize
	input:
	expdat : Experiment
	output:
	newexp : Experiment
		experiment without sequences containing "N"
	"""
	params=locals()

	keeplist=[]
	for idx,cseq in enumerate(expdat.seqs):
		if "N" in cseq:
			continue
		if "n" in cseq:
			continue
		keeplist.append(idx)
	newexp=hs.reorderbacteria(expdat,keeplist)
	newexp=hs.normalizereads(newexp)
	newexp.filters.append('Filter sequences containing N')
	hs.addcommand(newexp,"filtern",params=params,replaceparams={'expdat':expdat})
	hs.Debug(6,'%d sequences before filtering, %d after' % (len(expdat.seqs),len(newexp.seqs)))
	return newexp


def cleantaxonomy(expdat,mitochondria=True,chloroplast=True,bacteria=True,unknown=True,exclude=True):
	"""
	remove common non-16s sequences from the experiment and renormalize

	input:
	expdat : Experiment
	mitochondria : bool
		remove mitochondrial sequences
	chloroplast : bool
		remove chloroplast sequences
	bacteria : bool
		remove sequences only identified as "Bacteria" (no finer identification)
	unknown : bool
		remove unknown sequences
	exclude : bool
		True (default) to remove these sequecnes, False to keep them and throw other

	output:
	newexp : Experiment
		the renormalized experiment without these bacteria
	"""
	params=locals()

	newexp=hs.copyexp(expdat)
	if mitochondria:
		if exclude:
			newexp=hs.filtertaxonomy(newexp,'mitochondria',exclude=True)
		else:
			ne1=hs.filtertaxonomy(newexp,'mitochondria',exclude=False)
	if chloroplast:
		if exclude:
			newexp=hs.filtertaxonomy(newexp,'Streptophyta',exclude=True)
			newexp=hs.filtertaxonomy(newexp,'Chloroplast',exclude=True)
		else:
			ne2=hs.filtertaxonomy(newexp,'Streptophyta',exclude=False)
			ne3=hs.filtertaxonomy(newexp,'Chloroplast',exclude=False)
	if unknown:
		if exclude:
			newexp=hs.filtertaxonomy(newexp,'Unknown',exclude=True)
			newexp=hs.filtertaxonomy(newexp,'Unclassified;',exclude=True,exact=True)
		else:
			ne4=hs.filtertaxonomy(newexp,'Unknown',exclude=False)
			ne5=hs.filtertaxonomy(newexp,'Unclassified;',exclude=False,exact=True)
	if bacteria:
		if exclude:
			newexp=hs.filtertaxonomy(newexp,'Bacteria;',exclude=True,exact=True)
		else:
			ne6=hs.filtertaxonomy(newexp,'Bacteria;',exclude=False,exact=True)
	if exclude:
		newexp=hs.normalizereads(newexp)
	else:
		allseqs=[]
		allseqs+=(ne1.seqs)
		allseqs+=(ne2.seqs)
		allseqs+=(ne3.seqs)
		allseqs+=(ne4.seqs)
		allseqs+=(ne5.seqs)
		allseqs+=(ne6.seqs)
		allseqs=list(set(allseqs))
		newexp=hs.filterseqs(newexp,allseqs)
	newexp.filters.append('Clean Taxonomy (remove mitochondria etc.)')
	hs.addcommand(newexp,"cleantaxonomy",params=params,replaceparams={'expdat':expdat})
	hs.Debug(6,'%d sequences before filtering, %d after' % (len(expdat.seqs),len(newexp.seqs)))
	return newexp



def filterfieldwave(expdat,field,val1,val2=False,mineffect=1,method='mean',uselog=True):
	"""
	find all sequences which show an effect size of at least mineffect between val1 and val2 samples in field
	no statistical significance testing is performed

	input:
	expdat : Experiment
	field : string
		name of field to use for group separation
	val1 : string
		value in field for group1
	val2 : string
		value in field for group2 or False for all the other samples except val1
	mineffect : float
		min difference between groups per OTU in order to keep
	method: string
		'ranksum'
	uselog : bool
		True to log transform the data

	output:
	newexp : Experiment
		only with sequences showing a mineffect difference
	"""
	params=locals()

	numseqs=len(expdat.seqs)
	numsamples=len(expdat.samples)
	dat=expdat.data
	if uselog:
		dat[dat<1]=1
		dat=np.log2(dat)
	if method=='ranksum':
		for idx in range(numseqs):
			dat[idx,:]=stats.rankdata(dat[idx,:])

	pos1=hs.findsamples(expdat,field,val1)
	if val2:
		pos2=hs.findsamples(expdat,field,val2)
	else:
		pos2=np.setdiff1d(np.arange(numsamples),pos1,assume_unique=True)

	outpos=[]
	odif=[]
	for idx in range(numseqs):
		cdif=np.mean(dat[idx,pos1])-np.mean(dat[idx,pos2])
		if abs(cdif)>=mineffect:
			outpos.append(idx)
			odif.append(cdif)

	si=np.argsort(odif)
	outpos=hs.reorder(outpos,si)
	newexp=hs.reorderbacteria(expdat,outpos)
	newexp.filters.append('filterfieldwave field %s val1 %s val2 %s' % (field,val1,val2))
	hs.addcommand(newexp,"filterfieldwave",params=params,replaceparams={'expdat':expdat})
	return newexp


def filterwinperid(expdat,idfield,field,val1,val2,mineffect=1):
	"""
	do filterfieldwave on each individual (based on idfield) and join the resulting bacteria
	"""
	params=locals()

	iseqs=[]
	uids=hs.getfieldvals(expdat,idfield,ounique=True)
	for cid in uids:
		cexp=hs.filtersamples(expdat,idfield,cid)
		texp=hs.filterfieldwave(cexp,field,val1,val2,mineffect=mineffect)
		iseqs+=texp.seqs
	iseqs=list(set(iseqs))
	newexp=hs.filterseqs(expdat,iseqs)
	return newexp


def filternans(expdat,minpresence):
	"""
	filter an experiment containing nans in the table, keeping only bacteria with at least minpresence non-nan values
	input:
	expdat : Experiment
	minpresence: int
		minimal number of non-nan samples (keep only if >=)
	output:
	newexp : Experiment
	"""
	params=locals()

	numok=np.sum(np.isfinite(expdat.data),axis=1)
	keep=np.where(numok>=minpresence)[0]
	print(len(keep))
	newexp=hs.reorderbacteria(expdat,keep)
	newexp.filters.append('filternans keep only with >=%d non nan reads' % minpresence)
	hs.addcommand(newexp,"filternans",params=params,replaceparams={'expdat':expdat})
	return newexp
