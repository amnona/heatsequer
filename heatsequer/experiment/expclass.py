#!/usr/bin/env python


"""
heatsequer experiment class
"""

# amnonscript

__version__ = "0.9"

import heatsequer as hs

import os
import copy
import numpy as np
from pdb import set_trace as XXX
import time


class Experiment:
	'''
	experiment class holds the read data and metadata about the experiment, as well
	as the command history
	'''

	# the static unique experiment id
	experimentid=0

	def __init__(self):
		# the data matrix (non sparse)
		self.data=[]

		# the sample dictionary (double hash - sampleid and then mapping file field)
		self.smap=[]

		# name of all the fields in the mapping data
		self.fields=[]

		# list of sampleids ordered according to the data matrix
		self.samples=[]

		# the sequences in the table
		self.seqs=[]

		# dictionary holding all the sequences and their position (for fast lookup)
		self.seqdict={}

		# taxonomies
		self.tax=[]

		# the hashed ids for the sequences
		self.sids=[]

		# the original name for each otu (from the biom table)
		self.origotunames=[]

		# original table name
		self.tablefilename=''

		# original mapping file
		self.mapfilename=''

		# name of the study (or the table file name without path)
		self.studyname=''

		# number of reads for each sample in the biom table
		self.origreads=[]

		# the history of actions performed
		self.filters=[]

		# and the command list
		self.commands=[]

		# the complete sequence database
		self.seqdb=None

		# the cool sequences (manually curated) database
		self.cdb=None

		# the list of annotations to add to plot (for addplotmetadata)
		self.plotmetadata=[]

		# the tree structure of the sequences (from loadexptree)
		self.tree=False

		# the experiment type ('biom' or 'meta' for metabolite)
		self.datatype=''

		# the unqiue experiment id
		self.uniqueid=0

		# the md5 of the original data and mapping files loaded
		# used for a unique id for the data in the manual curation database
		self.datamd5=''
		self.mapmd5=''

		hs.Debug(0,'New experiment initialized')

	# get a unique identifier and increase by 1
	def getexperimentid(self):
		Experiment.experimentid+=1
		return Experiment.experimentid



def copyexp(expdat):
	"""
	copy an experiment (duplicating the important fields)
	but give it a unique identifier

	Parameters
	----------
	expdat : Experiment
		the experiment to copy
	output:
	newexp : Experiment
		a deep copy of expdat
	"""

	newexp=copy.copy(expdat)
	newexp.data=copy.deepcopy(expdat.data)
	newexp.smap=copy.deepcopy(expdat.smap)
	newexp.fields=copy.deepcopy(expdat.fields)
	newexp.samples=copy.deepcopy(expdat.samples)
	newexp.seqs=copy.deepcopy(expdat.seqs)
	newexp.seqdict=copy.deepcopy(expdat.seqdict)
	newexp.tax=copy.deepcopy(expdat.tax)
	newexp.sids=copy.deepcopy(expdat.sids)
	newexp.origotunames=copy.deepcopy(expdat.origotunames)
	newexp.tablefilename=copy.deepcopy(expdat.tablefilename)
	newexp.mapfilename=copy.deepcopy(expdat.mapfilename)
	newexp.studyname=copy.deepcopy(expdat.studyname)
	newexp.origreads=copy.deepcopy(expdat.origreads)
	newexp.filters=copy.deepcopy(expdat.filters)
	newexp.commands=copy.deepcopy(expdat.commands)
	newexp.plotmetadata=copy.deepcopy(expdat.plotmetadata)
#	nmewexp.tree=copy.deepcopy(expdat.tree)
	newexp.datatype=copy.deepcopy(expdat.datatype)

	# get a unique identifier for this experiment
	newexp.uniqueid=newexp.getexperimentid()
	return newexp


def hashseq(seq):
	'''
	calculate the hash value for a given sequence (for an almost unique sequence identifier)
	used for the sid field in experiment
	input:
	seq : str
		the sequence to hash

	output:
	hval : int
		the hash value for the sequence
	'''
	hval=hs.mlhash(seq, emod=10000000)
	return hval



def addcommand(expdat,command,params={},replaceparams={}):
	'''
	append a command string to the experiment command list from the command and the unique experiment id
	"expXXX=command" where XXX is the unique experiment id
	if params is supplied, use them as the function parameters, otherwise just use command

	input:
	expdat : experiment
		the experiment for which to prepare the command
	command : str
		the command
	params : dict
		if empty, just append the command
		if dict, append command+"("+params+")"
	replaceparams : dict
		a dict of parameters whos values need to be replaced by an experimentid.
		key is parameter, value is experiment, from where the experimentid will be taken
	'''
	newcommand='exp%d=hs.%s' % (expdat.uniqueid,command)
	if len(params)>0:
#		if replaceparams:
#			for rk,rv in replaceparams.items():
#				if rk not in params:
#					hs.Debug(9,'replacement parameter %s not in params' % rk)
#				params[rk]='exp%d' % rv.uniqueid
		newcommand+='('
		for k,v in params.items():
			if k in replaceparams:
				v='exp%d' % v.uniqueid
			else:
				v=repr(v)
			newcommand+='%s=%s,' % (k,str(v))
		newcommand=newcommand[:-1]+')'
	expdat.commands.append(newcommand)


def reordersamples(exp,newpos,inplace=False):
	"""
	reorder the samples of the experiment
	input:
	exp - the experiment
	newpos - array - the new positions (can skip positions to delete them)
	output:
	newexp - the new experiment
	"""

	if inplace:
		newexp=exp
	else:
		newexp=copyexp(exp)
#		newexp=copy.deepcopy(exp)
	newexp.data=newexp.data[:,newpos]
	newexp.samples=hs.reorder(newexp.samples,newpos)
	newexp.origreads=hs.reorder(newexp.origreads,newpos)
	return newexp

def reorderbacteria(exp,order,inplace=False):
	"""
	reorder the bacteria in an experiment (can delete if bacteria not in new order)
	input:
	exp - the experiment
	order - the new order
	output:
	newexp
	"""
	if inplace:
		newexp=exp
	else:
		newexp=copyexp(exp)
#		newexp=copy.deepcopy(exp)
	newexp.data=newexp.data[order,:]
	newexp.seqs=hs.reorder(newexp.seqs,order)
	newexp.seqdict={}
	for idx,cseq in enumerate(newexp.seqs):
		newexp.seqdict[cseq]=idx
	newexp.tax=hs.reorder(newexp.tax,order)
	newexp.sids=hs.reorder(newexp.sids,order)
	return newexp


def getfieldvals(expdat,field,ounique=False):
	"""
	get a list of the field values in all samples
	input:
	expdat : Experiment
	field : string
		name of the field to get the values from
	ounique : bool
		True to get unique values, False to get all
	"""
	vals=[]
	for cid in expdat.samples:
		vals.append(expdat.smap[cid][field])
	if ounique:
		vals=list(set(vals))
	return vals


def joinfields(expdat,field1,field2,newfield):
	"""
	join 2 fields to create a new field for each sample
	input:
	expdat : Experiment
	field1,field2 : string
		name of the 2 fields to join
	newfield : string
		name of new field to add
	"""
	params=locals()

	for csamp in expdat.samples:
		expdat.smap[csamp][newfield]=expdat.smap[csamp][field1]+';'+expdat.smap[csamp][field2]
	expdat.fields.append(newfield)
	expdat.filters.append("join fields %s, %s to new field %s" % (field1,field2,newfield))
	hs.addcommand(expdat,"joinfields",params=params,replaceparams={'expdat':expdat})
	return expdat


def joinexperiments(exp1,exp2,missingval='NA',origfieldname='origexp'):
	"""
	join 2 experiments into a new experiment. adding a new field origfieldname
	input:
	exp1,exp2 - the experiments to join
	missingval - string to put when field not in mapping file of one of the experiments
	origfieldname - name of the new field to add which contains the original experiment name
	"""
	params=locals()

	# join the sequences of both experiments
	# ASSUMING SAME SEQ LENGTH!!!!
	allseqs=list(set(exp1.seqs) | set(exp2.seqs))
	alldict={}
	alltax=[]
	allids=[]
	for idx,cseq in enumerate(allseqs):
		alldict[cseq]=idx
	# make the new joined data for each experiment
	dat1=np.zeros((len(allseqs),np.size(exp1.data,1)))
	for idx,cseq in enumerate(allseqs):
		if cseq in exp1.seqdict:
			dat1[idx,:]=exp1.data[exp1.seqdict[cseq],:]
			alltax.append(exp1.tax[exp1.seqdict[cseq]])
			allids.append(exp1.sids[exp1.seqdict[cseq]])
		else:
			alltax.append(exp2.tax[exp2.seqdict[cseq]])
			allids.append(exp2.sids[exp2.seqdict[cseq]])

	dat2=np.zeros((len(allseqs),np.size(exp2.data,1)))
	for idx,cseq in enumerate(allseqs):
		if cseq in exp2.seqdict:
			dat2[idx,:]=exp2.data[exp2.seqdict[cseq],:]

	newexp=hs.copyexp(exp1)
	# concatenate the reads
	newexp.data=np.concatenate((dat1,dat2), axis=1)
	newexp.seqdict=alldict
	newexp.seqs=allseqs
	newexp.tax=alltax
	newexp.sids=allids
	newexp.sids=newexp.seqs
	newexp.samples = list(exp1.samples) + list(exp2.samples)
	newexp.origreads=exp1.origreads+exp2.origreads
	newexp.fields=list(set(exp1.fields+exp2.fields))

	for cfield in newexp.fields:
		if cfield in exp1.fields:
			continue
		for csamp in exp1.samples:
			newexp.smap[csamp][cfield]=missingval

	for csamp in exp2.samples:
		newexp.smap[csamp]={}
		for cfield in newexp.fields:
			if cfield in exp2.fields:
				newexp.smap[csamp][cfield]=exp2.smap[csamp][cfield]
			else:
				newexp.smap[csamp][cfield]=missingval

	for csamp in exp1.samples:
		newexp.smap[csamp][origfieldname]=exp1.studyname
	for csamp in exp2.samples:
		newexp.smap[csamp][origfieldname]=exp2.studyname
	newexp.fields.append(origfieldname)

	newexp.filters.append('joined with %s' % exp2.studyname)
	hs.addcommand(newexp,"joinexperiments",params=params,replaceparams={'exp1':exp1,'exp2':exp2})
	return newexp



def clipseqs(expdat,startpos,addseq='TAC'):
	"""
	clip the first nucleotides in all sequences in experiment
	to fix offset in sequencing
	input:
	expdat
	startpos - the position to start from (0 indexed) or negative to add nucleotides
	addseq - the sequence to add (just a guess) if startpos is negative
	output:
	newexp - new experiment with all sequences clipped and joined identical sequences
	"""
	params=locals()

	newexp=copy.deepcopy(expdat)
	newseqs=[]
	newdict={}
	keeppos=[]
	for idx,cseq in enumerate(newexp.seqs):
		if startpos>=0:
			cseq=cseq[startpos:]
		else:
			cseq=addseq[:abs(startpos)]+cseq
			cseq=cseq[:len(expdat.seqs[0])]
		if cseq in newdict:
			newexp.data[newdict[cseq],:] += newexp.data[idx,:]
		else:
			newdict[cseq]=idx
			newseqs.append(cseq)
			keeppos.append(idx)
	newexp=reorderbacteria(newexp,keeppos)
	newexp.seqs=newseqs
	newexp.seqdict=newdict
	hs.addcommand(newexp,"clipseqs",params=params,replaceparams={'expdat':expdat})
	newexp.filters.append("trim %d nucleotides" % startpos)
	return newexp


def findsamples(expdat,field,value,exclude=False):
	"""
	return the positions of samples in expdat matching value in field
	similar to filtersamples but returns a list of indices (for the data matrix)
	input:
	expdat
	field - name of the field to test
	value - the value to look for (or a list of values)
	exclude - True to get positions without that value, False to get positions of the value

	output:
	pos - a list of positions matching the field/val (for use as indices in expdat.data)
	"""

	pos=[]
	if not isinstance(value,list):
		value=[value]

	for cidx,csamp in enumerate(expdat.samples):
		if expdat.smap[csamp][field] in value:
			if not exclude:
				pos.append(cidx)
		else:
			if exclude:
				pos.append(cidx)
	return pos



def insertbacteria(expdat,freqs=[],seq="unknown",tax="unknown",logit=True):
	"""
	insert a new bacteria to an experiment

	input:
	expdat
	freqs - the frequency of the bacteria in all samles of expdat or False to add zeros
	seq - the sequence of the new bacteria
	tax - taxonomy of the new bacteria
	logit - True to add command log/filter, False to not add (if called from other function)

	output:
	pos - position of the new bacteria
	"""
	params=locals()

	if len(freqs)==0:
		freqs=np.zeros([1,len(expdat.seqs)])

	expdat.data=np.vstack((expdat.data,freqs))
	expdat.tax.append(tax)

	if seq in expdat.seqdict:
		hs.Debug(6,'Sequence already in experiment',seq)
	# get a unique sequence
		cid=0
		while seq+str(cid) in expdat.seqdict:
			cid+=1
#		expdat.seqs.append()
		seq=seq+str(cid)

	expdat.seqs.append(seq)
	expdat.seqdict[seq]=len(expdat.seqs)-1
	expdat.sids.append(seq)
	if logit:
		expdat.filters.append("insert bacteria")
		hs.addcommand(expdat,"insertbacteria",params=params,replaceparams={'expdat':expdat})
	return expdat,len(expdat.seqs)-1


def addsubtrees(expdat,tree,inplace=False):
	"""
	add otus for all subtrees with the frequency being the sum of all bacteria in the subtree
	input:
	expdat - the experiment
	tree - the tree for the experiment
	inplace - if true, replace current experiment

	output:
	newexp - the new experiment with twice-1 number of otus
	"""
	params=locals()
#	if not expdat.tree:
#		hs.Debug(8,"No tree loaded for experiment")
#		return False

	if inplace:
		newexp=expdat
	else:
		newexp=hs.copyexp(expdat)

	subtrees=tree.subsets()
	for csubtree in subtrees:
		newname=""
		newtax=""
		numuse=0
		newfreq=np.zeros([1,len(newexp.samples)])
		for cbact in csubtree:
			if cbact not in newexp.seqdict:
				hs.Debug(4,'sequence not in seqdict',cbact)
				continue
			numuse+=1
			cpos=newexp.seqdict[cbact]
			newfreq+=newexp.data[cpos,:]
			newname+='%d,' % cpos
			if newtax=='':
				newtax=newexp.tax[cpos]
			else:
				newtax=hs.common_start(newtax,newexp.tax[cpos])
		# add only if we have 2 bacteria or more
		if numuse>1:
			if newname not in newexp.seqdict:
				newexp,newpos=insertbacteria(newexp,freqs=newfreq,seq=newname,tax=newtax,logit=False)

	newexp.filters.append("Add subtrees")
	hs.addcommand(newexp,"addsubtrees",params=params,replaceparams={'expdat':expdat})
	return(newexp)



def findseqsinexp(expdat,seqs):
	"""
	find sequences from seqs in expdat sequences and return the indices
	input:
	expdat
	seqs - a list of sequences

	output:
	res - a list of indices where seqs are in expdat sequences
	"""
	res=[]
	for cseq in seqs:
		res.append(expdat.seqdict[cseq])
	return res


def samplemeanpervalue(expdat,field):
	"""
	BETTER TO USE filtersimilarsamples!!!!
	create a new experiment, with 1 sample per value in field, containing the mean of all samples with that value

	input:
	expdat : Experiment
	field : string
		the field to use (i.e. 'ENV_MATTER')

	output:
	newexp : Experiment
		The new experiment with 1 sample per unique value of field
	"""
	params=locals()

	uvals=hs.getfieldvals(expdat,field,ounique=True)
	vals=hs.getfieldvals(expdat,field,ounique=False)

	vdict=hs.listtodict(vals)
	nsamps=[]
	for cval in uvals:
		nsamps.append(vdict[cval][0])
	newexp=hs.reordersamples(expdat,nsamps)
	for idx,cval in enumerate(uvals):
		cdat=expdat.data[:,vdict[cval]]
		mv=np.mean(cdat,axis=1)
		newexp.data[:,idx]=mv
	newexp.filters.append('samplemeanpervalue for field %s' % field)
	hs.addcommand(newexp,"samplemeanpervalue",params=params,replaceparams={'expdat':expdat})
	return(newexp)


def convertdatefield(expdat,field,newfield,timeformat='%m/%d/%y %H:%M'):
	"""
	convert a field containing date/time to a numeric (seocds since epoch) field (create a new field for that)
	input:
	expdat : Experiment
		the experiment to add the field to
	field : string
		name of the field containing the date/time format
	newfield : string
		name of the new field (with seconds since epoch)
	timeformat : string
		format of the date/time field (based on time format)
	output:
	newexp : Experiment
		the experiment with the added time since epoch field
	"""
	params=locals()

	newexp=hs.copyexp(expdat)
	newexp.fields.append(newfield)
	for csamp in newexp.samples:
		newexp.smap[csamp][newfield]=time.mktime(time.strptime(newexp.smap[csamp][field],timeformat))
	newexp.filters.append('add time field %s (based on field %s)' % (newfield,field))
	hs.addcommand(newexp,"convertdatefield",params=params,replaceparams={'expdat':expdat})
	return(newexp)


def fieldtobact(expdat,field,bactname='',meanreads=1000,cutoff=0):
	"""
	convert values in a map file field to a new bacteria (to facilitate numeric analysis)
	input:
	expdat : Experiment
	field : string
		name of the field to convert
	bactname : string
		name of the new bacteria (empty to have similar to field name)
	meanreads : int
		the mean number of reads for the new field bacteria
	cutoff : int
		the minimal value of the field per sample (otherwise replace with meanreads)

	output:
	newexp : Experiment
		with added bacteria with the field vals as reads
	"""
	params=locals()

	if len(bactname)==0:
		bactname=field
	fv=hs.getfieldvals(expdat,field)
	vals=np.array(hs.tofloat(fv))
	okpos=np.where(vals>=cutoff)[0]
	badpos=np.where(vals<cutoff)[0]
	scalefactor=np.mean(vals[okpos])
	vals[okpos]=(vals[okpos]/scalefactor)*meanreads
	vals[badpos]=meanreads
	newexp=hs.copyexp(expdat)
	hs.insertbacteria(newexp,vals,bactname,bactname,logit=False)
	newexp.filters.append('add bacteria from map field %s' % field)
	hs.addcommand(newexp,"fieldtobact",params=params,replaceparams={'expdat':expdat})
	return(newexp)


def get_data_path(fn, subfolder='data'):
	"""
	Return path to filename ``fn`` in the data folder.
	returns the joining of the heatsequerdir variable (set in __init__) and the subfolder and fn
	"""

	return os.path.join(hs.heatsequerdir,subfolder,fn)

