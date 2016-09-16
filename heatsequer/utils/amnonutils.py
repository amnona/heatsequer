#!/usr/bin/env python

"""
amnonscript
amnonutils.py

heatsequer
various utility functions
"""

import numpy as np
import scipy as sp

import inspect
import os
import logging
import time

__version__ = "0.2"


def start_log(level=logging.DEBUG, filename='log.hs.log'):
	"""start the logger for the run

	Parameters
	----------
	level : int, optional
		logging.DEBUG, logging.INFO etc. for the log level (between 0-50).
	filename : str, optional
		name of the filename to save the log to or None (default) to print to screen
	"""
	logging.basicConfig(filename=filename, level=level,format='%(asctime)s:%(message)s')
	logger = logging.getLogger(__name__)
	logger.info('*************************')
	logger.info('logging started')


def Debug(dlevel,*args):
	if dlevel>=DebugLevel:
		logger = logging.getLogger(__name__)
		logger.debug(args)
		print (args)


def SetDebugLevel(dlevel):
	"""
	set the debug level for output
	0 - all info (debug)
	5 - warnings / info
	9 - critical

	input:
	dlevel : int (0-10)
		the minimum message level to show
	"""
	global DebugLevel

	DebugLevel = dlevel


def reverse(seq):
	oseq=''
	for a in seq:
		oseq=a+oseq
	return oseq


def complement(seq):
	seq=seq.upper()
	oseq=''
	for a in seq:
		if a=='A':
			oseq+='T'
		elif a=='C':
			oseq+='G'
		elif a=='G':
			oseq+='C'
		elif a=='T':
			oseq+='A'
		else:
			oseq+='N'
	return oseq


def revcomp(seq):
	return reverse(complement(seq))


def iterfastaseqs(filename):
	"""
	iterate a fasta file and return header,sequence
	input:
	filename - the fasta file name

	output:
	seq - the sequence
	header - the header
	"""

	fl=open(filename,"rU")
	cseq=''
	chead=''
	for cline in fl:
		if cline[0]=='>':
			if chead:
				yield(cseq,chead)
			cseq=''
			chead=cline[1:].rstrip()
		else:
			cseq+=cline.strip()
	if cseq:
		yield(cseq,chead)
	fl.close()



def readfastaseqs(filename):
	"""
	read a fasta file and return a list of sequences
	input:
	filename - the fasta file name

	output:
	seqs - a list of sequences
	headers - a list of the headers
	"""
	fl=open(filename,"rU")
	cseq=''
	seqs=[]
	headers=[]
	for cline in fl:
		if cline[0]=='>':
			headers.append(cline[1:].rstrip())
			if cseq:
				seqs.append(cseq)
				cseq=''
		else:
			cseq+=cline.strip()
	if cseq:
		seqs.append(cseq)
	return seqs,headers


def isort(clist,reverse=False):
	"""
	matlab style sort
	returns both sorted list and the indices of the sort
	input:
	clist: a list to sort
	reverse - true to reverse the sort direction
	output:
	(svals,sidx)
	svals - the sorted values
	sidx - the sorted indices
	"""
	res=sorted(enumerate(clist), key=lambda x:x[1],reverse=reverse)
	svals=[i[1] for i in res]
	sidx=[i[0] for i in res]

	return svals,sidx


def tofloat(clist):
	"""
	convert a list of strings to a list of floats
	input:
	clist - list of strings
	output:
	res - list of floats
	"""
	res=[]
	for s in clist:
		try:
			cval=float(s)
		except:
			cval=0
		if np.isnan(cval):
			cval=0
		res.append(cval)
	return res


def reorder(clist,idx):
	""""
	reorder a list according to idx
	"""
	return [clist[i] for i in idx]


def delete(clist,idx):
	"""
	delete elements from list
	"""
	for i in sorted(idx, reverse=True):
		del clist[i]
	return clist


def clipstrings(clist,maxlen,reverse=False):
	"""
	clip all strings in a list to maxlen
	input:
	clist - list of strings
	maxlen - maximal length for each string
	reverse - if true - clip from end (otherwise from beginning)
	"""
	retlist=[]
	for cstr in clist:
		clen=min(maxlen,len(cstr))
		if reverse:
			retlist.append(cstr[-clen:])
		else:
			retlist.append(cstr[0:clen])
	return retlist


def mlhash(cstr,emod=0):
	"""
	do a hash function on the string cstr
	based on the matlab hash function string2hash
	input:
	cstr - the string to hash
	emod - if 0, don't do modulu, otherwise do modulo
	"""
	chash = 5381
	pnum=pow(2,32)-1
	for cc in cstr:
		chash=np.mod(chash*33+ord(cc),pnum)
	if emod>0:
		chash=np.mod(chash,emod)
	return(chash)


def nicenum(num):
	"""
	get a nice string representation of the numnber
	(turn to K/M if big, m/u if small, trim numbers after decimal point)
	input:
	num - the number
	output:
	numstr - the nice string of the number
	"""

	if num==0:
		numstr="0"
	elif abs(num)>1000000:
		numstr="%.1fM" % (float(num)/1000000)
	elif abs(num)>1000:
		numstr="%.1fK" % (float(num)/1000)
	elif abs(num)<0.000001:
		numstr="%.1fu" % (num*1000000)
	elif abs(num)<0.001:
		numstr="%.1fm" % (num*1000)
	else:
		numstr=int(num)
	return numstr



def SeqToArray(seq):
	""" convert a string sequence to a numpy array"""
	seqa=np.zeros(len(seq),dtype=np.int8)
	for ind,base in enumerate(seq):
		if base=='A':
			seqa[ind]=0
		elif base=='a':
			seqa[ind]=0
		elif base=='C':
			seqa[ind]=1
		elif base=='c':
			seqa[ind]=1
		elif base=='G':
			seqa[ind]=2
		elif base=='g':
			seqa[ind]=2
		elif base=='T':
			seqa[ind]=3
		elif base=='t':
			seqa[ind]=3
		elif base=='-':
			seqa[ind]=4
		else:
			seqa[ind]=5
	return(seqa)


def ArrayToSeq(seqa):
	""" convert a numpy array to sequence (upper case)"""
	seq=''
	for cnuc in seqa:
		if cnuc==0:
			seq+='A'
		elif cnuc==1:
			seq+='C'
		elif cnuc==2:
			seq+='G'
		elif cnuc==3:
			seq+='T'
		else:
			seq+='N'
	return(seq)



def fdr2(p):
	"""Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
	p = np.asfarray(p)
	by_descend = p.argsort()[::-1]
	by_orig = by_descend.argsort()
	steps = float(len(p)) / np.arange(len(p), 0, -1)
	q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
	return q[by_orig]



def fdr(pvalues, correction_type = "Benjamini-Hochberg"):
	"""
	consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1])
	"""

	pvalues = np.array(pvalues)
	n = float(pvalues.shape[0])
	new_pvalues = np.empty(n)
	if correction_type == "Bonferroni":
		new_pvalues = n * pvalues
	elif correction_type == "Bonferroni-Holm":
		values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
		values.sort()
		for rank, vals in enumerate(values):
			pvalue, i = vals
			new_pvalues[i] = (n-rank) * pvalue
	elif correction_type == "Benjamini-Hochberg":
		values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
		values.sort()
		values.reverse()
		new_values = []
		for i, vals in enumerate(values):
			rank = n - i
			pvalue, index = vals
			new_values.append((n/rank) * pvalue)
		for i in range(0, int(n)-1):
			if new_values[i] < new_values[i+1]:
				new_values[i+1] = new_values[i]
		for i, vals in enumerate(values):
			pvalue, index = vals
			new_pvalues[index] = new_values[i]
	return new_pvalues


def common_start(sa,sb):
	"""
	returns the longest common substring from the beginning of sa and sb
	from http://stackoverflow.com/questions/18715688/find-common-substring-between-two-strings
	"""

	def _iter():
		for a, b in zip(sa, sb):
			if a == b:
				yield a
			else:
				return
	return ''.join(_iter())


DebugLevel=5


def listdel(dat,todel):
	"""
	delete elements with indices from list todel in the list dat
	input:
	dat - the list to remove elements from
	todel - indices of the items to remove

	output:
	dat - the new deleted list
	"""

	for cind in sorted(todel, reverse=True):
		del dat[cind]
	return dat



def listtodict(dat):
	"""
	convert a list into a dict with keys as elements, values the position in the list
	input:
	dat - the list

	output:
	thedict
	"""

	thedict={}
	for idx,cdat in enumerate(dat):
		if cdat in thedict:
			thedict[cdat].append(idx)
		else:
			thedict[cdat]=[idx]
	return thedict


def savelisttofile(dat,filename,delimiter='\t'):
	"""
	save a list to a (tab delimited) file
	inputL
	dat - the list to save
	filename - the filename to save to
	delimiter - the delimiter to use
	"""

	with open(filename,'w') as fl:
		fl.write(delimiter.join(dat))


def dictupper(dat):
	"""
	turn dict keys to upper case
	input:
	dat - a dict with string keys
	output:
	newdat - a dict with the upper case keys
	"""

	newdat = {k.upper(): v for k,v in dat.iteritems()}
	return newdat


def listupper(dat):
	"""
	turn a list of strings to upper case
	input:
	dat : list of strings
	output:
	newdat : list of strings
		- in uppercase
	"""

	newdat = [cstr.upper() for cstr in dat]
	return newdat


def get_current_data_path(fn, subfolder='data'):
	"""Return path to filename ``fn`` in the data folder.
	During testing it is often necessary to load data files. This
	function returns the full path to files in the ``data`` subfolder
	by default.
	Parameters
	----------
	fn : str
		File name.
	subfolder : str, defaults to ``data``
		Name of the subfolder that contains the data.
	Returns
	-------
	str
		Inferred absolute path to the test data for the module where
		``get_data_path(fn)`` is called.
	Notes
	-----
	The requested path may not point to an existing file, as its
	existence is not checked.

	Taken from scikit-bio (Thanks!)
	"""
	# getouterframes returns a list of tuples: the second tuple
	# contains info about the caller, and the second element is its
	# filename
	callers_filename = inspect.getouterframes(inspect.currentframe())[1][1]
	path = os.path.dirname(os.path.abspath(callers_filename))
	data_path = os.path.join(path, subfolder, fn)
	return data_path


def findn(text,substr,num):
	"""
	find num-th occurance of substr in text
	input:
	text : string
		the string to search in
	substr : string
		the substring to search for in text
	num : int
		the occurance number (1 is the first occurance, etc)

	output:
	index : int
		the position of the start of the num-th substring in text, or -1 if not present
	"""
	index=0
	while index < len(text):
		index = text.find(substr, index)
		if index == -1:
			break
		num-=1
		if num==0:
			return index
		index+=len(substr)
	return -1



def getnicetax(name,separator=';'):
	"""
	get the last non empty string (separated by separator)
	used to get a nice taxonomy name from a taxonomy string

	input:
	name : str
		the taxonomy string
	separator: str
		the separator between taxonomic levels (i.e. ';')

	output:
	nicename : str
		only the last non empty part of name
	"""
	nicename='unknown'
	s=name.split(separator)
	for cstr in s:
		if len(cstr)>0:
			nicename=cstr
	return nicename


def sum(data,axis=None):
	"""
	Get the sum of the numpy array/scipy sparse matrix

	input:
	data : numpy array or scipy matrix
		the data to sum
	axis : int or None
		equivalent to axis parameter on numpy sum

	output:
	csum : a 1d numpy array
		the sum
	"""

	# if sparse matrix, the output is a matrix
	# so we need to convert to 1d array
	if sp.sparse.isspmatrix(data):
		csum=data.sum(axis=axis).A.flatten()
	else:
		csum=np.sum(data,axis=axis)
	return csum


def mean(data,axis=None):
	"""
	Get the mean of the numpy array/scipy sparse matrix

	input:
	data : numpy array or scipy matrix
		the data to sum
	axis : int or None
		equivalent to axis parameter on numpy sum

	output:
	cmean : a 1d numpy array
		the mean
	"""

	# if sparse matrix, the output is a matrix
	# so we need to convert to 1d array
	if sp.sparse.isspmatrix(data):
		cmean=data.mean(axis=axis).A.flatten()
	else:
		cmean=np.mean(data,axis=axis)
	return cmean


def median(data,axis=None):
	"""
	Get the mean of the numpy array/scipy sparse matrix

	input:
	data : numpy array or scipy matrix
		the data to sum
	axis : int or None
		equivalent to axis parameter on numpy sum

	output:
	cmedian : a 1d numpy array
		the mean
	"""

	# if sparse matrix, the output is a matrix
	# so we need to convert to 1d array
	if sp.sparse.isspmatrix(data):
		Debug(9,'Median not supported for sparse matrix. come back one year')
		cdat=data.todense().A[0]
		cmedian=np.median(cdat,axis=axis)
	else:
		cmedian=np.median(data,axis=axis)
	return cmedian


def divvec(data,vec):
	"""
	divide the data matrix data by the vector vec
	works for sparse and non sparse data types

	input:
	data : sparse matrix or numpy array
		the data to divide each colum by the vector element
	vec : numpy array
		the vector to devide by

	output:
	data : same as data
		each element in each column in data divided by the corresponding vec element
	"""
	if sp.sparse.isspmatrix(data):
		# for cpos,cval in data.items():
		# 	data[cpos]=cval/vec[cpos[1]]
		numcols=data.shape[1]
		b=sp.sparse.lil_matrix( (numcols,numcols) )
		for idx in range(numcols):
			b[idx,idx]=1.0/vec[idx]
		b=b.tocsr()
		data=data*b
	else:
		data=data/vec

	return data


def log2(data,minthresh=2):
	"""
	calculate the log2 of the data in dense or sparse format. numbers below minthresh are rounded to minthresh

	input:
	data : numpy array of sparse matrix
		the 2d array to calculate the log2
	minthresh: float
		the minimal number of reads
		(data numbers below minthresh are rounded to minthresh)

	output:
	ldat : numpy array or sparse matrix
		log2 of the data matrix after thresholding. output format is similar to input
	"""
	if sp.sparse.isspmatrix(data):
		pass
	else:
		data[data<minthresh]=minthresh
		ldat=np.log2(data)
	return ldat
