#!/usr/bin/env python

"""
amnonscript
amnonutils.py

various utility functions
"""

import sys
import numpy as np

from sys import getsizeof, stderr
from itertools import chain
from collections import deque
try:
	from reprlib import repr
except ImportError:
	pass

__version__ = "0.2"

def Debug(dlevel,*args):
	if dlevel>=DebugLevel:
		print (args)


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
			res.append(float(s))
		except:
			res.append(0)
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
		for i in xrange(0, int(n)-1):
			if new_values[i] < new_values[i+1]:
				new_values[i+1] = new_values[i]
		for i, vals in enumerate(values):
			pvalue, index = vals
			new_pvalues[index] = new_values[i]
	return new_pvalues



DebugLevel=5

