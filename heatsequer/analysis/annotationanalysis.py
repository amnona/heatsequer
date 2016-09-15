#!/usr/bin/env python


"""
heatsequer manual database analysis module
"""

# amnonscript

__version__ = "0.9"

import heatsequer as hs

import collections
import numpy as np
from scipy import stats
from scipy import spatial
import matplotlib as mpl
#mpl.use('Qt4Agg')
import matplotlib.pyplot as plt
#from matplotlib.pyplot import *
import sklearn.metrics
import sklearn.cross_validation
from sklearn.cluster import AffinityPropagation
from sklearn.preprocessing import scale
from sklearn import metrics
import copy
from pdb import set_trace as XXX


def getexpannotations(expdat,usecooldb=True,usesupercooldb=True):
	"""
	init the curations field of the experiment by getting curations for all bacteria in the experiment
	NOTE: inplace

	input:
	expdat : Experiment
	usecooldb : bool
		True (default) to get cooldb annotations, False to skip
	usesupercooldb : bool
		True (default) to get supercooldb annotations, False to skip

	output:
	expdat : Experiment
		with the seqannotations field initialized to dict {sampled: [annotation1, annotation2, ...]}
		also initialize the annotationseqs dict to {annotation: [seq1,seq2,...]}
	"""
	hs.Debug(2,'Adding annotations')
	seqannotations={}
	annotationseqs=collections.defaultdict(list)
	totcount=0
	for cseq in expdat.seqs:
		seqannotations[cseq]=[]
		info=[]
		if usecooldb:
			info.extend(hs.cooldb.getseqinfo(hs.cdb,cseq))
		if usesupercooldb:
			info.extend(hs.supercooldb.getcurationstrings(hs.scdb,cseq))
		if len(info)>0:
			seqannotations[cseq].extend(info)
			for cinfo in info:
				annotationseqs[cinfo].append(cseq)
				totcount+=1
	hs.Debug(5,'found %d total annotations, %d unique' % (totcount,len(annotationseqs)))
	expdat.seqannotations=seqannotations
	expdat.annotationseqs=annotationseqs
	return(expdat)


def annotationenrichment(expdat,seqs,compareseqs=None,fdrval=0.1):
	"""
	get a list of annotations enriched in seqs compared to random draw from expdat

	input:
	expdat : Experiment
	seqs : list of sequences ('ACGT')
		the sequences in which to test the enrichment
	compareseqs : list of sequences ('ACGT')
		the sequences to compare to
		None (default) to compare to all the experiment

	output:
	newplist - a sorted list of dict for annotaions which are below fdr ('description','pval','observed','expected')
	"""
	# if annotations not initialized - get them (to save time)
	if expdat.seqannotations is None:
		expdat=hs.getexpannotations(expdat)
	if expdat.annotationseqs is None:
		expdat=hs.getexpannotations(expdat)

	# count the number of annotations for each term in the group
	# into a dict {term:total number observed in group}
	groupannotationcount=collections.defaultdict(int)
	totgroup=0
	for cseq in seqs:
		for cinfo in expdat.seqannotations[cseq]:
			groupannotationcount[cinfo]+=1
			totgroup+=1

	# count the number of annotations for each term in the comparison group
	# into a dict {term:total number observed in comparison group}
	if compareseqs is None:
		compareseqs=expdat.seqs
	compgroupannotationcount=collections.defaultdict(int)
	totcompgroup=0
	for cseq in compareseqs:
		for cinfo in expdat.seqannotations[cseq]:
			compgroupannotationcount[cinfo]+=1
			totcompgroup+=1

	hs.Debug(6,'%d annotations in group, %d in all' % (totgroup,totcompgroup))
	# calculate the probability per term
	# note: we use a bad calculation (can choose the same term twice for a single bacteria). need to improve (permutations?)
	pvals={}
	allp=[]
	pv=[]
	for cinfo in expdat.annotationseqs.keys():
		pcompgroup=float(compgroupannotationcount[cinfo])/totcompgroup
		pval1=stats.binom.cdf(groupannotationcount[cinfo],totgroup,pcompgroup)
		pval2=stats.binom.cdf(totgroup-groupannotationcount[cinfo],totgroup,1-pcompgroup)
		p=np.min([pval1,pval2])
		# p=pval1
		pvals[cinfo]=p

		allp.append(p)
		cpv={}
		cpv['pval']=p
		cpv['observed']=groupannotationcount[cinfo]
		cpv['expected']=pcompgroup*totgroup
		cpv['description']=cinfo
		pv.append(cpv)

	fval=hs.fdr(allp)
	keep=np.where(np.array(fval)<=fdrval)
	plist=[]
	rat=[]
	for cidx in keep[0]:
		plist.append(pv[cidx])
		rat.append(np.abs(float(pv[cidx]['observed']-pv[cidx]['expected']))/np.mean([pv[cidx]['observed'],pv[cidx]['expected']]))
	si=np.argsort(rat)
	si=si[::-1]
	newplist=[]
	for idx,crat in enumerate(rat):
		newplist.append(plist[si[idx]])
	for cp in newplist:
		if cp['observed']>cp['expected']:
			hs.Debug(6,cp)

	return(newplist)
