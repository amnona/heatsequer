#!/usr/bin/env python


"""
heatsequer analysis metrics
"""

# amnonscript

__version__ = "0.9"

import heatsequer as hs

import numpy as np
from collections import defaultdict
import matplotlib as mpl
mpl.use('Qt4Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import *

def calcdist(sample1,sample2,distmetric='bc'):
	"""
	calculate the distance between 2 samples using a given distance metric
	input:
	sample1,sample2 : numpy 1d arrays of floats
		the column arrays for the 2 samples
	distmetric - string
		the distance metric to use:
			'bc' - bray curtis
			'bj' - binary jaccard
			'fbj' - filtered binary jaccard (using a low read threshold)
	output:
	dist : float
		the distance between the 2 samples
	"""
	if distmetric=='bc':
		dist=float(np.sum(np.abs(sample1-sample2)))/(np.sum(sample1)+np.sum(sample2))
	elif distmetric=='bj':
		dist=1-float(np.sum((sample1>0)*(sample2>0)))/np.sum((sample1+sample2)>0)
	elif distmetric=='fbj':
		thresh=1.9
		dist=1-float(np.sum((sample1>thresh)*(sample2>0)+(sample2>thresh)*(sample1>0)))/np.sum((sample1>thresh)+(sample2>thresh))
	else:
		hs.Debug(10,'Distance meteric %s not supported' % distmetric)
		raise
	return dist


def calcdistmat(expdat,distmetric='bc'):
	"""
	calculate the distance matrix between all samples of the experiment
	input:
	expdat : Experiment
	distmetric : string
		the name of the distance metric (see calcdist for options)
	output:
	dist : numpy 2d array
		the pairwise distance matrix
	dsamp : dict
		sample name (key) and position in distance matrix (val)
	"""

	msize=len(expdat.samples)
	dist=np.zeros([msize,msize])
	dsamp={}
	for cs1 in range(msize):
		for cs2 in range(msize):
			dist[cs1,cs2]=calcdist(expdat.data[:,cs1],expdat.data[:,cs2],distmetric)
	for idx,csamp in enumerate(expdat.samples):
		dsamp[csamp]=idx
	return dist,dsamp


def loaddistmat(expdat,dmfilename):
	"""
	load a distance matrix (from qiime) for analysis
	input:
	expdat : Experiment
	dmfilename : string
		name of the qiime distance matrix file

	output:
	distmat : numpy 2d array
		the distance matrix
	dsamp : dict
		the mapping to position in the experiment for each distmat entry
	"""

	fl=open(dmfilename,'rU')
	# get the column ids
	head=fl.readline().strip('\n')
	ids=head.split('\t')
	ids=ids[1:]
	dist=np.array([])
	snames={}
	for idx,cline in enumerate(fl):
		cline=cline.strip("\n")
		vals=cline.split('\t')
		dist=np.vstack((dist,hs.tofloat(vals[1:]))) if dist.size else np.array(hs.tofloat(vals[1:]))
		if not vals[0]==ids[idx]:
			hs.Debug(9,"strange! line %d row head %s but col head %s" % (idx,vals[0],ids[idx]))
		snames[vals[0]]=idx
	fl.close()
	expkeep=[]
	distorder=[]
	dsamp={}
	for idx,csamp in enumerate(expdat.samples):
		if csamp in snames:
			distorder.append(snames[csamp])
			expkeep.append(idx)
			dsamp[csamp]=snames[csamp]
	hs.Debug(6,"%d samples in dist mat, %d samples in experiment" % (len(ids),len(expdat.samples)))
	hs.Debug(6,"%d samples to keep from dist mat, %d samples to keep from experiment" % (len(distorder),len(expkeep)))
	return dist,dsamp


def getgroupdist(expdat,field,distmat,dsamp,plotit=True,plottype='heatmap',uvals=False):
	"""
	calculate the distance matrix based on groups of samples according to field
	using a distance matrix and mapping
	input:
	expdat : Experiment
	field : string
		name of the field to group by
	distmat : numpy 2d arrau
		the distance matrix (from calcdistmat or loaddistmat)
	dsamp : dict
		the mapping of each sample id to the distance matrix position (from calcdistmat or loaddistmat)
	plotit : bool
		True to plot heatmap, False to no plot
	plottype: string
		'heatmap' - to plot heatmap of pairwise values
		'hist' - to plot histogram of pairwise values
	uvals : string
		empty to plot all values, or a list of values to plot only them (in field)
	output:
	gdist : numpy 2d array
		the group distance matrix
	uvals : list
		group names in the matrix (ordered)
	"""

	vals=hs.getfieldvals(expdat,field)
	if not uvals:
		uvals=list(set(vals))
	gdist=np.empty([len(uvals),len(uvals)])
	gdist.fill(np.NaN)
	gmap=defaultdict(list)
	distdict={}
	for idx,cval in enumerate(vals):
		gmap[cval].append(idx)
	for idx1,cg1 in enumerate(uvals):
		pos1=gmap[cg1]
		for idx2,cg2 in enumerate(uvals):
			pos2=gmap[cg2]
			adist=[]
			for p1 in pos1:
				if expdat.samples[p1] not in dsamp:
					continue
				for p2 in pos2:
					if expdat.samples[p2] not in dsamp:
						continue
					if p1==p2:
						continue
					adist.append(distmat[dsamp[expdat.samples[p1]],dsamp[expdat.samples[p2]]])
			distdict[(cg1,cg2)]=adist
			gdist[idx1,idx2]=np.mean(adist)
	if plotit:
		figure()
		if plottype=='heatmap':
			iax=imshow(gdist,interpolation='nearest',aspect='auto',vmin=0,vmax=1)
			ax=iax.get_axes()
			ax.set_xticks(range(len(uvals)))
			ax.set_xticklabels(uvals,rotation=90)
			ax.set_yticks(range(len(uvals)))
			ax.set_yticklabels(uvals)
			title(expdat.studyname+' '+field)
		elif plottype=='hist':
			pl=[]
			pairs=[]
			names=[]
			for k,v in distdict.items():
				ks=set(k)
				if ks in pairs:
					continue
				pl.append(v)
				pairs.append(ks)
				names.append(k)
			hist(pl,alpha=0.5,normed=True,bins=50,range=[0,1])
			legend(names)
	return gdist,uvals

