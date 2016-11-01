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
#mpl.use('Qt4Agg')
import matplotlib.pyplot as plt


def calcdist(sample1,sample2,distmetric='bc',thresh=2):
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
			'logratio' - sum abs log ratio between bacteria with a lower threshold
	thresh : float (optional)
		the cutoff to use for low reads for fbj and logratio metrics (default=2)
	output:
	dist : float
		the distance between the 2 samples
	"""
	if distmetric=='bc':
		dist=float(np.sum(np.abs(sample1-sample2)))/(np.sum(sample1)+np.sum(sample2))
	elif distmetric=='bj':
		dist=1-float(np.sum((sample1>0)*(sample2>0)))/np.sum((sample1+sample2)>0)
	elif distmetric=='fbj':
		dist=1-float(np.sum((sample1>thresh)*(sample2>0)+(sample2>thresh)*(sample1>0)))/np.sum((sample1>thresh)+(sample2>thresh))
	elif distmetric=='logratio':
		sample1=np.copy(sample1)
		sample2=np.copy(sample2)
		sample1[sample1<thresh]=thresh
		sample2[sample2<thresh]=thresh
		sample1=np.log2(sample1)
		sample2=np.log2(sample2)
		dist=float(np.sum(np.abs(sample1-sample2)))
	else:
		hs.Debug(10,'Distance meteric %s not supported' % distmetric)
		raise
	return dist


def calcdistmat(expdat,distmetric='bc',thresh=2):
	"""
	calculate the distance matrix between all samples of the experiment
	input:
	expdat : Experiment
	distmetric : string
		the name of the distance metric (see calcdist for options)
	thresh : float (optional)
		the cutoff to use for low reads for fbj and logratio metrics (default=2)

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
			dist[cs1,cs2]=calcdist(expdat.data[:,cs1],expdat.data[:,cs2],distmetric,thresh=thresh)
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


def getgroupdist(expdat,field,distmat,dsamp,plotit=True,plottype='heatmap',uvals=False,showlegend=True,numbins=50):
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
	showlegend:
		True (default) to show legend for hist, False to hide
	numbins : int
		number of bins for the histogram

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
					# print('missing %s' % expdat.samples[p1])
					continue
				for p2 in pos2:
					if expdat.samples[p2] not in dsamp:
						# print('missing2 %s' % expdat.samples[p1])
						continue
					if p1==p2:
						continue
					cdist=distmat[dsamp[expdat.samples[p1]],dsamp[expdat.samples[p2]]]
					adist.append(cdist)
			distdict[(cg1,cg2)]=adist
			gdist[idx1,idx2]=np.nanmean(adist)
	if plotit:
		plt.figure()
		if plottype=='heatmap':
			plotdistheatmap(gdist,uvals)
			plt.title(expdat.studyname+' '+field)
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
			plt.hist(pl,alpha=0.5,normed=True,bins=numbins,range=[0,1])
			if showlegend:
				plt.legend(names)
	return gdist,uvals


def plotdistheatmap(gdist,uvals,neworder=False,vmin=0,vmax=1):
	"""
	plot a distance heat map and add axis labels
	input:
	gdist : numpy array of float
		the distance matrix (from getgroupdist)
	uvals : list of strings
		the names of the categories (from getgroupdist)
	neworder : list of integers of False
		if not False, the order by which to sort the matrix and labels prior to plotting
	vmin,vmax : int (optional)
		the range for the heatmap
	"""
	if neworder:
		gdist=gdist[neworder,:]
		gdist=gdist[:,neworder]
		uvals=hs.reorder(uvals,neworder)
	plt.figure()
	iax=plt.imshow(gdist,interpolation='nearest',aspect='auto',vmin=vmin,vmax=vmax)
	ax=iax.get_axes()
	ax.set_xticks(range(len(uvals)))
	ax.set_xticklabels(uvals,rotation=90)
	ax.set_yticks(range(len(uvals)))
	ax.set_yticklabels(uvals)
	plt.tight_layout()
	plt.draw()
	plt.colorbar()


def plotdistbar(gdist,uvals,crow=0,neworder=False):
	"""
	plot a distance heat map and add axis labels
	input:
	gdist : numpy array of float
		the distance matrix (from getgroupdist)
	uvals : list of strings
		the names of the categories (from getgroupdist)
	neworder : list of integers of False
		if not False, the order by which to sort the matrix and labels prior to plotting
	"""
	if neworder:
		gdist=gdist[neworder,:]
		gdist=gdist[:,neworder]
		uvals=hs.reorder(uvals,neworder)

	plt.figure()
	plt.bar(np.arange(len(uvals))-0.5,gdist[crow,:])
	ax=plt.gca()
	ax.set_xticks(range(len(uvals)))
	ax.set_xticklabels(uvals,rotation=90)
	plt.tight_layout()
	plt.draw()



def getgroupgroupdist(expdat,field,distmat,dsamp,uvals=False,subfield='host_subject_id',vmin=0,vmax=1):
	"""
	calculate the distance matrix based on groups of samples according to field but calculate seperately for each individual and then combine
	using a distance matrix and mapping
	for the Amina skin cosmetics study

	input:
	expdat : Experiment
	field : string
		name of the field to group by
	distmat : numpy 2d arrau
		the distance matrix (from calcdistmat or loaddistmat)
	dsamp : dict
		the mapping of each sample id to the distance matrix position (from calcdistmat or loaddistmat)
	uvals : string
		empty to plot all values, or a list of values to plot only them (in field)
	subfield : str
		name of the subfield so all distances are calculated seperately for each subfield value (i.e. 'host_subject_id')
	"""
	vals=hs.getfieldvals(expdat,field)
	if not uvals:
		uvals=list(set(vals))
	svals=hs.getfieldvals(expdat,subfield,ounique=True)
	omat=np.zeros([len(uvals),len(uvals)])
	numok=0
	for cval in svals:
		newexp=hs.filtersamples(expdat,subfield,cval)
#		dmap,dmapd=hs.loaddistmat(newexp,'amnon/bray_curtis_armpit-diff-log.txt')
		gdist,uvals=hs.getgroupdist(newexp,field,distmat,dsamp,plotit=False,uvals=uvals)
		gdist[np.isnan(gdist)]=0
		# print(cval)
		# print(gdist)
		if np.isnan(np.sum(np.sum(gdist))):
			continue
		omat=omat+gdist
		numok+=1
	omat=omat/numok
	# print('-----')
	# print(omat)
	plotdistheatmap(omat,uvals,vmin=vmin,vmax=vmax)
	return omat
