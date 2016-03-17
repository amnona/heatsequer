#!/usr/bin/env python


"""
heatsequer experiment sorting module
"""

# amnonscript

__version__ = "0.9"

import heatsequer as hs

import numpy as np
import matplotlib as mpl
#mpl.use('Qt4Agg')
import matplotlib.pyplot as plt


def findmislabels(expdat,field,distmetric='bc'):
	""""
	find mislabelled samples according to field
	input:
	expdat : Experiment
	field : string
		name of the field to examine (i.e. subjectid)
	distmetric : string
		the distance meteric to use (see calcdist)
	"""

	expdat=hs.sortsamples(expdat,field)
	fvals=hs.getfieldvals(expdat,field)
	ufvals=list(set(fvals))
	onames=[]
	for idx,csamp in enumerate(expdat.samples):
		onames.append(csamp+';'+fvals[idx])
	omat=np.zeros([len(fvals),len(ufvals)])
	for groupidx,groupval in enumerate(ufvals):
		cexp=hs.filtersamples(expdat,field,groupval,exact=True)
		for aidx,aval in enumerate(expdat.samples):
			cdist=[]
			for gidx,gval in enumerate(cexp.samples):
				# don't measure distance to ourselves
				if gval==aval:
					continue
				cdist.append(hs.calcdist(cexp.data[:,gidx],expdat.data[:,aidx],distmetric=distmetric))
			omat[aidx,groupidx]=np.mean(cdist)
	plt.figure()
	iax=plt.imshow(omat,interpolation='nearest',aspect='auto')
	ax=iax.get_axes()
	ax.set_xticks(range(len(ufvals)))
	ax.set_xticklabels(ufvals,rotation=90)
	ax.set_yticks(range(len(onames)))
	ax.set_yticklabels(onames)


def plotnumotushist(expdat,newfig=True,threshold=0.001):
	"""
	normalized histogram of number of OTUs per sample (alpha diversity)
	input:
	expdat : Experiment
	newfig : bool
		True to plot a new figure, False to plot on current figure
	threshold : float
		minimal number of reads for an otu to be present
	"""

	notus=[]
	for idx,csamp in enumerate(expdat.samples):
		cnotus=np.sum(expdat.data[:,idx]>=threshold)
		notus.append(cnotus)
	if newfig:
		plt.figure()
	plt.hist(notus,bins=20,range=[0,200],normed=True)
	plt.title('num OTUs for %s' % expdat.tablefilename)
