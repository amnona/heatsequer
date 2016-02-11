#!/usr/bin/env python


"""
heatsequer experiment sorting module
"""

# amnonscript

__version__ = "0.9"

import heatsequer as hs

import numpy as np
import copy
from sklearn.preprocessing import scale
from scipy import cluster
from scipy import spatial



def sortbacteria(exp,inplace=False,logit=True):
	"""
	sort bacteria according to taxonomy (alphabetically)

	input:
	exp : experiment
		the experiment to sort
	inplace : bool
		True to sort in place (replace current experiment), False to create a new experiment
	logit : bool
		True to add to command log, False to skip (if called from other heatsequer function)

	output:
		newexp : experiment
			The sorted experiment (by taxonomy name)
	"""
	params=locals()

	tax=exp.tax
	svals,sidx=hs.isort(tax)
	newexp=hs.reorderbacteria(exp,sidx,inplace=inplace)
	if logit:
		newexp.filters.append('sorted bacteria by taxonomy')
		hs.addcommand(newexp,"sortbacteria",params=params,replaceparams={'exp':exp})
	return newexp


def clusterbacteria(exp,minreads=0,uselog=True):
	"""
	cluster bacteria in an experiment according to similar behavior
	input:
	exp : Experiment
	minreads : int
		the minimal number of reads to keep before clustering (to make faster)
	uselog : bool
		True to log transform reads for clustering (before normalizing), false to use full reads

	output:
	newexp : Experiment
		the filtered and clustered experiment
	"""
	params=locals()

	newexp=hs.filterminreads(exp,minreads,logit=False)
	# normalize each row (bacteria) to sum 1
	dat=copy.copy(newexp.data)
	if uselog:
		dat[dat<=2]=2
		dat=np.log2(dat)
	dat=scale(dat,axis=1,copy=False)
	# cluster
	dm=spatial.distance.pdist(dat,metric='euclidean')
	ll=cluster.hierarchy.single(dm)
	order=cluster.hierarchy.leaves_list(ll)

	newexp=hs.reorderbacteria(newexp,order)
	hs.addcommand(newexp,"clusterbacteria",params=params,replaceparams={'exp':exp})
	newexp.filters.append("cluster bacteria minreads=%d" % minreads)
	return newexp


def clustersamples(exp,minreads=0):
	"""
	cluster samples in an experiment according to similar behavior
	input:
	exp :Experiment
	minreads : int
		the minimal original number of reads per sample to keep it
	output:
	newexp : Experiment
		the filtered and clustered experiment
	"""
	params=locals()

	newexp=hs.filterorigreads(exp,minreads)
	# normalize each row (bacteria) to sum 1
	dat=copy.copy(newexp.data)
	dat=np.transpose(dat)
	dat[dat<=2]=2
	dat=np.log2(dat)
	# cluster
	dm=spatial.distance.pdist(dat,metric='braycurtis')
	ll=cluster.hierarchy.single(dm)
	order=cluster.hierarchy.leaves_list(ll)

	newexp=hs.reordersamples(newexp,order)
	hs.addcommand(newexp,"clustersamples",params=params,replaceparams={'exp':exp})
	newexp.filters.append("cluster samples minreads=%d" % minreads)
	return newexp


def sortsamples(exp,field,numeric=False,logit=True):
	"""
	sort samples according to field
	input:
	exp : Experiment
	field : string
		name of the field to sort by
	numeric : bool
		True for numeric values in field, false for text
	output:
	newexp : Experiment
		the sorted experiment
	"""
	params=locals()

	fvals=hs.getfieldvals(exp,field)
	if numeric:
		fvals=hs.tofloat(fvals)
	svals,sidx=hs.isort(fvals)
	newexp=hs.reordersamples(exp,sidx)

	if logit:
		hs.addcommand(newexp,"sortsamples",params=params,replaceparams={'exp':exp})
		newexp.filters.append('sorted samples by field %s' % field)
	return newexp


def sortbyfreq(expdat,field=False,value=False,exact=False,exclude=False,logscale=True):
	"""
	sort bacteria in experiment according to frequency
	sorting is performed based on a subset of samples (field/val/exact) and then
	all the experiment is sorted according to them
	input:
	expdat : Experiment
	field : string
		name of the field to filter samples for freq. sorting or False for all samples
	value : string
		value of samples to use for the freq. sorting
	exact : bool
		is the value exact or partial string
	exclude : bool
		True to sort on all samples except the field/value ones, False to sort only on field/value samples (default=False)
	logscale : bool
		True (default) to use log2 transform for frequencies before mean and sorting, False to use original values

	output:
	newexp : Experiment
		the experiment with bacteria sorted according to subgroup freq.
	"""
	params=locals()

	if field:
		texp=hs.filtersamples(expdat,field,value,exact=exact,exclude=exclude)
	else:
		texp=copy.deepcopy(expdat)
	if logscale:
		texp.data[texp.data<2]=2
		texp.data=np.log2(texp.data)
	meanvals=np.mean(texp.data,axis=1)
	svals,sidx=hs.isort(meanvals)

	newexp=hs.reorderbacteria(expdat,sidx)
	newexp.filters.append("sort by freq field=%s value=%s" % (field,value))
	hs.addcommand(newexp,"sortbyfreq",params=params,replaceparams={'expdat':expdat})
	return newexp

def sortbyvariance(expdat,field=False,value=False,exact=False,norm=False):
	"""
	sort bacteria by their variance
	sorting is performed based on a subset of samples (field/val/exact) and then
	all the experiment is sorted according to them
	input:
	expdat : Experiment
	field : string
		name of the field to filter samples for freq. sorting or False for all samples
	value : string
		value of samples to use for the freq. sorting
	exact : bool
		is the value exact or partial string
	norm : bool
		- False to sort by varinace, True to sort by variance/mean
	output:
	newexp : Experiment
		the experiment with bacteria sorted according to subgroup freq.
	"""
	params=locals()

	if field:
		texp=hs.filtersamples(expdat,field,value,exact=exact)
	else:
		texp=copy.deepcopy(expdat)

	svals=np.std(texp.data,axis=1)
	if norm:
		svals=svals/np.mean(texp.data,axis=1)
	svals,sidx=hs.isort(svals)

	newexp=hs.reorderbacteria(expdat,sidx)
	newexp.filters.append("sort by variance field=%s value=%s normalize=%s" % (field,value,norm))
	hs.addcommand(newexp,"sortbyvariance",params=params,replaceparams={'expdat':expdat})
	return newexp


def sortbygroupdiff(expdat,field,val1,val2):
	"""
	sort bacteria in the experiment by the difference in the mean between the 2 groups (val1,val2 in field)
	input:
	expdat
	field - the name of the field for the 2 groups
	val1,val2 - the values for the 2 groups

	output:
	newexp - the experiment with sorted bacteria
	"""
	params=locals()

	exp1=hs.filtersamples(expdat,field,val1,exact=True)
	exp2=hs.filtersamples(expdat,field,val2,exact=True)
	m1=np.mean(np.log2(exp1.data+2),axis=1)
	m2=np.mean(np.log2(exp2.data+2),axis=1)
	diff=(m1-m2)/(m1+m2+20)
	sv,si=hs.isort(diff)
	newexp=hs.reorderbacteria(expdat,si)
	newexp.filters.append("sort by group difference field=%s val1=%s val2=%s" % (field,val1,val2))
	hs.addcommand(newexp,"sortbygroupdiff",params=params,replaceparams={'expdat':expdat})
	return newexp


def sortcorrelation(expdat,method='all'):
	"""
	EXPERIMENTAL
	sort bacteria according to highest correlation/anti-correlation

	input:
	expdat
	method:
		pres - use correlation only on samples where present in one of the two sequnences
		all - use correlation on all samples (default)

	output:
	newexp - the experiment with bacteria sorted by correlation (each time next bacteria the most abs(corr) to the current bacteria)
	"""
	params=locals()

	cdat=copy.copy(expdat.data)
	cdat[cdat<=2]=2
	cdat=np.log2(cdat)
	cdat=scale(cdat,axis=1,copy=False,with_mean=False)
	hs.Debug(6,"Calculating correlation matrix")
	cmat=np.corrcoef(cdat)
	hs.Debug(6,"sorting bacteria")
	cmat=np.abs(cmat)
	cmat-=np.identity(len(expdat.seqs))
	maxpos=np.argmax(cmat)
	maxpos=np.unravel_index(maxpos,np.shape(cmat))
	order=[maxpos[0]]
	ubact=np.arange(len(expdat.seqs))
	ubact=np.delete(ubact,maxpos[0])
	maxpos=maxpos[0]
	while len(ubact)>0:
		cdat=cmat[ubact,maxpos]
		cdat=cdat.flatten()
		maxpos=np.argmax(cdat)
		order.append(ubact[maxpos])
		ubact=np.delete(ubact,maxpos)
	newexp=hs.reorderbacteria(expdat,order)
	newexp.filters.append("correlation sort")
	hs.addcommand(newexp,"sortcorrelation",params=params,replaceparams={'expdat':expdat})
	return newexp


############
# add sort by center of mass (for time/1d series)
###########


def sortbycentermass(expdat,field=False,numeric=True,uselog=True):
	"""
	sort bacteria in the experiment according to a 1d gradient by calculating the center of mass
	input:
	expdat
	field : string
		the name of the field to sort by or False to skip sorting
	numeric : bool
		True if the sort field is numeric (ignored if no sort field)
	uselog : bool
		True to log transform the data before mass center calculation
	output:
	newexp - the experiment with sorted bacteria
	"""
	params=locals()

	if field:
		newexp=hs.sortsamples(expdat,field,numeric=numeric)
	else:
		newexp=hs.copyexp(expdat)
	dat=newexp.data
	if uselog:
		dat[dat<1]=1
		dat=np.log2(dat)
	cm=[]
	multpos=np.arange(len(newexp.samples))
	for cseqind in range(len(newexp.seqs)):
		cm.append(np.dot(dat[cseqind,:],multpos)/np.sum(dat[cseqind,:]))
	sv,si=hs.isort(cm)
	newexp=hs.reorderbacteria(expdat,si)
	newexp.filters.append("sort by center of mass field=%s, uselog=%s" % (field,uselog))
	hs.addcommand(newexp,"sortbycentermass",params=params,replaceparams={'expdat':expdat})
	return newexp
