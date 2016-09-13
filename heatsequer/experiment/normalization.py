#!/usr/bin/env python


"""
heatsequer experiment normalization module
"""

# amnonscript

__version__ = "0.9"


import heatsequer as hs

import numpy as np


def normalizeprctile(expdat,percent=80):
	"""
	normalize reads per experiment so percentile (rather than mean) will be normalized
	used to reduce effect of outliers (compositionality correction)
	note normalization is done on the same set of bacteria for all samples
	input:
	expdat : Experiment
	percent : float
		the percentile to normalize (0-100)

	output:
	newexp : Experiment
		the new normalized experiment
	"""
	params=locals()

	# select the bacteria to use - don't want to include very low freq. bacteria
	newexp=hs.filterminreads(expdat,1*len(expdat.samples))

	percvals=np.percentile(newexp.data,percent,axis=0)
#	plt.figure()
#	plt.plot(percvals)
	percvals=percvals/np.mean(percvals)
	newexp=hs.copyexp(expdat)
	for idx,samp in enumerate(expdat.samples):
		newexp.data[:,idx]=newexp.data[:,idx]*percvals[idx]
	newexp.filters.append("normalize percentile %f" % percent)
	hs.addcommand(newexp,"normalizeprctile",params=params,replaceparams={'expdat':expdat})

	return newexp


def normalizereads(expdat,numreads=10000,fixorig=False,inplace=False):
	"""
	normalize the number of reads per sample to 10k
	input:
	expdat
	numreads - the number of reads to normalize to
	fixorig - True to fix origreads with the same ratio, False to keep as before
	inplace - true to replace orig experiment, false to create a new experiment

	output:
	newexp - the normalized experiment
	"""
	params=locals()

	if inplace:
		newexp=expdat
	else:
		newexp=hs.copyexp(expdat)

	for idx,csamp in enumerate(newexp.samples):
		totreads=np.sum(newexp.data[:,idx])
		if totreads==0:
			continue
		ratio=float(numreads)/totreads
		newexp.data[:,idx]=newexp.data[:,idx]*ratio
		if fixorig:
			hs.Debug(2,'fixing original frequencies')
			newexp.origreads[idx]=float(newexp.origreads[idx])/ratio
	newexp.filters.append("renormalized reads to sum %d" % numreads)
	hs.addcommand(newexp,"normalizereads",params=params,replaceparams={'expdat':expdat})
	return newexp


def normalizebyseqs(expdat,seqs,exclude=False,fixorig=True):
	"""
	normalize experiment by making the sum of frequencies in seqs constant in each sample
	input:
	expdat
	seqs - the sequences to use as the normalization factor (sum of the sequences)
	exclude - true to use all sequences except in seqs as the normalization factor, False to use seqs
	fixorig - True to modify the origreads field, false to leave it as it was
	"""
	params=locals()

	newexp=hs.copyexp(expdat)
	spos=[]
	for cseq in seqs:
		spos.append(expdat.seqdict[cseq])
	if exclude:
		spos=np.setdiff1d(np.arange(len(expdat.seqs)),spos)
	ssum=np.sum(expdat.data[spos,:],axis=0)+0.0
	ssum[ssum==0]=1
	frat=ssum/np.mean(ssum)
	for idx in range(len(expdat.samples)):
		newexp.data[:,idx]=newexp.data[:,idx]/frat[idx]
		if fixorig:
			newexp.origreads[idx]=newexp.origreads[idx]/frat[idx]
	filt='Normalize By Seqs '
	if len(spos)==1:
		filt+=newexp.tax[spos[0]]
	else:
		filt+=str(len(spos))
	if exclude:
		filt+=' Exclude'
	newexp.filters.append(filt)
	hs.addcommand(newexp,"normalizebyseqs",params=params,replaceparams={'expdat':expdat})
	return newexp


def toorigreads(expdat,inplace=False):
	"""
	convert the number of reads to absolute using the origreads field
	input:
	expdat
	inplace - True to replace current exp, false to create a new one

	output:
	newexp - each sample has origreads reads (instead of 10k)
	"""
	params=locals()

	if inplace:
		newexp=expdat
	else:
		newexp=hs.copyexp(expdat)

	colsum=hs.sum(newexp.data,axis=0)
	# for idx,csamp in enumerate(newexp.samples):
	# 	totreads=colsum[idx]
	# 	origreads=newexp.origreads[idx]
	# 	if totreads==0:
	# 		continue
	# 	ratio=float(origreads)/totreads
	# 	newexp.data[:,idx]=newexp.data[:,idx]*ratio
	newexp.data=hs.divvec(newexp.data,colsum/np.array(expdat.origreads))
	newexp.data=newexp.data.astype(int)

	newexp.filters.append("changed reads to origread value")
	hs.addcommand(newexp,"toorigreads",params=params,replaceparams={'expdat':expdat})
	return newexp


def subsample(expdat,numreads=10000,inplace=False):
	"""
	subsample (rarify) reads from all samples in an experiment
	input:
	expdat
	numreads - number of reads to subsample to
	inplace - true to replace current experiment

	output:
	newexp - the new subsampled experiment
	"""
	import biom

	params=locals()

	newexp=hs.filterorigreads(expdat,numreads,inplace)
	newexp=hs.toorigreads(newexp,inplace=True)

	table=biom.table.Table(newexp.data,newexp.seqs,newexp.samples)
	table=table.subsample(numreads,axis='observation')
	tids=table.ids(axis='sample')
	for idx,cid in enumerate(tids):
		if not cid==newexp.samples[idx]:
			print('problem with sample ids!!!!')
	newpos=[]
	for cseq in table.ids(axis='observation'):
		newpos.append(newexp.seqdict[cseq])
	newexp=hs.reorderbacteria(newexp,newpos,inplace=True)
	newexp.data=table.matrix_data.todense().A
	newexp=normalizereads(newexp,numreads=10000,inplace=True,fixorig=False)
	for cidx in range(len(newexp.samples)):
		newexp.origreads[cidx]=numreads
	newexp=updateorigreads(newexp)
	newexp.filters.append("subsample to %d" % numreads)
	hs.addcommand(newexp,"subsample",params=params,replaceparams={'expdat':expdat})
	return newexp


def updateorigreads(expdat,logit=True):
	params=locals()

	for idx,csamp in enumerate(expdat.samples):
		expdat.smap[csamp]['origReads']=expdat.origreads[idx]
	if logit:
		expdat.filters.append("Update orig reads")
		hs.addcommand(expdat,"updateorigreads",params=params,replaceparams={})
	return expdat



def findcandidatecompositional(expdat,fracappear=0.1,minlevel=3000):
	"""
	find the sequences that may affect compositional data
	for georg nose analysis
	input:
	expdat : Experiment
	fracappear : float
		fraction of samples where the otu is above the minlevel
	minlevel : float
		the minimal number of reads in a sample to count it

	output:
	oseqs : list of sequences
		a list of candidate sequences
	"""

	dat=np.sum(expdat.data>=minlevel,1)
	iseqs=np.where(dat>=fracappear*len(expdat.samples))[0]
	oseqs=[]
	for cseq in iseqs:
		oseqs.append(expdat.seqs[cseq])
	return oseqs
