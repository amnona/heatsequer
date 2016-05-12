#!/usr/bin/env python


"""
heatsequer analysis module
"""

# amnonscript

__version__ = "0.9"

import heatsequer as hs

import numpy as np
from scipy import stats
from scipy import spatial
from collections import defaultdict
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


def getdiffsigall(expdat,field,val1,val2=False,numperm=1000,maxfval=0.1,nofreqpres=False):
	"""
	Get the differentially abundant bacteria (using getdiffsig) using all methods possible.
	Sort results according to combined effect order
	input:
	see getdiffsig()

	output:
	newexp - the new experiment with bacteria significantly differentiating the 2 groups by at least 1 method
	"""
	params=locals()

	# do the 4 tests
	if nofreqpres:
		methods=['mean','binary','ranksum']
	else:
		methods=['mean','binary','ranksum','freqpres']
	res=[]
	for cmethod in methods:
		res.append(hs.getdiffsig(expdat,field=field,val1=val1,val2=val2,method=cmethod,numperm=numperm,maxfval=maxfval))

	# combine the results from the different tests
	keep=[]
	keeporder=[]
	for cidx,cseq in enumerate(expdat.seqs):
		pos=[]
		edir=[]
		for cres in res:
			if cseq in cres.seqdict:
				pos.append(float(cres.seqdict[cseq])/len(cres.seqs))
				edir.append(np.sign(cres.odif[cres.seqdict[cseq]]))
		if len(pos)>0:
			keep.append(cidx)
			keeporder.append(np.mean(pos)*np.mean(edir))
	keep=np.array(keep)
	if len(keep)>0:
		si=np.argsort(keeporder)
		newexp=hs.reorderbacteria(expdat,keep[si])
		newexp.diffdir=hs.reorder(np.sign(keeporder),si)
		newexp.filters.append('differential expression (all) in %s between %s and %s' % (field,val1,val2))
		hs.addcommand(newexp,"getdiffsigall",params=params,replaceparams={'expdat':expdat})
		return newexp
	else:
		hs.Debug(6,'No bacteria found')
		newexp=hs.reorderbacteria(expdat,[])
		return newexp


def getdiffsig(expdat,field,val1,val2=False,method='mean',numperm=1000,maxfval=0.1):
	"""
	test the differential expression between 2 groups (val1 and val2 in field field)
	for bacteria that have a high difference.
	input:
	expdat
	field - the field for the 2 categories
	val1 - values for the first group
	val2 - value for the second group or false to compare to all other
	method - the test to compare the 2 groups:
		mean - absolute difference in mean frequency
		binary - abs diff in binary presence/absence
		ranksum - abs diff in rank order (to ignore outliers)
		freqpres - abs diff in frequency only in samples where bacteria is present
	numperm - number of random permutations to run
	maxfval - the maximal f-value (FDR) for a bacteria to keep

	output:
	newexp - the experiment with only significant (FDR<=maxfval) difference, sorted according to difference
	"""
	params=locals()

	minthresh=2
	exp1=hs.filtersamples(expdat,field,val1,exact=True)
	if val2:
		exp2=hs.filtersamples(expdat,field,val2,exact=True)
	else:
		exp2=hs.filtersamples(expdat,field,val1,exact=True,exclude=True)
	cexp=hs.joinexperiments(exp1,exp2)
	len1=len(exp1.samples)
	len2=len(exp2.samples)
	dat=cexp.data
	dat[dat<minthresh]=minthresh
	dat=np.log2(dat)
	numseqs=len(cexp.seqs)

	eps=0.000001

	if method=='mean':
		pass
	elif method=='ranksum':
		for idx in range(numseqs):
			dat[idx,:]=stats.rankdata(dat[idx,:])
	elif method=='binary':
		dat=(dat>np.log2(minthresh))
	elif method=='freqpres':
		dat[dat<=minthresh]=np.nan
		# remove sequences that don't have >thresh in one group
		keeplist=[]
		for cbact in range(len(cexp.seqs)):
			keepit=True
			if np.sum(np.isnan(dat[cbact,0:len1]))==len1:
				keepit=False
			if np.sum(np.isnan(dat[cbact,len1:]))==len2:
				keepit=False
			if keepit:
				keeplist.append(cbact)
		hs.Debug(6,'keeping %d bacteria (out of %d) for freqpres analysis' % (len(keeplist),len(cexp.seqs)))
		cexp=hs.reorderbacteria(cexp,keeplist)
		dat=cexp.data
		dat[dat<minthresh]=minthresh
		dat=np.log2(dat)
		dat[dat<=minthresh]=np.nan
	else:
		hs.Debug(9,"Method not supported!",method)
		return


	m1=np.nanmean(dat[:,0:len1],axis=1)
	m2=np.nanmean(dat[:,len1:],axis=1)
	odif=(m1-m2)/(m1+m2+eps)
	odif[np.isnan(odif)]=0
	alldif=np.zeros([len(cexp.sids),numperm])
	for x in range(numperm):
		rp=np.random.permutation(len1+len2)
		m1=np.nanmean(dat[:,rp[0:len1]],axis=1)
		m2=np.nanmean(dat[:,rp[len1:]],axis=1)
		diff=(m1-m2)/(m1+m2+eps)
#		diff[np.isnan(diff)]=0
		alldif[:,x]=diff
	pval=[]
	for crow in range(len(odif)):
		cdat=alldif[crow,:]
		cdat=cdat[np.logical_not(np.isnan(cdat))]
		cnumperm=len(cdat)
		if cnumperm==0:
			pval.append(1)
			continue
		cpval=float(np.sum(np.abs(cdat)>=np.abs(odif[crow])))/cnumperm
		# need to remember we only know as much as the number of permutations - so add 1 as upper bound for fdr
		cpval=min(cpval+(1.0/cnumperm),1)
		pval.append(cpval)
	# NOTE: maybe use better fdr (this one not tested...)
	fval=hs.fdr(pval)
	keep=np.where(np.array(fval)<=maxfval)
	seqlist=[]
	for cidx in keep[0]:
		seqlist.append(cexp.seqs[cidx])

	# do we need this or is reorder enough?
	newexp=hs.filterseqs(expdat,seqlist,logit=False)
	odif=odif[keep[0]]
	sv,si=hs.isort(odif)
	hs.Debug(6,'method %s. number of higher in %s : %d. number of higher in %s : %d. total %d' % (method,val1,np.sum(odif>0),val2,np.sum(odif<0),len(odif)))
	newexp=hs.reorderbacteria(newexp,si)
	newexp.odif=sv
	hs.addcommand(newexp,"getdiffsigall",params=params,replaceparams={'expdat':expdat})
	newexp.filters.append('differential expression (%s) in %s between %s and %s' % (method,field,val1,val2))
	return newexp


def getsigcorr(expdat,field=False,numeric=True,numperm=1000,maxfval=0.1,method='ranksum'):
	"""
	find bacteria significantly associated (correlated) with a mapping file sortable field
	input:
	expdat
	field : string
		the field to sort experiment by, or false to use current order
	numeric : bool
		True if values in field are numeric, False if not (for ascii sort)
	method - the test to compare the 2 groups:
		'ranksum' - ranksum correlation
	numperm - number of random permutations to run
	maxfval - the maximal f-value (FDR) for a bacteria to keep

	output:
	newexp - the experiment with only significant (FDR<=maxfval) difference, sorted according to difference
	"""
	params=locals()

	minthresh=2
	if field:
		cexp=hs.sortsamples(expdat,field,numeric,logit=False)

	numseqs=len(cexp.seqs)
	numsamps=len(cexp.samples)

	dat=cexp.data
	dat[dat<minthresh]=minthresh
	dat=np.log2(dat)
	dat=dat-np.mean(dat,axis=1,keepdims=True)

	eps=0.000001

	if method=='ranksum':
		compareto=np.arange(numsamps)
		compareto=compareto-np.mean(compareto,keepdims=True)
#		compareto=compareto-np.mean(compareto,axis=1,keepdims=True)
		# calculate the rank
		for idx in range(numseqs):
			dat[idx,:]=stats.rankdata(dat[idx,:])
	else:
		hs.Debug(9,"Method not supported!",method)
		return
	odif=np.zeros(numseqs)
	for idx in range(numseqs):
		odif[idx]=np.sum(dat[idx,:]*compareto)
#		odif[idx]=spatial.distance.correlation(dat[idx,:],compareto)

	alldif=np.zeros([len(cexp.sids),numperm])
	for x in range(numperm):
		rp=np.random.permutation(numsamps)
		diff=np.zeros(numseqs)
		for idx in range(numseqs):
#			diff[idx]=spatial.distance.correlation(dat[idx,rp],compareto)
			diff[idx]=np.sum(dat[idx,rp]*compareto)
		alldif[:,x]=diff

	pval=[]
	for crow in range(len(odif)):
		cdat=alldif[crow,:]
		cdat=cdat[np.logical_not(np.isnan(cdat))]
		cnumperm=len(cdat)
		if cnumperm==0:
			pval.append(1)
			continue
		cpval=float(np.sum(np.abs(cdat)>=np.abs(odif[crow])))/cnumperm
		# need to remember we only know as much as the number of permutations - so add 1 as upper bound for fdr
		cpval=min(cpval+(1.0/cnumperm),1)
		pval.append(cpval)
	# NOTE: maybe use better fdr (this one not tested...)
	fval=hs.fdr(pval)
	keep=np.where(np.array(fval)<=maxfval)
	seqlist=[]
	for cidx in keep[0]:
		seqlist.append(cexp.seqs[cidx])

	# do we need this or is reorder enough?
	newexp=hs.filterseqs(expdat,seqlist,logit=False)
	odif=odif[keep[0]]
	sv,si=hs.isort(odif)
	newexp=hs.reorderbacteria(newexp,si)
	hs.addcommand(newexp,"getdiffsigall",params=params,replaceparams={'expdat':expdat})
	newexp.filters.append('significant correlation %s for field %s' % (method,field))
	return newexp


def bicluster(expdat,numiter=5,startb=False,starts=False,method='binary',sampkeep=0.5,bactkeep=0.25,justcount=False,numruns=1):
	"""
	EXPERIMENTAL
	cluster bacteria and samples from subgroup
	input:
	expdat
	numiter - number of iterations to run the biclustering
	startb - start list of bacteria [acgt] of False for random
	method: the method for choosing which sample/bacteria to keep. options areL
		'zscore'
		'ranksum'
		'binary' - only one working currently!!!
	sampkeep - the minimal fraction of bacteria to be present in a sample in order to keep the sample (for binary) or 0 for random
	bactkeep - the minimal difference in the number of samples a bacteria apprears in order to keep the bacteria (for binary) or 0 for random
	justcount - True to not reorder the experiment - just the bacteria & samples (to run faster)
	numruns - number of times to run

	output:
	newexp - the reordered experiment
	seqs - the sequences in the cluster
	samples - the samples in the cluster (position in experiment)
	"""

	dat=copy.copy(expdat.data)
	if method=='zscore':
		dat[dat<20]=20
		dat=np.log2(dat)
#		dat=(dat>1)
	elif method=='binary':
		dat=(dat>1)
	bdat=dat
#	bdat=scale(dat,axis=1,copy=True)
	nbact=np.size(dat,0)
	nsamp=np.size(dat,1)
	allsamp=np.arange(nsamp)
	allbact=np.arange(nbact)


	allseqs=[]
	allsamples=[]
	bthresh=0.25
	sthresh=0
	for crun in range(numruns):
		if bactkeep==0:
			bactkeep=np.random.uniform(0.1,0.5)
			hs.Debug(6,"bactkeep %f" % bactkeep)
		if sampkeep==0:
			sampkeep=np.random.uniform(0.25,0.75)
			hs.Debug(6,"sampkeep %f" % sampkeep)
		if startb:
			ubact=[]
			for cbact in startb:
				ubact.append(expdat.seqdict[cbact])
		else:
			ubact=[np.random.randint(nbact)]
		if starts:
			usamp=starts
		else:
			usamp=np.arange(nsamp)
		for citer in range(numiter):
			if method=='zscore':
				# find samples
				meanin=np.mean(bdat[ubact,:],axis=0)
	#			print(meanin[0:10])
				sdiff=meanin-np.mean(np.mean(bdat[ubact,:]))
	#			print(sdiff[0:10])
				if len(ubact)>1:
					usamp=allsamp[sdiff>sthresh*np.std(np.mean(bdat[ubact,:],axis=0))]
				else:
					usamp=allsamp[sdiff>sthresh*np.std(bdat[ubact,:])]
				print("num samples %d" % len(usamp))

				meanin=np.mean(dat[:,usamp],axis=1)
				nusamp=np.setdiff1d(allsamp,usamp)
				sdiff=meanin-np.mean(dat[:,nusamp],axis=1)
	#			sdiff=meanin-np.mean(np.mean(bdat[:,usamp]))
				if len(usamp)>1:
	#				ubact=allbact[sdiff>bthresh*np.std(np.mean(bdat,axis=1))]
					ubact=allbact[sdiff>bthresh]
				else:
	#				ubact=allbact[sdiff>bthresh*np.std(bdat)]
					ubact=allbact[sdiff>bthresh]
				print("num bacteria %d" % len(ubact))

			elif method=='binary':
				# find samples
				meanin=np.mean(bdat[ubact,:],axis=0)
				sdiff=meanin
				if len(ubact)>1:
					usamp=allsamp[sdiff>=sampkeep]
				else:
					usamp=allsamp[sdiff>=sampkeep]
				print("num samples %d" % len(usamp))

				meanin=np.mean(dat[:,usamp],axis=1)
				nusamp=np.setdiff1d(allsamp,usamp)
				sdiff=meanin-np.mean(dat[:,nusamp],axis=1)
				if len(usamp)>1:
	#				ubact=allbact[sdiff>bthresh*np.std(np.mean(bdat,axis=1))]
					ubact=allbact[sdiff>=bactkeep]
				else:
	#				ubact=allbact[sdiff>bthresh*np.std(bdat)]
					ubact=allbact[sdiff>=bactkeep]
				print("num bacteria %d" % len(ubact))

			elif method=='ranksum':
				nubact=np.setdiff1d(allbact,ubact)
				keepsamp=[]
				apv=[]
				astat=[]
				for idx,csamp in enumerate(expdat.samples):
					g1=bdat[ubact,idx]
					g2=bdat[nubact,idx]
					if len(g1)>1:
						g1=np.squeeze(g1)
					if len(g2)>1:
						g2=np.squeeze(g2)
					stat,pv=stats.mannwhitneyu(g2,g1)
					apv.append(pv)
					astat.append(stat)
					if pv<0.05:
						keepsamp.append(idx)
				# figure()
				# hist(apv,100)
				# show()
				# figure()
				# hist(astat,100)
				# show()
				usamp=keepsamp
				print('number of samples: %d' % len(usamp))
				nusamp=np.setdiff1d(allsamp,usamp)
				keepbact=[]
				for idx,cbact in enumerate(expdat.sids):
					g1=np.squeeze(bdat[idx,usamp])
					g2=np.squeeze(bdat[idx,nusamp])
					try:
						stat,pv=stats.mannwhitneyu(g2,g1)
						if pv<0.001:
							keepbact.append(idx)
					except:
						pass
				ubact=keepbact
				print('number of bacteria: %d' % len(ubact))

			else:
				hs.Debug(9,"biclustering method %s not supported")
				return

		x=np.setdiff1d(allsamp,usamp)
		sampo=np.concatenate((usamp,x))
		bacto=np.concatenate((ubact,np.setdiff1d(allbact,ubact)))

		seqs=[]
		for cbact in ubact:
			seqs.append(expdat.seqs[cbact])
		samples=[]
		for csamp in usamp:
			samples.append(csamp)

		if not justcount:
			newexp=hs.reordersamples(expdat,sampo)
			newexp=hs.reorderbacteria(newexp,bacto,inplace=True)
			newexp.filters.append('biclustering')
		else:
			newexp=False

		allseqs.append(seqs)
		allsamples.append(samples)
	return newexp,allseqs,allsamples


def testbicluster(expdat,numiter=5,numruns=100):
	"""
	show the sizes of clusters in data, random and random normalized
	"""
	plt.figure()
	cc,seqs,samps=bicluster(expdat,method='binary',justcount=True,numruns=numruns,numiter=numiter)
	for idx,cseqs in enumerate(seqs):
		plt.plot(len(samps[idx]),len(cseqs),'xr')
	rp=randomizeexp(expdat,normalize=False)
	cc,seqs,samps=bicluster(rp,method='binary',justcount=True,numruns=numruns,numiter=numiter)
	for idx,cseqs in enumerate(seqs):
		plt.plot(len(samps[idx]),len(cseqs),'xk')
	rp=randomizeexp(expdat,normalize=True)
	cc,seqs,samps=bicluster(rp,method='binary',justcount=True,numruns=numruns,numiter=numiter)
	for idx,cseqs in enumerate(seqs):
		plt.plot(len(samps[idx]),len(cseqs),'xk')
	plt.xlabel('# Samples in cluster')
	plt.ylabel('# Bacteria in cluster')


def getbigfreqsource(expdat,seqdb):
	"""
	EXPERIMENTAL
	get the most common automatic db annotation for each sequence in the experiment

	input:
	expdat : Experiment
		the experiment with sequences to annotate
	seqdb : BactDB
		the automatic bacterial database (from bactdb.load())

	output:
	newep : Experiment
		the annotated experiment (source instead of taxonomy for each bacteria)
	"""
	hs.Debug(6,"Getting db info for all sequences")
	dat=hs.bactdb.GetDBSource(seqdb,expdat.seqs)

	# get the metadata info about each sample
	hs.Debug(6,"Getting metadata info")
	sampinfo=hs.bactdb.GetAllSampleInfo(seqdb,['ENV_MATTER','HOST_TAXID'])
	hs.Debug(6,"Got info. Processing")
	# count number of appearance in each category for each bacteria
	numcat=np.zeros([len(expdat.seqs),0])
	cats=[]
	for k,v in sampinfo.items():
		if len(v)>1:
			numcat=np.hstack([numcat,(np.sum(dat[:,v]>0.001,axis=1,keepdims=True)+0.0)/max(10,len(v))])
		else:
			numcat=np.hstack([numcat,((dat[:,v]>0.001)+0.0)/10])
		cats.append(k)

	newexp=hs.copyexp(expdat)
	for idx,cbact in enumerate(expdat.seqs):
		si=np.argmax(numcat[idx,:])
		sv=numcat[idx,si]
		mcat=cats[si]
#		newexp.tax[idx]=mcat+'-%f' % sv
		newexp.tax[idx]=mcat

	return newexp


def getexpdbsources(expdat,seqdb=False):
	'''
	EXPERIMENTAL
	predict the probable source for each bacteria
	based on the bactdb automatic bacteria database,
	by each time taking the sample with the highest number of common bacteria with the experiment
	and defining that sample as the source for these bacteria.
	input:
	expdat : Experiment
	seqdb : the bactdb automatic database ( from bactdb.load() )
	output:
	newexp : Experiment
		same as expdat, but taxonomy now modified to include predicted source for each bacteria
	'''
	params=locals()

	if not seqdb:
		if not expdat.seqdb:
			hs.Debug(9,'No sequence database loaded')
			return
		else:
			seqdb=expdat.seqdb

	dat=hs.bactdb.GetDBSource(seqdb,expdat.seqs)

	newexp=hs.copyexp(expdat)

	THRESH=0.001
	used=np.arange(np.size(dat,0))
	done=False
	while not done:
		npsamp=np.sum(dat[used,:]>=THRESH,axis=0)
		pos=np.argmax(npsamp)
		print('sample is %d, size is %d' % (pos,npsamp[pos]))
		sid,sname=hs.bactdb.GetSampleStudy(expdat.seqdb,pos+1)
		print('studyid %d name %s' % (sid,sname))
		ubact=np.where(dat[used,pos]>=THRESH)[0]
		for cpos in ubact:
			newexp.tax[used[cpos]]+=sname
		used=np.setdiff1d(used,used[ubact])
		if len(used)<10:
			print('got them all')
			done=True
		if len(ubact)<=2:
			print('piti')
			done=True
		if npsamp[pos]<=2:
			print('pata')
			done=True
	newexp.filters.append("Assign sources from bactdb (max sample binary overlap)")
	hs.addcommand(newexp,"getexpdbsources",params=params,replaceparams={'expdat':expdat})
	return newexp



def BaysZeroClassifyTest(oexpdat,field,val1,val2=False,n_folds=10,istreeexpand=False,randseed=False,numiter=1,method='kfold'):
	"""
	Test the baysian zero inflated classifier by doing n_fold cross validation on a given dataset
	input:
	expdat
	field - the field to use for classification
	val1 - the value of group1 in field
	val2 - value of group2 in field, or False to use all non val1 as group2
	n_folds - number of groups to divide to for crossvalidation
	istreeexpand - True if we want to use the tree shrinkage on each training set (if the exp is from addsubtrees)
	randseed - if non False, use the specified random seed for the test/validation division
	numiter - the number of times to run the cross validation
	method : string
		- 'kfold' - k-fold cross validation (using n_folds variable) (default)
		- 'loo' - leave one out cross validation

	output:
	auc - the auc of each iteration
	"""
	if randseed:
		np.random.seed(randseed)

	# we want only to keep the 2 classes
	if val2:
		oexpdat=hs.filtersamples(oexpdat,field,[val1,val2])

	# remove the validation samples and keep them in a different experiment valexp
	ind1=hs.findsamples(oexpdat,field,val1)
	types=np.zeros(len(oexpdat.samples))
	types=types>10
	types[ind1]=True
	aucres=[]
	for citer in range(numiter):
		if method=='kfold':
			rs = sklearn.cross_validation.StratifiedKFold(types, n_folds=n_folds,shuffle=True)
			for trainidx,testidx in rs:
				valexp=hs.reordersamples(oexpdat,testidx)
				expdat=hs.reordersamples(oexpdat,trainidx)
				# classify
				lrscores=hs.BayZeroClassify(expdat,valexp,field,val1,val2,istreeexpand)

				# prepare the correct answer list
				typeis1=[]
				for vsamp in valexp.samples:
					if valexp.smap[vsamp][field]==val1:
						vtype=1
					else:
						vtype=2
					typeis1.append(vtype==1)
				cauc=sklearn.metrics.roc_auc_score(typeis1, lrscores, average='macro', sample_weight=None)
				hs.Debug(4,"auc=%f" % cauc)
				aucres.append(cauc)
		elif method=='loo':
			rs = sklearn.cross_validation.LeaveOneOut(len(types))
			predlist=[]
			groundtruthlist=[]
			for trainidx,testidx in rs:
				valexp=hs.reordersamples(oexpdat,testidx)
				expdat=hs.reordersamples(oexpdat,trainidx)
				# classify
				lrscores=hs.BayZeroClassify(expdat,valexp,field,val1,val2,istreeexpand)

				# prepare the correct answer list
				for vsamp in valexp.samples:
					if valexp.smap[vsamp][field]==val1:
						vtype=1
					else:
						vtype=2
					groundtruthlist.append(vtype==1)
#				predlist.append(lscores[0])
				predlist.append(lrscores[0])
			cauc=sklearn.metrics.roc_auc_score(typeis1, lrscores, average='macro', sample_weight=None)
			hs.Debug(4,"auc=%f" % cauc)
			aucres.append(cauc)
		else:
			hs.Debug(9,'Crossvalidation Method %s not supported' % method)
			return

	hs.Debug(8,"mean is : %f" % np.mean(aucres))
	hs.Debug(8,"s.err is : %f" % (np.std(aucres)/np.sqrt(len(aucres))))
	return(aucres)


def BayZeroClassify(expdat,valexp,field,val1,val2=False,istreeexpand=False):
	"""
	Do the Zero inflated Naive Bayes Classifier
	Does a p-value based on presence/absence if bacteria is 0 in sample, otherwise do non-parametric permutaion test p-value
	combine p-values for all bacteria as if independent and look at log ratio of 2 categories as prediction
	input:
	expdat - the training set
	valexp - the validation set (to be classified)
	field - the field to use for classification
	val1 - the value of group1 in field
	val2 - value of group2 in field, or False to use all non val1 as group2
	istreeexpand - True if we want to use the tree shrinkage on each training set (if the exp is from addsubtrees)

	output:
	pred - the log2(ratio) prediction score for each sample in the validation experiment (>0 means from val1, <0 means from val2)
	"""

	# if val2 is not empty, keep only samples with val1 or val2
	if val2:
		expdat=hs.filtersamples(expdat,field,[val1,val2])
	ind1=hs.findsamples(expdat,field,val1)
	types=np.zeros(len(expdat.samples))
	types=types>10
	types[ind1]=True

	# if an expanded tree, keep the best subtrees
	if istreeexpand:
		expdat=hs.keeptreebest(expdat,field,val1,val2)
		valexp=hs.filterseqs(valexp,expdat.seqs)

	# prepare the claissifier
	g1ind=hs.findsamples(expdat,field,val1)
	if val2:
		g2ind=hs.findsamples(expdat,field,val2)
	else:
		g2ind=hs.findsamples(expdat,field,val1,exclude=True)

	tot1=len(g1ind)
	tot2=len(g2ind)
	dat=expdat.data
	zero1=np.sum(dat[:,g1ind]==0,axis=1)
	zero2=np.sum(dat[:,g2ind]==0,axis=1)

	# p value for getting a 0 in each sample type
	# we do max(zero1,1) to remove effect of sampling error
	MINTHRESH=1
	pmissing1=np.divide((np.maximum(MINTHRESH,zero1)+0.0),tot1)
	pmissing2=np.divide((np.maximum(MINTHRESH,zero2)+0.0),tot2)

	ppres1=np.divide((np.maximum(MINTHRESH,tot1-zero1)+0.0),tot1)
	ppres2=np.divide((np.maximum(MINTHRESH,tot2-zero2)+0.0),tot2)
	# and the log ratio of proability 1 to probability 2
	lograt0=np.log2(pmissing1/pmissing2)
	logratn0=np.log2(ppres1/ppres2)

	# the prediction log ratio scores
	lrscores=[]
	for vidx,vsamp in enumerate(valexp.samples):
		hs.Debug(2,"Classifying sample %s" % vsamp)
		cvdat=valexp.data[:,vidx]
		vzero=np.where(cvdat==0)[0]
		crat0=np.sum(lograt0[vzero])
		vnzero=np.where(cvdat>0)[0]
		cratn0=np.sum(logratn0[vnzero])
		# need to choose bigger or smaller (direction of test)
		# we test both and take the more extreme p-value
		# the probability to be bigger
		ratnz=[]
		for cnzpos in vnzero:
			allz=np.where(dat[cnzpos,:]>0)[0]
			nz1=np.intersect1d(allz,g1ind)
			if len(nz1)<5:
				continue
			nz2=np.intersect1d(allz,g2ind)
			if len(nz2)<5:
				continue
			d1=dat[cnzpos,nz1]
			d2=dat[cnzpos,nz2]
			p1b=(0.0+max(1,np.sum(d1>=cvdat[cnzpos])))/len(d1)
			p2b=(0.0+max(1,np.sum(d2>=cvdat[cnzpos])))/len(d2)
			ratb=np.log2(p1b/p2b)
			p1s=(0.0+max(1,np.sum(d1<=cvdat[cnzpos])))/len(d1)
			p2s=(0.0+max(1,np.sum(d2<=cvdat[cnzpos])))/len(d2)
			rats=np.log2(p1s/p2s)

			if np.abs(ratb)>=np.abs(rats):
				ratnz.append(ratb)
			else:
				ratnz.append(rats)
		cratfreq=np.sum(ratnz)
		totratio=crat0+cratfreq+cratn0+np.log2((tot1+0.0)/tot2)
		lrscores.append(totratio)
		hs.Debug(2,"LRScore %f, zscore %f, nzscore %f, freqscore %f" % (totratio,crat0,cratn0,cratfreq))
	hs.Debug(3,"Finished classifying")
	return lrscores



def keeptreebest(expdat,field,val1,val2,method="meandif"):
	"""
	EXPERIMENTAL
	keep only the best combinations wrt given criteria
	use after addsubtrees()

	input:
	expdat - after addsubtrees
	field - the field to use for comparison
	val1 - value for group1
	val2 - value for group2 or False for all except group1
	method:
		meandif - keep the largest mean difference between groups / total mean
	"""
	params=locals()

	pos1=hs.findsamples(expdat,field,val1)
	if val2:
		pos2=hs.findsamples(expdat,field,val2)
	else:
		pos2=hs.findsamples(expdat,field,val1,exclude=True)
	allpos=list(set(pos1+pos2))

	mean1=np.mean(expdat.data[:,pos1],axis=1)
	mean2=np.mean(expdat.data[:,pos2],axis=1)
	meanall=np.mean(expdat.data[:,allpos],axis=1)
	minval=1.0/len(allpos)
	meanall[meanall<1.0/minval]=minval
	difval=np.abs((mean1-mean2)/meanall)

	si=np.argsort(difval)
	si=si[::-1]

	dontuse={}
	keep=[]
	for cidx in si:
		cseq=expdat.seqs[cidx]
		poss=cseq.split(',')
		if len(poss)==1:
			if not str(expdat.seqdict[cseq]) in dontuse:
				keep.append(cidx)
				dontuse[str(expdat.seqdict[cseq])]=True
			continue
		keepit=True
		for cpos in poss:
			if cpos=='':
				continue
			if cpos in dontuse:
				keepit=False
				break
		if keepit:
			keep.append(cidx)
			for cpos in poss:
				dontuse[cpos]=True
	newexp=hs.reorderbacteria(expdat,keep)
	newexp.filters.append("keeptreebest field %s val1 %s val2 %s" % (field,val1,str(val2)))
	hs.addcommand(newexp,"keeptreebest",params=params,replaceparams={'expdat':expdat})
	return newexp


def randomizeexp(expdat,normalize=False):
	"""
	randomly permute each bacteria in the experiment indepenedently (permute the samples where it appears)
	input:
	expdat
	normalize - True to renormalize each sample to constant sum, False to not normalize

	output:
	newexp - the permuted experiment
	"""
	params=locals()

	newexp=hs.copyexp(expdat)
	numsamps=len(newexp.samples)
	for idx,cseq in enumerate(newexp.seqs):
		rp=np.random.permutation(numsamps)
		newexp.data[idx,:]=newexp.data[idx,rp]
	if normalize:
		newexp=hs.normalizereads(newexp,inplace=True,fixorig=False)

	newexp.filters.append("RANDOMIZED!!! normalize = %s" % normalize)
	hs.addcommand(newexp,"randomizeexp",params=params,replaceparams={'expdat':expdat})
	return newexp


def testmdenrichment(expdat,samples,field,numeric=False):
	"""
	test for enrichment in a subset of samples of the experiment for metadata field field
	input:
	expdat
	samples - the samples (positions) for the enrichment testing
	field - the field to test
	numeric - True if the field is numeric (test mean)
	"""

	vals=hs.getfieldvals(expdat,field)
	numsamps=len(vals)
	numgroup=len(samples)
	uvals=list(set(vals))
	gmap=defaultdict(list)
	for idx,cval in enumerate(vals):
		gmap[cval].append(idx)

	pv={}
	for cval in uvals:
		glen=float(len(gmap[cval]))
		numin=float(len(np.intersect1d(samples,gmap[cval])))
		pnull=glen/numsamps
		p1=stats.binom.cdf(numin,numgroup,pnull)
		p2=stats.binom.cdf(numgroup-numin,numgroup,1-pnull)
		p=min(p1,p2)
		pv[cval]={}
		pv[cval]['pval']=p
		pv[cval]['observed']=numin
		pv[cval]['expected']=pnull*numgroup

#		if p<0.05:
#			print("cval %s numin %f groupsize %d pnull %f p1 %f p2 %f" % (cval,numin,numgroup,pnull,p1,p2))

	return pv


def testmdenrichmentall(expdat,samples,maxpv=0.001,fdr=0.05):
	"""
	test enrichment in all metadata fields/values
	input:
	expdat
	samples - a list of samples to test (positions, not sample names)
	fdr - the false discovery rate in order to show a category

	output:
	upv - a list of dict of pvalues for significant fields/values ('pval','expected','observed','field','val')
	"""

	justp=[]
	allpv=[]
	for cfield in expdat.fields:
		vals=hs.getfieldvals(expdat,cfield)
		uvals=list(set(vals))
		if len(uvals)>10:
			continue
		pv=testmdenrichment(expdat,samples,cfield)
		for k,v in pv.items():
			justp.append(v['pval'])
			v['field']=cfield
			v['val']=k
			allpv.append(v)
#			if v['pval']<=maxpv:
#				print("field %s, val %s, pv %f (observed %d, expected %f)" % (cfield,k,v['pval'],v['observed'],v['expected']))

	# do the fdr if needed
	if fdr:
		fval=hs.fdr(justp)
		keep=np.where(np.array(fval)<=fdr)
		keep=keep[0]
	else:
		keep=np.arange(len(justp))

	if len(keep)==0:
		hs.Debug(6,'No significant cateogries found')

	upv=[]
	for ckeep in keep:
		upv.append(allpv[ckeep])
		hs.Debug(6,allpv[ckeep])

	upv=sortenrichment(upv)
	return upv



def sortenrichment(enrich,method='bidirectional',epsilon=2):
	"""
	sort an enrichment list (with 'observed', 'expected' dict values) according to effect size
	the effect size is abs(log(obs/(expected+EPS)))
	input:
	enrich - a list of dict with 'observed' abd 'expected' keys (i.e.e from testmdenrichmentall)
	method:
		bidirectional - use abs(log(o+EPS)/log(E+EPS))
		single - use log(o)/log(E)
		val - use o
	epsilon - the value used to reduce effect of low counts (a+eps)/(b+eps)

	output:
	newenrich - the sorted list
	"""

	# get the effect size
	effects=[]
	for citem in enrich:
		if method=='bidirectional':
			lograt=np.log2((citem['observed']+epsilon)/(epsilon+citem['expected']))
			effects.append(np.abs(lograt))
		elif method=='single':
			lograt=np.log2((citem['observed'])/np.log2(citem['expected'])+epsilon)
			effects.append(lograt)
		elif method=='val':
			effects.append(citem['observed'])
		else:
			hs.Debug('method %s not supported' % method)
	si=np.argsort(effects)
	newenrich=hs.reorder(enrich,si[::-1])
	return newenrich


def testenrichment(data,group,method='binary',fdr=0.05,twosided=False,printit=True):
	"""
	test for enrichment for samples in groupind in the dict of arrays data
	input:
	data - a dict (by category value) of numpy arrays (each of length totseqs) of the value of each sequence
	group - the indices of the group elements
	method - the test to apply:
		'binary' - presence/abscence
		'ranksum' - not implemented yet
	fdr - the false discovery rate value or false for no fdr
	twosided - True to test both lower and higher, False to test just higher in group
	printit - True to print the significant, False to not print

	output:
	plist - a list of dict entries ('pval','observed','expected','name')
	"""

	grouplen=len(group)
	allpv=[]
	justp=[]
	for k,v in data.items():
		if method=='binary':
			gvals=v[group]
			gnz=np.count_nonzero(gvals)
			anz=np.count_nonzero(v)
			pnull=float(anz)/len(v)
			p1=stats.binom.cdf(grouplen-gnz,grouplen,1-pnull)
			if twosided:
				p2=stats.binom.cdf(gnz,grouplen,pnull)
				p=min(p1,p2)
			else:
				p=p1
			pv={}
			pv['pval']=p
			pv['observed']=gnz
			pv['expected']=pnull*grouplen
			pv['name']=k
			allpv.append(pv)
			justp.append(p)
		elif method=='ranksum':
			rdat=v
			notgroup=np.setdiff1d(np.arange(len(v)),group)
			u,p=stats.mannwhitneyu(rdat[group],rdat[notgroup])
			pv={}
			pv['pval']=p
			pv['observed']=np.mean(rdat[group])
			pv['expected']=np.mean(rdat)
			pv['name']=k
			allpv.append(pv)
			justp.append(p)
		else:
			hs.Debug(9,'testenrichment method not supported',method)
			return False
	if fdr:
		fval=hs.fdr(justp)
		keep=np.where(np.array(fval)<=fdr)
		keep=keep[0]
	else:
		keep=np.arange(len(justp))
	plist=[]
	rat=[]
	for cidx in keep:
		plist.append(allpv[cidx])
		rat.append(np.abs(float(allpv[cidx]['observed']-allpv[cidx]['expected']))/np.mean([allpv[cidx]['observed'],allpv[cidx]['expected']]))
	si=np.argsort(rat)
	si=si[::-1]
	if printit:
		for idx,crat in enumerate(rat):
			print(plist[si[idx]])
	return(plist)


def testbactenrichment(expdat,seqs,cdb=False,bdb=False,dbexpres=False,translatestudy=False):
	"""
	EXPERIMENTAL
	test for enrichment in bacteria database categories for the bacteria in the list seqs
	enrichment is tested against manual curation (if cdb not False) and automatic curation (if bactdb not false)

	input:
	expdat : Experiment
	seqs - the sequences in the cluster
	cdb - the cooldb (manual curation) or false to skip
	bactdb - the automatic database or false to skip
	dbexpres - the assignment of values to all bacteria in the experiment (for bactdb) or false to calulate it. it is the output of bactdb.GetSeqListInfo()

	output:
	dbexpres - new if calculated
	"""

	# maybe need to keep similar freq bacteria?
	if cdb:
		hs.cooldb.testenrichment(cdb,expdat.seqs,seqs)
	if bdb:
		if not dbexpres:
			dbexpres=hs.bactdb.GetSeqListInfo(bdb,expdat.seqs,info='studies')
		seqpos=hs.findseqsinexp(expdat,seqs)
		plist=testenrichment(dbexpres,seqpos,printit=False)
		for cpv in plist:
			cname=cpv['name']
			if translatestudy:
				studyname=hs.bactdb.StudyNameFromID(bdb,cname)
			else:
				studyname=cname
			print("%s - observed %f, expected %f, pval %f" % (studyname,cpv['observed'],cpv['expected'],cpv['pval']))
	return dbexpres


def getdiffsummary(expdat,seqs,field,val1,val2=False,method='mean',threshold=0.1):
	"""
	get the fold change between 2 groups in each of the sequences in seqs
	for zech chinese ibd paper
	input:
	expdat
	seqs - the sequences to examine
	field - name of the field dividing the 2 groups
	val1 - value of the field for group 1
	val2 - value of the field for group 2 or False for all the rest (not val1)
	method:
		- 'mean' - calculate the difference in the mean of the 2 groups
		- 'binary' - calculate the difference in the presence/abscence of the 2 groups
	threshold : float
		the detection limit (out of 10k reads) - round up means below this number. should be ~10k/rarefaction

	output:
	diff - a list of the difference between the 2 groups for each sequence
	"""

	pos1=hs.findsamples(expdat,field,val1)
	if val2:
		pos2=hs.findsamples(expdat,field,val2)
	else:
		pos2=hs.findsamples(expdat,field,val1,exclude=True)

	diff=[]
	for cseq in seqs:
		if cseq in expdat.seqdict:
			seqpos=expdat.seqdict[cseq]
		else:
			diff.append[np.nan]
			continue
		if method=='mean':
			cval1=np.mean(expdat.data[seqpos,pos1])
			cval2=np.mean(expdat.data[seqpos,pos2])
		elif method=='binary':
			cval1=np.mean(expdat.data[seqpos,pos1]>0)
			cval2=np.mean(expdat.data[seqpos,pos2]>0)
		else:
			hs.Debug(9,"Unknown method %s for getdiff" % method)
			return False
		if cval1<=threshold and cval2<=threshold:
			diff.append(np.nan)
			continue
		if cval1<threshold:
			cval1=threshold
		if cval2<threshold:
			cval2=threshold
		cdiff=np.log2(cval1/cval2)
		diff.append(cdiff)
	return diff



def diffexpfastpermute(expdat,field,val1,val2=False,method='mean',numperm=[100,1000],maxfval=0.1,mineffect=0.1,permutefield=False):
	"""
	test the differential expression between 2 groups (val1 and val2 in field field)
	for bacteria that have a high difference.
	using the fast permutation test (see fastpermutepv)
	input:
	expdat
	field - the field for the 2 categories
	val1 - values for the first group
	val2 - value for the second group or false to compare to all other
	method - the test to compare the 2 groups:
		mean - absolute difference in mean frequency
		binary - abs diff in binary presence/absence
		ranksum - abs diff in rank order (to ignore outliers)
		freqpres - abs diff in frequency only in samples where bacteria is present (NOT SUPPORTED YET)
	numperm - a list of the number of random permutations to run in each step (small to large)
	maxfval - the maximal f-value (FDR) for a bacteria to keep
	mineffect : float
		the minimal effect size to keep (abs(diff)/mean)
	permutefield : string or False
		if False (default) permute all samples, otherwise permute separately on samples with same value on field permutefield

	output:
	newexp - the experiment with only significant (FDR<=maxfval) difference, sorted according to difference
	"""
	params=locals()

	exp1=hs.filtersamples(expdat,field,val1,exact=True)
	if val2:
		exp2=hs.filtersamples(expdat,field,val2,exact=True)
	else:
		exp2=hs.filtersamples(expdat,field,val1,exact=True,exclude=True)
	cexp=hs.joinexperiments(exp1,exp2)

	pos1=hs.findsamples(cexp,field,val1)
	if val2:
		pos2=hs.findsamples(cexp,field,val2)
	else:
		pos2=hs.findsamples(cexp,field,val1,exclude=True)

	minthresh=2
	dat=cexp.data
	dat[dat<minthresh]=minthresh
	dat=np.log2(dat)
	numseqs=len(cexp.seqs)

	if method=='mean':
		pass
	elif method=='ranksum':
		for idx in range(numseqs):
			dat[idx,:]=stats.rankdata(dat[idx,:])
	elif method=='binary':
		dat=(dat>np.log2(minthresh))
	elif method=='freqpres':
		hs.Debug(9,"Method not supported yet")
		return
#		dat[dat<=minthresh]=np.nan
	else:
		hs.Debug(9,"Method not supported!",method)
		return

	if permutefield:
		gvals=hs.getfieldvals(cexp,permutefield)
		ugvals=list(set(gvals))
		permutegroups=[]
		for cgval in ugvals:
			permutegroups.append(np.array(hs.findsamples(cexp,permutefield,cgval)))
		print(permutegroups)
	else:
		permutegroups=False

	pval,odif=fastpermutepv(dat,meanfunc,pos1,pos2,numperm=numperm,maxpval=maxfval,mineffect=mineffect,permutegroups=permutegroups)

	# do fdr
	fval=hs.fdr(pval)
	keep=np.where(np.array(fval)<=maxfval)
	seqlist=[]
	for cidx in keep[0]:
		seqlist.append(cexp.seqs[cidx])

	# do we need this or is reorder enough?
	newexp=hs.filterseqs(expdat,seqlist,logit=False)
	odif=odif[keep[0]]
	sv,si=hs.isort(odif)
	newexp=hs.reorderbacteria(newexp,si)
	hs.addcommand(newexp,"diffexpfastpermute",params=params,replaceparams={'expdat':expdat})
	newexp.filters.append('fast differential expression (%s) in %s between %s and %s, min effect=%d, fdr=%f' % (method,field,val1,val2,mineffect,maxfval))
	return newexp


def fastpermutepv(dat,tfunc,group1,group2,numperm=[100,1000],maxpval=0.05,mineffect=0.1,permutegroups=False):
	"""
	do a fast iterative permutation test for differential expression dat (each row is tested independently)
	by skipping the row if it is clearly above the p-value or if the effect size is not big enough
	input:
	dat : 2d np array
		row per bacteria, column per sample. NOTE it needs to be preprocessed for the appropriate method (i.e. binary,ranksum, etc)
	tfunc : the test function (takes dat and returns a vector of values (1 per row))
	maxpv : float
		the maximal p-value to test
	group1,group2 : list of integers
		the indices of group1,group2 in the data
	numperm : list of integers
		number of permutations in each step
	mineffect - the minimal effect size (abs(group difference)/mean) - keep only bacteria with at least this effect size or 0 to use all
	output:
	pvals - np array of floats
		p-value for each row (approximate if high)
	odif - np array of floats
		the statistic for the original (non-permuted) groups
	permutegroups : list of lists
		if False (default) permute all samples, otherwise permute separately on samples within each group
	"""

	numbact=np.shape(dat)[0]
	numg1=len(group1)
	numg2=len(group2)
	# calculate the true difference
	val1=tfunc(dat[:,group1])
	val2=tfunc(dat[:,group2])
	odif=val1-val2
	oeffect=np.abs(odif)/((val1+val2)/2)

	odat=copy.copy(dat)

	# init for the first permutations
	pval=np.ones(numbact)
	pval[oeffect>=mineffect]=0

	# go over permutation numbers
	for cnumperm in numperm:
		alldif=np.empty([numbact,cnumperm])
		alldif[:]=np.nan
		usebact=np.where(pval<2*maxpval)[0]
		notusebact=np.where(pval>=2*maxpval)[0]
		print("cnumperm %d, numbact %d, numnotuse %d" % (cnumperm,len(usebact),len(notusebact)))
		dat=odat[usebact,:]
		# do permutations
		for x in range(cnumperm):
			if permutegroups:
				rp=np.arange(numg1+numg2)
				for cgroup in permutegroups:
					crp=np.random.permutation(len(cgroup))
					rp[cgroup]=rp[cgroup[crp]]
				pass
			else:
				rp=np.random.permutation(numg1+numg2)
			val1=tfunc(dat[:,rp[0:numg1]])
			val2=tfunc(dat[:,rp[numg1:]])
			diff=val1-val2
			alldif[usebact,x]=diff

		# calculate the p-values
		pval=np.ones([numbact])
		for crow in range(numbact):
			cdat=alldif[crow,:]
			cdat=cdat[np.logical_not(np.isnan(cdat))]
			ccnumperm=len(cdat)
			if ccnumperm==0:
				pval[crow]=1
				continue
			cpval=float(np.sum(np.abs(cdat)>=np.abs(odif[crow])))/ccnumperm
			# need to remember we only know as much as the number of permutations - so add 1 as upper bound for fdr
			cpval=min(cpval+(1.0/ccnumperm),1)
			pval[crow]=cpval
	return pval,odif


def fastpermutecorr(expdat,field,numeric=True,method='ranksum',numperm=[100,1000],maxfval=0.05,mineffect=0.1):
	"""
	do a fast iterative permutation test for correlation to continuous field in the experiment (each row is tested independently)
	by skipping the row if it is clearly above the p-value or if the effect size is not big enough
	input:
	dat : 2d np array
		row per bacteria, column per sample. NOTE it needs to be preprocessed for the appropriate method (i.e. binary,ranksum, etc)
	tfunc : the test function (takes dat and returns a vector of values (1 per row))
	maxfv : float
		the maximal f-value to keep
	group1,group2 : list of integers
		the indices of group1,group2 in the data
	numperm : list of integers
		number of permutations in each step
	mineffect - the minimal effect size (abs correlation of the true data) - keep only bacteria with at least this effect size or 0 to use all
	output:
	pvals - np array of floats
		p-value for each row (approximate if high)
	odif - np array of floats
		the statistic for the original (non-permuted) groups
	"""
	params=locals()

	minthresh=2
	if field:
		cexp=hs.sortsamples(expdat,field=field,numeric=numeric,logit=False)

	numbact=len(cexp.seqs)
	numsamps=len(cexp.samples)

	dat=cexp.data
	dat[dat<minthresh]=minthresh
	dat=np.log2(dat)
	dat=dat-np.mean(dat,axis=1,keepdims=True)

	if method=='ranksum':
		compareto=np.arange(numsamps)
		compareto=compareto-np.mean(compareto,keepdims=True)
		# calculate the rank
		for idx in range(numbact):
			dat[idx,:]=stats.rankdata(dat[idx,:])
	else:
		hs.Debug(9,"Method not supported!",method)
		return
	odif=corrfunc(dat,compareto)
	oeffect=np.abs(odif)

	odat=copy.copy(dat)

	# init for the first permutations
	pval=np.ones(numbact)
	pval[oeffect>=mineffect]=0

	# go over permutation numbers
	for cnumperm in numperm:
		alldif=np.empty([numbact,cnumperm])
		alldif[:]=np.nan
		usebact=np.where(pval<2*maxfval)[0]
		notusebact=np.where(pval>=2*maxfval)[0]
		print("cnumperm %d, numbact %d, numnotuse %d" % (cnumperm,len(usebact),len(notusebact)))
		dat=odat[usebact,:]
		# do permutations
		for x in range(cnumperm):
			rp=np.random.permutation(numsamps)
			diff=corrfunc(dat[:,rp],compareto)
			alldif[usebact,x]=diff

		# calculate the p-values
		pval=np.ones([numbact])
		for crow in range(numbact):
			cdat=alldif[crow,:]
			cdat=cdat[np.logical_not(np.isnan(cdat))]
			ccnumperm=len(cdat)
			if ccnumperm==0:
				pval[crow]=1
				continue
			cpval=float(np.sum(np.abs(cdat)>=np.abs(odif[crow])))/ccnumperm
			# need to remember we only know as much as the number of permutations - so add 1 as upper bound for fdr
			cpval=min(cpval+(1.0/ccnumperm),1)
			pval[crow]=cpval

	# calculate the fdr
	fval=hs.fdr(pval)
	keep=np.where(np.array(fval)<=maxfval)
	seqlist=[]
	for cidx in keep[0]:
		seqlist.append(cexp.seqs[cidx])

	# do we need this or is reorder enough?
	newexp=hs.filterseqs(expdat,seqlist,logit=False)
	odif=odif[keep[0]]
	sv,si=hs.isort(odif)
	newexp=hs.reorderbacteria(newexp,si)
	hs.addcommand(newexp,"fastpermutecorr",params=params,replaceparams={'expdat':expdat})
	newexp.filters.append('fast significant correlation %s for field %s' % (method,field))
	return newexp




def meanfunc(dat,axis=1):
	"""
	calculate the mean of each row
	for the permutation tests
	"""
	return np.mean(dat,axis=axis)


def corrfunc(dat,compareto):
	"""
	calculate the partial correlation between dat and compareto
	for the permutation tests
	"""
	numseqs=np.shape(dat)[0]
	odif=np.zeros([numseqs])
	for idx in range(numseqs):
		odif[idx]=np.sum(dat[idx,:]*compareto)
	return odif


def corrfunc2(dat,compareto):
	"""
	calculate the partial correlation between dat and compareto
	for the permutation tests
	TESTING. NOT WORKING!!!!
	"""
	numseqs=np.shape(dat)[0]
	ct=np.tile(compareto,[numseqs,1])
	return dat*ct


def getsimilarfreqseqs(expdat,seqs):
	"""
	NOT READY YET!!!!
	get a list if sequences with similar frequency distribution to seqs

	input:
	expdat : Experiment
		the experiment containing the sequences
	seqs : list of sequences ('ACGT')
		the list of sequences to which we are looking for simiar sequences

	output:
	simseqs : list of sequences ('ACGT')
		a list of sequences with frequency distribution similar to seqs
	"""

	# get the original sequence frequency histogram
	ofreqs=np.zeros(len(seqs))
	for idx,cseq in enumerate(seqs):
		cseqpos=expdat.seqdict[cseq]
		ofreqs[idx]=np.mean(expdat.data[cseqpos,:])
	ffrac,edges=np.histogram(ofreqs,normed=False)
	plt.hist(ofreqs,normed=False)

	allfreqs=np.mean(expdat.data,axis=1)
	cpos=[]
	for idx,cfreq in enumerate(ffrac):
		cpos.append(np.where(np.logical_and(allfreqs>edges[idx],allfreqs<=edges[idx+1]))[0])
	done=False
	simseqs=[]
	for idx,cc in enumerate(cpos):
		print("num in seqs %d, num in all %d" % (ffrac[idx], len(cpos[idx])))
	print(edges)


def FindSimValBact(expdat,field):
	"""
	EXPERIMENTAL
	find bacteria showing similar behavior on same value in field (i.e. same family)
	input:
	expdat : Experiment
	field : string
		the field containing the values to group by

	output:
	newexp : Experiment
		sorted by effect size
	"""
	cdat=expdat.data>0
	totgvar=np.zeros([len(expdat.seqs)])
	totmean=np.mean(cdat,axis=1)
	uvals=hs.getfieldvals(expdat,field,ounique=True)
	for cval in uvals:
		cpos=hs.findsamples(expdat,field,cval)
		if not len(cpos)==2:
			continue
		cvar=(cdat[:,cpos[0]]-totmean)*(cdat[:,cpos[1]]-totmean)
		totgvar+=cvar

	totvar=np.var(cdat,axis=1)
	rat=totgvar/totvar
	sind=np.argsort(rat)
	newexp=hs.reorderbacteria(expdat,sind)
	return newexp


def fastpermutegroupsim(expdat,field,method='binary',numperm=[100,1000],maxfval=0.05,mineffect=2):
	"""
	do a fast iterative permutation test for group similarity (like twins)
	by skipping the row if it is clearly above the p-value or if the effect size is not big enough
	input:
	expdat : Experiment
	field : string
		Name of the field containing the per-group identifier
	numperm : list of integers
		number of permutations in each step
	output:
	newexp : Experiment
	"""
	params=locals()

	numbact=len(expdat.seqs)
	numsamps=len(expdat.samples)
	if method=='binary':
		dat=expdat.data>0
	else:
		hs.Debug(9,"Method not supported!",method)
		return

	uvals=hs.getfieldvals(expdat,field,ounique=True)
	groups=[]
	for cval in uvals:
		cpos=hs.findsamples(expdat,field,cval)
		if len(cpos)<2:
			continue
		if len(cpos)>2:
			cpos=cpos[:2]
		groups.append(cpos)
	hs.Debug(6,'Found %d group pairs' % len(groups))

	odif=groupfunc(dat,groups)
	oeffect=np.abs(odif)

	odat=copy.copy(dat)

	# init for the first permutations
	pval=np.ones(numbact)
	pval[oeffect>=mineffect]=0

	# go over permutation numbers
	for cnumperm in numperm:
		alldif=np.empty([numbact,cnumperm])
		alldif[:]=np.nan
		usebact=np.where(pval<2*maxfval)[0]
		notusebact=np.where(pval>=2*maxfval)[0]
		print("cnumperm %d, numbact %d, numnotuse %d" % (cnumperm,len(usebact),len(notusebact)))
		dat=odat[usebact,:]
		# do permutations
		for x in range(cnumperm):
			rp=np.random.permutation(numsamps)
			diff=groupfunc(dat[:,rp],groups)
			alldif[usebact,x]=diff

		# calculate the p-values
		pval=np.ones([numbact])
		for crow in range(numbact):
			cdat=alldif[crow,:]
			cdat=cdat[np.logical_not(np.isnan(cdat))]
			ccnumperm=len(cdat)
			if ccnumperm==0:
				pval[crow]=1
				continue
			cpval=float(np.sum(np.abs(cdat)>=np.abs(odif[crow])))/ccnumperm
			# need to remember we only know as much as the number of permutations - so add 1 as upper bound for fdr
			cpval=min(cpval+(1.0/ccnumperm),1)
			pval[crow]=cpval

	# calculate the fdr
	fval=hs.fdr(pval)
	keep=np.where(np.array(fval)<=maxfval)
	seqlist=[]
	for cidx in keep[0]:
		seqlist.append(expdat.seqs[cidx])

	# do we need this or is reorder enough?
	newexp=hs.filterseqs(expdat,seqlist,logit=False)
	odif=odif[keep[0]]
	sv,si=hs.isort(pval[keep[0]])
	newexp=hs.reorderbacteria(newexp,si)
	hs.addcommand(newexp,"fastpermutecorr",params=params,replaceparams={'expdat':expdat})
	newexp.filters.append('fast significant correlation %s for field %s' % (method,field))
	return newexp


def fastpermutetwogroupsim(expdat,twofield,twoval1,twoval2=False,field='',method='binary',numperm=[100,1000],maxfval=0.05,mineffect=2,twoway=True):
	"""
	do a fast iterative permutation test for group similarity (like twins)
	by skipping the row if it is clearly above the p-value or if the effect size is not big enough
	compare the values in the 2 groups and see where it is bigger in one group
	input:
	expdat1,expdat2 : Experiment
		the experiments to test (a division of 1 experiment)
	field : string
		Name of the field containing the per-group identifier
	numperm : list of integers
		number of permutations in each step
	twoway : bool
		True (default) to test significant differences in both directions, False to only test group1>group2
	output:
	newexp : Experiment
	"""
	params=locals()

	expdat1=hs.filtersamples(expdat,twofield,twoval1)
	if twoval2:
		expdat2=hs.filtersamples(expdat,twofield,twoval2)
	else:
		expdat2=hs.filtersamples(expdat,twofield,twoval1,exclude=True)


	numbact=len(expdat1.seqs)
	numsamps1=len(expdat1.samples)
	numsamps2=len(expdat2.samples)
	if method=='binary':
		dat1=expdat1.data>0
		dat2=expdat2.data>0
	else:
		hs.Debug(9,"Method not supported!",method)
		return

	uvals1=hs.getfieldvals(expdat1,field,ounique=True)
	groups1=[]
	for cval in uvals1:
		cpos=hs.findsamples(expdat1,field,cval)
		if len(cpos)<2:
			continue
		if len(cpos)>2:
			cpos=cpos[:2]
		groups1.append(cpos)
	hs.Debug(6,'Found %d group pairs for %s' % (len(groups1),twoval1))

	uvals2=hs.getfieldvals(expdat2,field,ounique=True)
	groups2=[]
	for cval in uvals2:
		cpos=hs.findsamples(expdat2,field,cval)
		if len(cpos)<2:
			continue
		if len(cpos)>2:
			cpos=cpos[:2]
		groups2.append(cpos)
	hs.Debug(6,'Found %d group pairs for %s' % (len(groups2),twoval2))

	odif1=(groupfunc(dat1,groups1)+0.0)/len(groups1)
	odif2=(groupfunc(dat2,groups2)+0.0)/len(groups2)
	odif=odif1-odif2
	oeffect=np.abs(odif)

	odat1=copy.copy(dat1)
	odat2=copy.copy(dat2)

	# init for the first permutations
	pval=np.ones(numbact)
	pval[oeffect>=mineffect]=0

	# go over permutation numbers
	for cnumperm in numperm:
		alldif=np.empty([numbact,cnumperm])
		alldif[:]=np.nan
		usebact=np.where(pval<2*maxfval)[0]
		notusebact=np.where(pval>=2*maxfval)[0]
		print("cnumperm %d, numbact %d, numnotuse %d" % (cnumperm,len(usebact),len(notusebact)))
		dat1=odat1[usebact,:]
		dat2=odat2[usebact,:]
		# do permutations
		for x in range(cnumperm):
			rp1=np.random.permutation(numsamps1)
			rp2=np.random.permutation(numsamps2)
			diff1=(groupfunc(dat1[:,rp1],groups1)+0.0)/len(groups1)
			diff2=(groupfunc(dat2[:,rp2],groups2)+0.0)/len(groups2)
			diff=diff1-diff2
			alldif[usebact,x]=diff

		# calculate the p-values
		pval=np.ones([numbact])
		for crow in range(numbact):
			cdat=alldif[crow,:]
			cdat=cdat[np.logical_not(np.isnan(cdat))]
			ccnumperm=len(cdat)
			if ccnumperm==0:
				pval[crow]=1
				continue
			if twoway:
				cpval=float(np.sum(np.abs(cdat)>=np.abs(odif[crow])))/ccnumperm
			else:
				cpval=float(np.sum(cdat>=odif[crow]))/ccnumperm
			# need to remember we only know as much as the number of permutations - so add 1 as upper bound for fdr
			cpval=min(cpval+(1.0/ccnumperm),1)
			pval[crow]=cpval

	# calculate the fdr
	fval=hs.fdr(pval)
	keep=np.where(np.array(fval)<=maxfval)
	seqlist=[]
	for cidx in keep[0]:
		seqlist.append(expdat.seqs[cidx])

	# do we need this or is reorder enough?
	newexp=hs.filterseqs(expdat,seqlist,logit=False)
	odif=odif[keep[0]]
	sv,si=hs.isort(pval[keep[0]])
	newexp=hs.reorderbacteria(newexp,si)
	hs.addcommand(newexp,"fastpermutetwogroupsim",params=params,replaceparams={'expdat':expdat})
	newexp.filters.append('fast permutation for 2 groups similarity %s for field %s' % (method,field))
	return newexp


def groupfunc(dat,groups):
	"""
	for group similarity
	count the number of 1s the appear in the 2 samples of the same group
	"""
	odat=np.zeros([np.shape(dat)[0]])
	for cgroup in groups:
		odat += dat[:,cgroup[0]]*dat[:,cgroup[1]]
	return odat


def getsubtaxdict(expdat,separator=';',removemin=2,removemax=250):
	"""
	get a dict of all sub taxonomies and sequence positions associated with it

	input:
	expdat : Experiment
	separator : string
		The sperator between different taxonomic levels
	removemin : int
		The minimal number of bacteria in a group in order to keep the group (0 to keep all)
	removemax : int
		The maximal number of bacteria in a group in order to keep the group (0 to keep all)


	output:
	taxgroups : dictionary of Name (string) ,pos (list of int)
		The taxonomic groups
	"""
	taxgroups={}
	for idx,ctax in enumerate(expdat.tax):
		if ctax in taxgroups:
			taxgroups[ctax].append(idx)
		else:
			taxgroups[ctax]=[idx]
		for cpos,cchar in enumerate(ctax):
			if cchar==separator:
				subtax=ctax[:cpos]
				if subtax in taxgroups:
					taxgroups[subtax].append(idx)
				else:
					taxgroups[subtax]=[idx]
	# and remove groups with 1 bacteria
	for k,v in taxgroups.items():
		taxgroups[k]=np.array(v)
		if len(v)<removemin:
			del taxgroups[k]
		elif removemax>0:
			if len(v)>removemax:
				del taxgroups[k]
	return taxgroups


def fastpermuteoverlap(expdat,tfunc,taxgroups={},numperm=[100,1000],maxpval=0.05,mineffect=0.1,method='binary'):
	"""
	do a fast iterative permutation test for taxonomic group overlap
	by skipping the row if it is clearly above the p-value or if the effect size is not big enough
	input:
	expdat : Experiment
	taxgroups : dictionary of Name (string) ,pos (list of int) or empty list
		The taxonomic groups to test for overlap - name and position in the data matrix. Empty (default) to calculate here
	tfunc : the test function (takes dat and a list of taxonomic group positions, returns value per group
	maxpv : float
		the maximal p-value to test
	numperm : list of integers
		number of permutations in each step
	mineffect - the minimal effect size (abs(group difference)/mean) - keep only bacteria with at least this effect size or 0 to use all
	output:
	taxgroups - dict of name,list of positions
	pvals - np array of floats
		p-value for each taxonomic group (approximate if high)
	odif - np array of floats
		the statistic for the original (non-permuted) groups
	"""

	dat=expdat.data
	# if no taxonomic group dict provided, make it ourselves
	if len(taxgroups)==0:
		taxgroups=getsubtaxdict(expdat)
	taxglist=[]
	taxgnames=[]
	for k,v in taxgroups.items():
		taxglist.append(v)
		taxgnames.append(k)

	odat=copy.copy(expdat.data)
	if method=='binary':
		odat=odat>0
	numgroups=len(taxgroups)
	numsamps=len(expdat.samples)
	numseqs=len(expdat.seqs)
	# calculate the true difference
	odif=tfunc(odat,taxglist)

	# init for the first permutations
	pval=np.ones(numgroups)
#	pval[oeffect>=mineffect]=0

	# go over permutation numbers
	for cnumperm in numperm:
		alldif=np.empty([numgroups,cnumperm])
		alldif[:]=np.nan
		usegroup=np.where(pval<2*maxpval)[0]
		notusegroup=np.where(pval>=2*maxpval)[0]
		print("cnumperm %d, numbact %d, numnotuse %d" % (cnumperm,len(usegroup),len(notusegroup)))
		dat=odat
		# do permutations
		for x in range(cnumperm):
			for idx in range(numseqs):
				rp=np.random.permutation(numsamps)
				dat[idx,:]=dat[idx,rp]
			diff=tfunc(dat,taxglist)
			alldif[:,x]=diff

		# calculate the p-values
		pval=np.ones([numgroups])
		for crow in range(numgroups):
			cdat=alldif[crow,:]
			cdat=cdat[np.logical_not(np.isnan(cdat))]
			ccnumperm=len(cdat)
			if ccnumperm==0:
				pval[crow]=1
				continue
			cpval=float(np.sum(cdat>=odif[crow]))/ccnumperm
			# need to remember we only know as much as the number of permutations - so add 1 as upper bound for fdr
			cpval=min(cpval+(1.0/ccnumperm),1)
			pval[crow]=cpval
	return pval,odif,taxgnames,taxglist


def getgroupcover(dat,taxgroups):
	"""
	for fastpermuteoverlap. Get the number of samples that have at least 1 bacteria out of the tax group
	input:
	dat - the experiment data matrix
	taxgroups : list of list of int
		a list of bacteria positions per group

	output:
	cover : numpy array of floats
		Number of samples that have coverage for each taxgroup
	"""
	cover=np.zeros([len(taxgroups)])
	for idx,cgroup in enumerate(taxgroups):
		persamp=np.sum(dat[cgroup,:],axis=0)
		cover[idx]=np.sum(persamp>0)
	return cover


def getgroupcoverstd(dat,taxgroups):
	"""
	for fastpermuteoverlap. Get the std of the sum of all bacteria in the group
	input:
	dat - the experiment data matrix
	taxgroups : list of list of int
		a list of bacteria positions per group

	output:
	cover : numpy array of floats
		1/std of each group
	"""

	cover=np.zeros([len(taxgroups)])
	for idx,cgroup in enumerate(taxgroups):
		persamp=np.sum(dat[cgroup,:],axis=0)
		cover[idx]=np.std(persamp)
	return cover




def get2expcorr(exp1,exp2,field,method='spearman',maxfval=0.1):
	"""
	test the correlation between the same bacteria in 2 experiments
	input:
	exp1,exp2 : Experiment
		the experiments containing the same bacteria
	field - the field containing the same unique value for both experiments (i.e. subject_id etc.)
	method - the test to compare the 2 groups:
		'spearman' - nonparametric correlation
		'binary' - spearman correlation of presence/absence
	maxfval - the maximal f-value (FDR) for a bacteria to keep

	output:
	newexp : Experiment
		the experiment with 2 lines per bacteria - from exp1 and exp2.
		with only significant (FDR<=maxfval) difference, sorted according to difference
	"""
	params=locals()

	# make sure 1 sample per value - otherwise join
	exp1=hs.filtersimilarsamples(exp1,field,method='mean')
	exp2=hs.filtersimilarsamples(exp2,field,method='mean')

	# reorder both experiments the same
	e2vals=hs.getfieldvals(exp2,field)
	exp1=hs.filtersamples(exp1,field,e2vals,exact=True)
	e1vals=hs.getfieldvals(exp1,field)
	neworder=[]
	for idx,cval in enumerate(e1vals):
		neworder.append(e2vals.index(cval))
	exp2=hs.reordersamples(exp2,neworder)

	minthresh=2
	len1=len(exp1.samples)
	len2=len(exp2.samples)
	dat1=copy.copy(exp1.data)
	dat2=copy.copy(exp2.data)
	dat1[dat1<minthresh]=minthresh
	dat2[dat2<minthresh]=minthresh
	dat1=np.log2(dat1)
	dat1=np.log2(dat1)
	numseqs=len(exp1.seqs)

	eps=0.000001

	if method=='binary':
		dat1=(dat1>np.log2(minthresh))
		dat2=(dat2>np.log2(minthresh))
	elif method=='spearman':
		pass
	else:
		hs.Debug(9,"Method not supported!",method)
		return

	pval=np.ones([numseqs])
	odif=np.zeros([numseqs])
	for pos1,cseq in enumerate(exp1.seqs):
		if cseq in exp2.seqdict:
			odif[pos1],pval[pos1]=stats.spearmanr(dat1[pos1,:],dat2[exp2.seqdict[cseq],:])

	# NOTE: maybe use better fdr (this one not tested...)
	fval=hs.fdr(pval)
	keep=np.where(np.array(fval)<=maxfval)
	seqlist=[]
	for cidx in keep[0]:
		seqlist.append(exp1.seqs[cidx])

	# do we need this or is reorder enough?
	newexp1=hs.filterseqs(exp1,seqlist,logit=False)
	newexp2=hs.filterseqs(exp2,seqlist,logit=False)
	odif=odif[keep[0]]
	sv,si=hs.isort(odif)
	newexp1=hs.reorderbacteria(newexp1,si)
	newexp2=hs.reorderbacteria(newexp2,si)
	newexp1.odif=sv
	newexp2.odif=sv
	newexp=hs.copyexp(newexp1)
	newseqs=[]
	newtax=[]
	newsids=[]
	newodif=[]
	newexp.data=np.vstack([newexp.data,newexp.data])
	for idx,cseq in enumerate(newexp1.seqs):
		newodif.append(newexp1.odif[idx])
		newodif.append(newexp1.odif[idx])
		newsids.append(newexp1.sids[idx])
		newsids.append(newexp1.sids[idx])
		newseqs.append(newexp1.seqs[idx])
		cseq2=newexp1.seqs[idx]+'2'
		newseqs.append(cseq2)
		newexp.seqdict[newexp1.seqs[idx]]=idx*2
		newexp.seqdict[cseq2]=idx*2+1
		newtax.append(newexp1.tax[idx]+'1')
		newtax.append(newexp1.tax[idx]+'2')
		newexp.data[idx*2,:]=newexp1.data[idx,:]
		newexp.data[idx*2+1,:]=newexp2.data[idx,:]
	newexp.seqs=newseqs
	newexp.tax=newtax
	newexp.sids=newsids
	newexp.odif=newodif
	hs.addcommand(newexp,"get2expcorr",params=params,replaceparams={'exp1':exp1,'exp2':exp2})
	newexp.filters.append('two experiment correlation (%s) in %s' % (method,field))
	return newexp


def groupbacteria(expdat,minreads=0,uselog=True):
	"""
	group the bacteria into distinct clusters

	input:
	expdat : Experiment
	minreads : int
		minimal reads for a bacteria to keep.
	uselog : bool
		True to log scale the data, False to not

	output:
	outexp : experiment
		with the sum of all bacteria in each group
	seqpergroup : list of list of seqs (ACGT)
		the sequences in each group
	taxpergroup : list of list of str
		the taxonomies in each group
	"""
	newexp=hs.filterminreads(expdat,minreads,logit=False)
	# normalize each row (bacteria) to sum 1
	dat=copy.copy(newexp.data)
	if uselog:
		dat[dat<=2]=2
		dat=np.log2(dat)
	dat=scale(dat,axis=1,copy=False)
	# cluster
	dm=spatial.distance.pdist(dat,metric='euclidean')
	dm=spatial.distance.squareform(dm)
	af = AffinityPropagation().fit(dm)
	cluster_centers_indices = af.cluster_centers_indices_
	labels = af.labels_

	n_clusters_ = len(cluster_centers_indices)

	hs.Debug(6,"number of clusters - %d" % n_clusters_)
	hs.Debug(1,labels)

	seqpergroup=[]
	taxpergroup=[]
	outexp=hs.zerobacteria(newexp)
	uids=list(set(labels))
	for cid in uids:
		cpos=np.where(labels==cid)[0]
		cseqs=[]
		ctax=[]
		for ccpos in cpos:
			cseqs.append(newexp.seqs[ccpos])
			ctax.append(newexp.tax[ccpos])
		seqpergroup.append(cseqs)
		taxpergroup.append(ctax)
		if len(cpos)>1:
			print('-------')
			print(cid)
			print(newexp.data[cpos,:4])
			newfreq=np.sum(newexp.data[cpos,:],axis=0)
		else:
			newfreq=newexp.data[cpos,:]
		hs.insertbacteria(outexp,freqs=newfreq,seq=str(cid),tax=str(cid),logit=False)

	for cidx,cnum in enumerate(labels):
		newexp.tax[cidx]=str(cnum)
	newexp=hs.sortbacteria(newexp)

	return outexp,seqpergroup,taxpergroup
