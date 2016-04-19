#!/usr/bin/env python


"""
heatsequer other plot functions module
"""

# amnonscript

__version__ = "0.9"

import heatsequer as hs

import numpy as np
import matplotlib as mpl
#mpl.use('Qt4Agg')
import matplotlib.pyplot as plt
#from matplotlib.pyplot import *
import copy

from pdb import set_trace as XXX


def plotseqfreq(expdat,seqs,toaxis=False,xfield=False,normalizey=False):
	"""
	plot the frequency of sequences in seq as a function of the sortfield
	input:
	expdat : Experiment
	seqs : list of sequence strings
		a list of sequnces (acgt) to plot
	toaxis : matplotlib axis
		if not empty - the axis to plot to, False plot a new figure
	xfield : string
		if not empty - space the points on the x axis according to (numeric) value of in xfield, False - just according to the sorted order
	normalizey : bool
		True: normalize all y values to 0-1, False: no normalization
	"""


	if xfield:
		xdat=hs.tofloat(hs.getfieldvals(expdat,xfield))
	else:
		xdat=range(len(expdat.samples))
	sv,si=hs.isort(xdat)
	if not toaxis:
		plt.figure()
		toaxis=plt.gca()
	labels=[]
	for cseq in seqs:
		if cseq not in expdat.seqdict:
			continue
		spos=expdat.seqdict[cseq]
		cdat=expdat.data[spos,si]
		if normalizey:
			cdat=cdat/sum(cdat)
		toaxis.plot(sv,cdat)
		labels.append(str(expdat.sids[spos])+'-'+expdat.tax[spos])
	labels=hs.clipstrings(labels,20,reverse=True)
	toaxis.legend(labels,prop={'size':6})
	toaxis.set_xticks(xdat)
#	toaxis.set_xticklabels(labs,rotation=45,ha='right')


def analyzenumreads(expdat,blanks=['blank','empty']):
	"""
	Plot a graph of cumulative number of reads per sample in the biom table.
	we assume all samples with an element of blanks in their name a blanks (for a different graph)
	input:
	expdat
	blanks - a list of words which if appearing in the sample name indicate they are blanks (case insensitive!)
	"""

	nreads=expdat.origreads
	nreadsb=[]
	nreadsnb=[]
	for idx,csamp in enumerate(expdat.samples):
		isblank=False
		for cblank in blanks:
			if cblank.lower() in csamp.lower():
				isblank=True
		if isblank:
			nreadsb.append(expdat.origreads[idx])
		else:
			nreadsnb.append(expdat.origreads[idx])

	tsamps=len(nreads)
	tsampsb=max(len(nreadsb),1)
	tsampsnb=max(len(nreadsnb),1)
	nreads=np.array(nreads)
	nreadsb=np.array(nreadsb)
	nreadsnb=np.array(nreadsnb)
	y=[]
	yb=[]
	ynb=[]
	x=[]
	print(max(nreads))
	for cx in xrange(0,int(max(nreads)),500):
		y.append(float(sum(nreads>=cx))/tsamps)
		yb.append(float(sum(nreadsb>=cx))/tsampsb)
		ynb.append(float(sum(nreadsnb>=cx))/tsampsnb)
		x.append(cx)
	plt.figure()
	plt.plot(x,y)
	plt.plot(x,yb)
	plt.plot(x,ynb)
	plt.title(expdat.studyname+' (%d samples)' % len(expdat.samples))
	plt.xlabel('number of reads')
	plt.ylabel('fraction of samples with >= reads')
	plt.legend(['all','blanks','non blanks'])
	plt.show()


def plotnucdistribution(expdat,position):
	"""
	EXPERIMENTAL
	get the distribution of nucleotides in the positions in position
	for conserved/nonconserved positions
	note this is unweighted
	input:
	expdat
	position - a list of positions (0 based) to test

	output:
	"""

	retv=np.zeros((len(position),6))
	for cseq in expdat.seqs:
		cseqn=hs.SeqToArray(cseq)
		for idx,cpos in enumerate(position):
			retv[idx,cseqn[cpos]]+=1
	plt.figure()
	for cnuc in range(np.size(retv,axis=1)):
		for crow in range(np.size(retv,axis=0)-1):
			plt.bar(np.arange(np.size(retv,axis=0)),retv[crow+1,:],bottom=retv[crow,:])
	return (retv)



def plotdiffsummary2(expdatlist1,expdatlist2,seqs1,seqs2,field,val1,val2=False,method='mean',sortit=True,threshold=0.1):
	"""
	plot a heat map for 2 results (i.e. 16s and KO predictions) for zech chinese ibd study
	input:
	expdatlist1,expdatlist2 - a list of experiments to plot (row per experiment - all must contain field and val1,val2 in it)
	seqs1,seqs2 - the sequences to examine
	field - name of the field dividing the 2 groups
	val1 - value of the field for group 1 (or a list of values 1 per experiment)
	val2 - value of the field for group 2 or False for all the rest (not val1) (or a list of values 1 per experiment)
	method:
		- mean - calculate the difference in the mean of the 2 groups
	sortit - True to sort according to difference in the first expdat, False to use the order in seqs
	"""
	diff1,names1=plotdiffsummary(expdatlist1,seqs1,field,val1,val2,method,sortit,threshold=threshold)
	diff2,names2=plotdiffsummary(expdatlist2,seqs2,field,val1,val2,method,sortit,threshold=threshold)
	print(np.shape(diff1))
	print(np.shape(diff2))
	diff=np.vstack([diff1,diff2])
	maxdiff=np.nanmax(np.abs(diff))
	plt.figure()
	plt.imshow(diff,interpolation='nearest',aspect='auto',cmap=plt.get_cmap("coolwarm"),clim=[-maxdiff,maxdiff])
	plt.colorbar()
	plt.title("log2 fold change between %s and %s in field %s" % (val1,val2,field))
	plt.plot([-0.5,len(expdatlist1)-0.5],[np.shape(diff1)[0]-0.5,np.shape(diff1)[0]-0.5],'k')
	plt.xticks(np.arange(len(names1)),names1,rotation=90)
	plt.autoscale(tight=True)


def plotdiffsummary(expdatlist,seqs,field,val1,val2=False,method='mean',sortit=True,threshold=0.1,ptitle=False,showallfirstexp=True):
	"""
	plot a heat map for the fold change in each experiment in expdatlist
	for the log2 foldchange between the 2 groups (val1,val2 values in field)
	for zech chinese ibd paper
	input:
	expdatlist - a list of experiments to plot (row per experiment - all must contain field and val1,val2 in it)
	seqs - the sequences to examine
	field - name of the field dividing the 2 groups
	val1 - value of the field for group 1 (or a list of values 1 per experiment)
	val2 - value of the field for group 2 or False for all the rest (not val1) (or a list of values 1 per experiment)
	method:
		- mean - calculate the difference in the mean of the 2 groups
	sortit - True to sort according to difference in the first expdat, False to use the order in seqs
	threshold - minimum value of stat for ratio calculation (otherwise rounded up to threshold)
	ptitle - name of figure of False for auto title
	showallfirstexp : bool
		True - show all sequences, False - show only sequences present in at least one other study except the first

	output:
	diffsum - the same as the plotted heatmap (row per otu, column per experiment)
	expnames - names (studyname) of the experiments plotted (for label)
	otus - the otu sequences for the rows
	"""

	if not(type(val1) is list):
		tval1=val1
		val1=[]
		for cexp in expdatlist:
			val1.append(tval1)
	if not(type(val2) is list):
		tval2=val2
		val2=[]
		for cexp in expdatlist:
			val2.append(tval2)
	diff=np.array(hs.getdiffsummary(expdatlist[0],seqs,field,val1[0],val2[0],method,threshold=threshold))
	odiff=copy.copy(diff)
	odiffnotnan=np.where(np.isfinite(odiff))[0]
	diffsum=[]
	for cidx,cexp in enumerate(expdatlist[1:]):
		cdiff=np.array(hs.getdiffsummary(cexp,seqs,field,val1[cidx+1],val2[cidx+1],method,threshold=threshold))
		diff=np.vstack([diff,cdiff])
		notnan=np.where(np.isfinite(cdiff))[0]
		notnan=np.intersect1d(notnan,odiffnotnan)
		if len(notnan)>0:
			cdiffsum=float(np.sum((cdiff[notnan]>0)==(odiff[notnan]>0)))/len(notnan)
		else:
			cdiffsum=np.nan
		diffsum.append(cdiffsum)

	# remove all NaN lines (not enough reads for threshold)
	if showallfirstexp:
		nanlines=np.where(~np.isnan(diff).all(axis=0))[0]
	else:
		nanlines=np.where(~np.isnan(diff[1:,:]).all(axis=0))[0]
	diff=diff[:,nanlines]
	otus=hs.reorder(seqs,nanlines)

	if sortit:
		si=np.argsort(diff[0,:])
		diff=diff[:,si]
		otus=hs.reorder(otus,si)
	plt.figure()
	maxdiff=np.nanmax(np.abs(diff))
	diff=np.transpose(diff)
	plt.imshow(diff,interpolation='nearest',aspect='auto',cmap=plt.get_cmap("coolwarm"),clim=[-maxdiff,maxdiff])
	plt.colorbar()
	if ptitle:
		plt.title(ptitle)
	else:
		plt.title("log2 fold change between %s and %s in field %s" % (val1,val2,field))
	expnames=[]
	for cexp in expdatlist:
		expnames.append(cexp.studyname)
	plt.xticks(np.arange(len(expnames)),expnames,rotation=45)
	plt.tight_layout()
	plt.show()
	return diff,expnames,otus



def showleakage(expdat,seq,wwpp=['1','2','3','4','5','6','7','8']):
	plt.figure()
	newexp=hs.filterseqs(expdat,[seq])
	for idx,cplate in enumerate(wwpp):
		print(cplate)
		cexp=hs.filtersamples(newexp,'primerplate_int',cplate,exact=True)
		print(len(cexp.samples))
		plt.subplot(3,3,idx+1)
		plt.title(cplate)
		cmat=np.empty([8,12])
		cmat[:] = np.NAN
		for idx2,csamp in enumerate(cexp.samples):
			crow=ord(cexp.smap[csamp]['Row'].lower())-ord('A'.lower())
			ccol=int(cexp.smap[csamp]['column_int'])-1
			cval=cexp.data[0,idx2]
			if cval>0:
				cval=np.log2(cval)
			else:
				cval=-5
			cmat[crow,ccol]=cval
		plt.imshow(cmat,interpolation='nearest',aspect='auto',clim=[-5,10],cmap=plt.get_cmap("coolwarm"))
		for idx2,csamp in enumerate(cexp.samples):
			crow=ord(cexp.smap[csamp]['Row'].lower())-ord('A'.lower())
			ccol=int(cexp.smap[csamp]['column_int'])-1
			isntc=int(cexp.smap[csamp]['NTC_bool'])
			if isntc:
				plt.text(ccol,crow,'x')


def plottimeseries(expdat,idfield,timefield,seq,toaxis=False,numeric=True):
	"""
	plot a line plot for a timeseries with a different line for each individual
	input:
	expdat : Experiment
	idfield : string
		name of the per individual id field
	timefield : string
		name of the timepoint field
	seq : string (ACGT)
		the sequence to plot
	toaxis : axis
		False (default) to draw in a new figure
	numeric : bool
		True (default) if timefield is numeric
	"""

	if not toaxis:
		plt.figure()
		toaxis=plt.gca()

	newexp=hs.sortsamples(expdat,timefield,numeric=numeric)
	times=hs.getfieldvals(expdat,timefield)
	utimes=list(set(times))
	if numeric:
		utimes=hs.tofloat(utimes)
	utimes=np.sort(utimes)
	timesd=hs.listtodict(times)
	ids=hs.getfieldvals(newexp,idfield)
	for cid in list(set(ids)):
		x=[]
		y=[]
		for ct in hs.findsamples(newexp,idfield,cid):
			y.append(newexp.data[newexp.seqdict[seq],ct])
			if numeric:
				x.append(float(times[ct]))
			else:
				x.append(timesd[times[ct]])
		toaxis.plot(x,y)
	toaxis.legend(list(set(ids)))
	toaxis.title('lala')



def plotgroupbar(expdat,field,seqs=[],type='meanse',uvals=[]):
	"""
	plot a bar graph of the mean frequency of each sequence in seqs as a function of field.
	for skin-gallo analysis

	input:
	expdat : Experiment
	field : str
		name of the field to group by
	seqs : list of str (ACGT)
		list of sequences to plot or [] for all sequences in expdat (default)
	type : str
		'meanstd' (default) - plot mean + standard error per group
	uvals : list of str
		the order of the unique field values to use, or [] for all (default)
	"""
	colors=['b','r','k','g','m','c','y']
	if len(uvals)==0:
		vals=hs.getfieldvals(expdat,field)
		uvals=list(set(vals))
	if len(seqs)==0:
		seqs=expdat.seqs
	else:
		expdat=hs.filterseqs(expdat,seqs)

	numseqs=len(seqs)
	numgroups=len(uvals)
	gmeans=np.zeros([numseqs,numgroups])
	gstd=np.zeros([numseqs,numgroups])
	ldat=[]

	width = 0.7/numgroups
	fig, ax = plt.subplots()
	for idx,cval in enumerate(uvals):
		cpos=hs.findsamples(expdat,field,cval)
		numnotnan=len(cpos)-np.sum(np.isnan(expdat.data[:,cpos]),1)
		gmeans[:,idx]=np.nanmean(expdat.data[:,cpos],axis=1)
		gstd[:,idx]=np.nanstd(expdat.data[:,cpos],axis=1)/np.sqrt(numnotnan)
		ind = np.arange(numseqs)+0.7-idx*width
#		rects = ax.bar(ind, gmeans[:,idx]/100, width, yerr=gstd[:,idx]/100,color=colors[np.mod(idx,len(colors))])
		rects = ax.barh(ind, gmeans[:,idx], width, xerr=gstd[:,idx],color=colors[np.mod(idx,len(colors))])
		ldat.append(rects[0])

	# # add some text for labels, title and axes ticks
	ax.set_xlabel('number of reads')
	# ax.set_title('Scores by group and gender')
	ax.set_yticks(np.arange(numseqs)+0.45)
	ylabs=[]
	for cseq in seqs:
		clab=hs.getnicetax(expdat.tax[expdat.seqdict[cseq]])
		ylabs.append(clab)
	ax.set_yticklabels(ylabs)
	ax.legend(ldat, uvals,loc='best')
	plt.tight_layout()
#	cx.ofig.canvas.draw()

#	ax.set_xscale('log')
