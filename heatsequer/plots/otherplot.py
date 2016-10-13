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
	for cx in np.arange(0,int(max(nreads)),500):
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
	hs.Debug(6,'total samples %d, blanks %d, other %d' % (tsamps,tsampsb,tsampsnb))
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


def plottimeseries2(expdat,idfield,timefield,seq,toaxis=False,numeric=True,showscatter=True,showlegend=True,colordict=None,showtitle=True,scaleto=None,ylabel='perceent of reads'):
	"""
	plot a line plot for a timeseries with a different line for each individual
	with mean of points in same timepoint and scatter
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
	showscatter : bool
		True (default) to plot the scatter, False to plot just the mean
	showlegend : bool
		True (default) to show the legend, False to not show
	colordict : None or dict of {label:color}
		if None, use default colors, if not - use color specified for each label
	showtitle : bool
		True (default) to show figure title, False to hide
	scaleto : int or None
		None (default) to keep the values, otherwise int to multiply each number by int
	ylabel : str or None
		y axis label or None to not show y label
	"""

	coloridx=0
	colors=['o-b','o-r','o-k','o--g','o--c','o--m','o:y']
	if not toaxis:
		plt.figure()
		toaxis=plt.gca()

	oldexp=hs.sortsamples(expdat,timefield,numeric=numeric)
	times=hs.getfieldvals(expdat,timefield)
	utimes=list(set(times))
	if numeric:
		sv,si=hs.isort(hs.tofloat(utimes))
		utimes=hs.reorder(utimes,si)
	else:
		utimes=np.sort(utimes)
	ids=hs.getfieldvals(oldexp,idfield)
	allids=list(set(ids))
	allids.sort()
	for cid in allids:
		newexp=hs.filtersamples(oldexp,idfield,cid)
		if colordict is None:
			ccolor=colors[coloridx]
			coloridx+=1
			if coloridx>=len(colors):
				coloridx=0
		else:
			ccolor=colordict[cid]
		x=[]
		y=[]
		xp=[]
		yp=[]
		for ctime in utimes:
			tpoints=hs.findsamples(newexp,timefield,ctime)
			if not tpoints:
				continue
			cx=[]
			cy=[]
			for cpoint in tpoints:
				cx.append(float(ctime))
				currenty=newexp.data[newexp.seqdict[seq],cpoint]/100
				if scaleto is not None:
					currenty=currenty*scaleto
				cy.append(currenty)
			xp+=cx
			yp+=cy
			x.append(np.mean(cx))
			y.append(np.mean(cy))
		toaxis.plot(x,y,ccolor)
		if showscatter:
			toaxis.plot(xp,yp,'o',color=ccolor)
	if showlegend:
		toaxis.legend(allids,loc='best')
	if showtitle:
		plt.title(expdat.tax[expdat.seqdict[seq]])
	else:
		print(expdat.tax[expdat.seqdict[seq]])
	plt.xlabel(timefield)
	if ylabel is not None:
		plt.ylabel(ylabel)


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


def plottwoexperimentreads(exp1,exp2):
	"""
	plot correlation between reads of 2 experiments
	input:
	exp1,exp2 : ExpClass
		the experiments to compare (should have approx. same samples)
	"""
	allsamples=exp1.samples
	allsamples.extend(exp2.samples)
	allsamples=list(set(allsamples))
	plt.figure()
	plt.xlabel(exp1.tablefilename)
	plt.ylabel(exp2.tablefilename)
	mv=0
	for csamp in allsamples:
		reads1=0
		reads2=0
		if csamp in exp1.samples:
			reads1=float(exp1.smap[csamp]['origReads'])
		if csamp in exp2.samples:
			reads2=float(exp2.smap[csamp]['origReads'])
		print('sample %s read1 %d read2 %d' % (csamp,reads1,reads2))
		plt.plot(reads1,reads2,'x')
		mv=max(mv,reads1)
		mv=max(mv,reads2)
	plt.plot([0,mv],[0,mv],'k')





def barplotexp(exp,bseq='',sortby=False,numeric=False,showxall=False,uselog=True,showxlabel=True,subline='',lowcutoff=1,ptitle='',showline=True,linewidth=1):
	"""
	Plot an experiment
	input:
	exp - from load()
	sortby - name of mapping file field to sort by or Flase to not sort
	numeric - True if the field is numeric
	minreads - minimum number of reads per bacteria in order to show it or 0 to show all
	rangeall - True to show all frequencies in image scale, false to saturate at 10%
	seqdb - the SRBactDB database (from bactdb.load)
	cdb - the cool sequences database (from cooldb.load), or None (default) to use the heatsequer loaded cdb
	showline - if True plot lines between category values
	ontofig - name of ontology to plot for bactdb or false to no plot
	usegui - True use a gui for otu summary, False just print
	showxall - True to show all sample names when not sorting, False to show no more than 10
	showcolorbar - True to plot the colorbar. False to not plot
	ptitle - name of the figure or False to show processing history as name
	lowcutoff - minimal value for read (for 0 log transform) - the minimal resolution - could be 10000*2/origreads
	showxlabel : bool
		True to show the x label (default), False to hide it
	colormap : string or False
		name of colormap or False (default) to use mpl default colormap
	colorrange : [min,max] or False
		[min,max] to set the colormap range, False to use data min,max (default) as specified in rangeall
	subline : str
		Name of category for subline plotting or '' (Default) for no sublines

	output:
	newexp - the plotted experiment (sorted and filtered)
	ax - the plot axis
	"""

	hs.Debug(1,"BarPlot experiment %s" % exp.studyname)
	hs.Debug(1,"Commands:")
	for ccommand in exp.commands:
		hs.Debug(1,"%s" % ccommand)
	exp=hs.filterseqs(exp,bseq)
	vals=[]
	if sortby:
		hs.Debug(1,"Sorting by field %s" % sortby)
		for csamp in exp.samples:
			vals.append(exp.smap[csamp][sortby])
		if numeric:
			hs.Debug(1,"(numeric sort)")
			vals=hs.tofloat(vals)
		svals,sidx=hs.isort(vals)
		newexp=hs.reordersamples(exp,sidx)
	else:
		hs.Debug(1,"No sample sorting")
		svals=hs.getfieldvals(exp,'#SampleID')
		newexp=hs.copyexp(exp)

	ldat=newexp.data
	if uselog:
		hs.Debug(1,"Using log, cutoff at %f" % lowcutoff)
		ldat[np.where(ldat<lowcutoff)]=lowcutoff
		ldat=np.log2(ldat)

	plt.figure()
	hs.Debug(1,"plotting bar plot")
	bdat=np.sum(ldat,axis=0)/100
	if subline:
		lsvals=hs.getfieldvals(newexp,subline)
		print(svals)
		ulsvals=list(set(lsvals))
		lsvals=np.array(lsvals)
		print(svals)
		for idx,cval in enumerate(ulsvals):
			ppos=np.where(lsvals==cval)[0]
			ccolor=getcolor(idx)
			plt.bar(ppos,bdat[ppos],color=ccolor)
		plt.legend(ulsvals,loc='best')
	else:
		plt.bar(np.arange(len(bdat)),bdat)
	ax=plt.gca()
	plt.ylabel('fraction of reads belonging to sequence')

	if not ptitle:
		hs.Debug(1,"Showing filters in title")
		if (len(newexp.filters))>4:
			cfilters=[newexp.filters[0],'...',newexp.filters[-2],newexp.filters[-1]]
		else:
			cfilters=newexp.filters
		cfilters=hs.clipstrings(cfilters,30)
		ptitle='\n'.join(cfilters)
	plt.title(ptitle,fontsize=10)

	if showline:
		hs.Debug(1,"Showing lines")
		labs=[]
		labpos=[]
		linepos=[]
		minpos=0
		svals.append('end')
		for idx,cval in enumerate(svals[:-1]):
			if cval==svals[idx+1]:
				continue
			labpos.append(minpos-0.5+float(idx+1-minpos)/2)
			minpos=idx+1
			linepos.append(idx)
			labs.append(cval)
		hs.Debug(1,"number of lines is %d" % len(linepos))
		if showxlabel:
			ax.set_xticks(labpos)
			ax.set_xticklabels(labs,rotation=45,ha='right')
		for cx in linepos:
			plt.plot([cx,cx],[-0.5,np.max(bdat)],'k',linewidth=linewidth)
	else:
		hs.Debug(1,"Not showing lines")
		if showxall or len(newexp.samples)<=10:
			hs.Debug(1,"less than 10 samples, showing all sample names")
			ax.set_xticklabels(svals,rotation=90)
			ax.set_xticks(range(len(newexp.samples)))
	plt.tight_layout()




def getcolor(idx):
	colors=['b','r','k','g','m','c','y']
	return(colors[np.mod(idx,len(colors))])



def compare2exp(exp1,exp2,unnormalize=True):
	"""
	compare 2 experiments to see differences in bacteria

	input:
	exp1,exp2 : Experiment
		the experiments to compare
	unnormalize : bool
		True (default) to revert to raw reads, False to use normalized reads instead (i.e. sum 10k)

	output:
	newexp : Experiment
		with the difference in each otu freq (exp1-exp2)
	"""
	if unnormalize:
		exp1=hs.toorigreads(exp1)
		exp2=hs.toorigreads(exp2)
	u1=hs.filterseqs(exp1,exp2.seqs,exclude=True)
	u2=hs.filterseqs(exp2,exp1.seqs,exclude=True)

	com=hs.filterseqs(exp1,exp2.seqs)
	e1=hs.filterseqs(exp1,com.seqs)
	e2=hs.filterseqs(exp2,com.seqs)

	print("unique 1 - %d, unique 2 - %d" % (len(u1.seqs),len(u2.seqs)))

	e1=hs.sortsamples(e1,'#SampleID')
	e2=hs.sortsamples(e2,'#SampleID')

	newexp=hs.copyexp(e1)
	newexp.data=e1.data-e2.data

	newexp=hs.filterminreads(newexp,0.00001,useabs=True)

	newexp=hs.sortbyfreq(newexp,logscale=False,useabs=True)
	hs.plotexp(newexp,'',uselog=False,rangeall=True,showcolorbar=True,minreads=0,showline=False)
	return newexp


def animatetimeseries(expdat,compfield,subfield,fieldname,numeric=True,seqs=[]):
	"""
	Do a bar plot animation
	for Rob change your microbes presentation on Turnbaugh diet dataset
	"""
	# import matplotlib.animation as manimation
	# FFMpegWriter = manimation.writers['ffmpeg']
	# metadata = dict(title='Movie Test', artist='Matplotlib',comment='Movie support!')
	# writer = FFMpegWriter(fps=3, metadata=metadata)

	if len(seqs)==0:
		seqs=expdat.seqs
	expdat=hs.sortsamples(expdat,subfield)
	sfvals=hs.getfieldvals(expdat,subfield,ounique=True)
	fvals=hs.getfieldvals(expdat,fieldname,ounique=True)
	compvals=hs.getfieldvals(expdat,compfield,ounique=True)
	if numeric:
		tvals=hs.tofloat(fvals)
	else:
		tvals=fvals
	svals,sidx=hs.isort(tvals)
	fvals=hs.reorder(fvals,sidx)
	fig=plt.figure()
	# with writer.saving(fig, "~/writer_test.mp4", 100):
	for cseq in seqs:
		for cval in fvals:
			outdat=np.zeros([len(compvals),len(sfvals)])
			fig=plt.figure()
			plt.clf()
			print(cval)
			cexp=hs.filtersamples(expdat,fieldname,cval)
			if len(cexp.samples)==0:
				continue
			print(len(cexp.samples))
			for idx,ccompval in enumerate(compvals):
				print(ccompval)
				ccexp=hs.filtersamples(cexp,compfield,ccompval)
				print(len(ccexp.samples))
				for idx2,csf in enumerate(sfvals):
					cccexp=hs.filtersamples(ccexp,subfield,csf)
					if len(cccexp.samples)==0:
						continue
					outdat[idx,idx2]=np.mean(cccexp.data[cccexp.seqdict[cseq]])
			plt.bar(np.arange(len(sfvals)),outdat[0,:],width=0.4)
			plt.bar(np.arange(len(sfvals))+0.5,outdat[1,:],width=0.4,color='r')
			plt.ylim([0,1000])
			plt.title(cval)
			# writer.grab_frame()


def plottaxonomybar(expdat,taxlevel=3,maxnum=8,normalize=True,sortfield=None,showlegend=True):
	"""
	plot a taxonomy bar plot of each sample at the taxonomic level taxlevel.
	input:
	expdat : Experiment
	taxlevel : int
		the taxonomic level (1=kingdom, 2=phylum etc.)
	maxnum : int
		the maximum number of unique taxonomy to show (others are grouped into 'Other')
	normalize : bool
		True (default) to show normalized data, False to use original reads
	sortfield : str
		Name of field to sort and show values in x axis, or None (default) to show sample names
	showlegend : bool
		True (default) to show the legend, False to hide it
	"""

	# collapse the taxonomy
	colors=['b','r','w','g','m','c','y','k']
	if not normalize:
		expdat=hs.toorigreads(expdat)
	expdat=hs.collapsetaxonomy(expdat,taxlevel)
	expdat=hs.sortbyfreq(expdat,reverse=True)

	if sortfield is not None:
		expdat=hs.sortsamples(expdat,sortfield)

	# clip at maxnum unique labels
	if len(expdat.seqs)>maxnum:
		expdat.data[maxnum-1,:]=np.sum(expdat.data[maxnum-1:,:],axis=0)
		expdat.tax[maxnum-1]='Other'
		expdat=hs.reorderbacteria(expdat,np.arange(maxnum))

	fig, ax = plt.subplots()
	csum=np.zeros(len(expdat.samples))
	xpos=np.arange(len(expdat.samples))
	for idx,cseq in enumerate(expdat.seqs):
		cdat=expdat.data[idx,:]
		clabel=expdat.tax[idx].split(';')
		clabel=clabel[-1]
		plt.bar(xpos,cdat,bottom=csum,color=colors[np.mod(idx,len(colors))],label=clabel)
		csum+=cdat

	# # add some text for labels, title and axes ticks
	ax.set_xlabel('number of reads')
	if showlegend:
		# need to reverse the legend (first is at the bottom)
		handles, labels = ax.get_legend_handles_labels()
		plt.legend(handles[::-1], labels[::-1], loc='upper left')
	ax.set_xticks(xpos)
	if sortfield is not None:
		xlab=hs.getfieldvals(expdat,sortfield)
	else:
		xlab=expdat.samples
	ax.set_xticklabels(xlab,rotation=90)
	plt.tight_layout()


def plotdendrogram(expdat,field=None,nosort=False,minreads=0,xwidth=0.5,zeroisnone=True,showtaxnames=True):
	"""
	plot an experiment with dendrogram for the sequences (based on sequence similarity)
	input:
	expdat : Experiment
		the experiment to plot the dendrogram for
	field : str
		the field to sort by for plotting

	output:
	newexp : Experiment
		sorted by the sequence similarity clustering
	"""
	import scipy.spatial.distance
	import scipy.cluster.hierarchy

	# calculate the sequence distance matrix
	dm=np.zeros([len(expdat.seqs),len(expdat.seqs)])
	for idx1,seq1 in enumerate(expdat.seqs):
		for idx2,seq2 in enumerate(expdat.seqs):
			cdist=sum(c1!=c2 for c1,c2 in zip(seq1,seq2))
			dm[idx1,idx2]=cdist
	dm[dm>0]+=0.5
	dm[dm>0]=np.log2(dm[dm>0])
	sdm=scipy.spatial.distance.squareform(dm)

	fig=plt.figure()
	[ax1_x, ax1_y, ax1_w, ax1_h] = [0.75,0.1,0.2,0.8]
	[axm_x, axm_y, axm_w, axm_h] = [ax1_x-xwidth,0.1,xwidth,0.8]
	# [ax1_x, ax1_y, ax1_w, ax1_h] = [0.05,0.1,0.2,0.8]
	# [axm_x, axm_y, axm_w, axm_h] = [0.25,0.1,0.7,0.8]
	ax1 = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], frame_on=True)
	linkmat = scipy.cluster.hierarchy.linkage(sdm, method='single')
	dend = scipy.cluster.hierarchy.dendrogram(linkmat,orientation='right',no_labels=True)
	tickdist=np.array([1,2,3,5,8,15,30])
	ticklabs=[]
	for ctick in tickdist:
		ticklabs.append(str(ctick))
	plt.xticks(np.log2(tickdist+0.5),ticklabs)
	idx=dend['leaves']
	print(len(idx))
	print(len(expdat.seqs))
	newexp=hs.reorderbacteria(expdat,idx)
	axm = fig.add_axes([axm_x, axm_y, axm_w, axm_h], frame_on=True)
	hs.plotexp(newexp,field,newfig=False,ptitle=None,nosort=nosort,minreads=minreads,zeroisnone=zeroisnone,showtaxnames=showtaxnames)
	return newexp



def compare2repeats(exp1,exp2,samplefield='#SampleID',sampleid=None):
	"""
	compare 2 experiments which are technical repeats of the same samples

	input:
	exp1,exp2 : Experiment
		the two experiments to compare. need to have same sampleids
	samplefield : str
		name of the field containing the sampleids for comparison (i.e. identical between the 2 experiments)
	sampleid : str
		the sample to test
	"""
	import scipy.stats

	idx1=hs.findsamples(exp1,samplefield,sampleid)
	if len(idx1)!=1:
		hs.Debug(9,'Not only 1 value for %s in exp1!!' % sampleid)
		return
	idx2=hs.findsamples(exp2,samplefield,sampleid)
	if len(idx2)!=1:
		hs.Debug(9,'Not only 1 value for %s in exp2!!' % sampleid)
		return
	idx1=idx1[0]
	idx2=idx2[0]
	dat1=exp1.data[:,idx1]
	dat2=exp2.data[:,idx2]

	minvals=[]
	maxvals=[]
	allseqs=set(exp1.seqs).union(exp2.seqs)
	print(len(allseqs))
	for cseq in allseqs:
		if cseq in exp1.seqdict:
			val1=dat1[exp1.seqdict[cseq]]
		else:
			val1=0
		if cseq in exp2.seqdict:
			val2=dat2[exp2.seqdict[cseq]]
		else:
			val2=0
		if val1==0 and val2==0:
			continue
		val1+=0.001
		val2+=0.001
		# minvals.append(np.min([val1,val2]))
		# maxvals.append(np.max([val1,val2]))
		minvals.append((val1+val2)/2)
		maxvals.append(np.abs(val1-val2)/(val1+val2))
	plt.figure()
	# plt.plot([0,5000],[0,5000],'k')
	plt.plot(minvals,maxvals,'.')
	plt.xlabel('min value')
	plt.ylabel('max value')

	# bins=np.logspace(-2,3.5,num=25)
	# binmedian,binedges,binnumber=scipy.stats.binned_statistic(minvals,maxvals,statistic='mean',bins=bins)
	# plt.plot(binedges[:-1],binmedian,'r')
	# plt.xscale('log')
	# plt.yscale('log')
	# plt.xlim([0.0001,1000])
	# plt.ylim([0.0001,1000])
	plt.xlim([-10,1000])
	plt.ylim([-0.1,1.1])


def plotstability(exp1,exp2,samplefield='#SampleID',sampleid=None,newfig=False,plotit=True,minpresratio=0.1,xpos=None):
	"""
	"""
	if sampleid is not None:
		idx1=hs.findsamples(exp1,samplefield,sampleid)
		idx2=hs.findsamples(exp2,samplefield,sampleid)
		if len(idx1)!=1:
			hs.Debug(7,'sampleid %s has %d samples in exp1' % (sampleid,len(idx1)))
			return [],[]
		if len(idx2)!=1:
			hs.Debug(7,'sampleid %s has %d samples in exp2' % (sampleid,len(idx1)))
			return [],[]
		# exp1=hs.filtersamples(exp1,samplefield,sampleid)
		# exp2=hs.filtersamples(exp2,samplefield,sampleid)

	allx=[]
	ally=[]
	if xpos is None:
		xpos=np.arange(5000)
	for minreads in xpos:
		allx.append(minreads)
		bact1=np.where(exp1.data[:,idx1]>=minreads)[0]
		if len(bact1)==0:
			ally.append(np.nan)
			continue
		numbact2=0.0
		for cseq in bact1:
			if exp1.seqs[cseq] in exp2.seqdict:
				if exp2.data[exp2.seqdict[exp1.seqs[cseq]],idx2]>minpresratio*exp1.data[cseq,idx1]:
					numbact2+=1
		ratio=numbact2/len(bact1)
		ally.append(ratio)
	if newfig:
		plt.figure()
	if plotit:
		plt.plot(allx,ally)
		plt.xscale('log')
		plt.xlim([0.1,10000])
		plt.ylim([-0.1,1.1])
	return allx,ally


def plotstabilitysummary(exp1,exp2,samplefield='#SampleID',newfig=False,minpresratio=0.0,color='r',removesingletons=False):
	xpos=np.logspace(-1,3,100)
	ally=np.zeros([0,len(xpos)])
	for cid in hs.getfieldvals(exp1,samplefield):
		cx,cy=hs.plotstability(exp1,exp2,samplefield,cid,newfig=False,plotit=False,xpos=xpos,minpresratio=minpresratio)
		ally=np.vstack([ally,cy])
	print(np.shape(ally))
	if newfig:
		plt.figure()
	cmean=np.nanmean(ally,axis=0)
	cstd=np.nanstd(ally,axis=0)/np.sqrt(np.sum(np.isfinite(ally),axis=0))
	nxpos=xpos/10000
	plt.fill_between(nxpos, cmean+cstd, cmean-cstd, alpha=0.25,color=color)
	plt.plot(nxpos,cmean,color=color)
	plt.xscale('log')
	plt.xlim([1.0/100000,0.1])
	plt.ylim([0,1.05])
	plt.xlabel('minimal sequence frequency')
	plt.ylabel('mean fraction of overlap between repeats')
