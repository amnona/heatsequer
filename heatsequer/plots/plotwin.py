#!/usr/bin/env python


"""
heatsequer heatmap plot module
"""

# amnonscript

__version__ = "0.9"


import heatsequer as hs

import numpy as np
import matplotlib as mpl
mpl.use('Qt4Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import *


def plotexp(exp,sortby=False,numeric=False,minreads=4,rangeall=False,seqdb=None,cdb=None,showline=True,ontofig=False,usegui=True,showxall=False,showcolorbar=False,ptitle=False,lowcutoff=1,uselog=True,showxlabel=True,colormap=False,colorrange=False):
	"""
	Plot an experiment
	input:
	exp - from load()
	sortby - name of mapping file field to sort by or Flase to not sort
	numeric - True if the field is numeric
	minreads - minimum number of reads per bacteria in order to show it or 0 to show all
	rangeall - True to show all frequencies in image scale, false to saturate at 10%
	seqdb - the SRBactDB database (from bactdb.load)
	cdb - the cool sequences database (from cooldb.load)
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

	output:
	newexp - the plotted experiment (sorted and filtered)
	ax - the plot axis
	"""

	hs.Debug(1,"Plot experiment %s" % exp.studyname)
	hs.Debug(1,"Commands:")
	for ccommand in exp.commands:
		hs.Debug(1,"%s" % ccommand)
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
	hs.Debug(1,"Filtering min reads. original bacteria - %d" % len(newexp.seqs))
	if minreads>0:
		newexp=hs.filterminreads(newexp,minreads,logit=uselog)
	hs.Debug(1,"New number of bacteria %d" % len(newexp.seqs))
	newexp.seqdb=seqdb
	newexp.cdb=cdb

#	ldat=ldat[:,sidx]
	ldat=newexp.data
	if uselog:
		hs.Debug(1,"Using log, cutoff at %f" % lowcutoff)
		ldat[np.where(ldat<lowcutoff)]=lowcutoff
		ldat=np.log2(ldat)
	oldparams=plt.rcParams
	mpl.rc('keymap',back='c, backspace')
	mpl.rc('keymap',forward='v')
	mpl.rc('keymap',all_axes='A')
	f=figure()
	# set the colormap to default if not supplied
	if not colormap:
		colormap=plt.rcParams['image.cmap']
	# plot the image
	if colorrange:
		hs.Debug(1,"colormap range is 0,10")
		iax=imshow(ldat,interpolation='nearest',aspect='auto',clim=colorrange,cmap=plt.get_cmap(colormap))
	elif rangeall:
		hs.Debug(1,"colormap range is all")
		iax=imshow(ldat,interpolation='nearest',aspect='auto',cmap=plt.get_cmap(colormap))
	else:
		hs.Debug(1,"colormap range is 0,10")
		iax=imshow(ldat,interpolation='nearest',aspect='auto',clim=[0,10],cmap=plt.get_cmap(colormap))

	if not ptitle:
		hs.Debug(1,"Showing filters in title")
		if (len(newexp.filters))>4:
			cfilters=[newexp.filters[0],'...',newexp.filters[-2],newexp.filters[-1]]
		else:
			cfilters=newexp.filters
		cfilters=hs.clipstrings(cfilters,30)
		ptitle='\n'.join(cfilters)
	title(ptitle,fontsize=10)

	ax=iax.get_axes()
	ax.autoscale(False)
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
			linepos.append(idx+0.5)
			labs.append(cval)
		hs.Debug(1,"number of lines is %d" % len(linepos))
		if showxlabel:
			ax.set_xticks(labpos)
			ax.set_xticklabels(labs,rotation=45,ha='right')
		for cx in linepos:
			plot([cx,cx],[-0.5,np.size(ldat,0)-0.5],'k',linewidth=2)
	else:
		hs.Debug(1,"Not showing lines")
		if showxall or len(newexp.samples)<=10:
			hs.Debug(1,"less than 10 samples, showing all sample names")
			ax.set_xticklabels(svals,rotation=90)
			ax.set_xticks(range(len(newexp.samples)))
	tight_layout()
	ax.set_ylim(-0.5,np.size(ldat,0)+0.5)

	if showcolorbar:
		hs.Debug(1,"Showing colorbar")
		cb=colorbar(ticks=list(np.log2([2,10,100,500,1000])))
		cb.ax.set_yticklabels(['<0.02%','0.1%','1%','5%','>10%'])

	# create the plot
	ax.expdat=newexp
	ax.lastselect=-1
	ax.sampline=''
	ax.ofig=f
	ax.labelson=False
	ax.labelnames=[]
	f.canvas.mpl_connect('button_press_event', onplotmouseclick)
	f.canvas.mpl_connect('key_press_event', onplotkeyclick)
#	show()
	plt.rcParams=oldparams

	# if want the ontology analysis for a given category:
	if ontofig:
		hs.Debug(1,"Ontofig is set")
		newexp.ontofigname=ontofig
	else:
		newexp.ontofigname=False

	# if we want gui, open it
	if usegui:
		hs.Debug(1,"Using the GUI window")
		import plotwingui
		guiwin = plotwingui.PlotGUIWindow(newexp)
		ax.guiwin=guiwin
		guiwin.plotfig=f
		guiwin.plotax=ax
		guiwin.show()
	else:
		ax.guiwin=False
		hs.Debug(7,'Not using gui')

	if newexp.plotmetadata:
		hs.Debug(1,"Experiment has metadata attached for plotting (%d points)" % len(newexp.plotmetadata))
		for cmet in newexp.plotmetadata:
			addplotmetadata(newexp,field=cmet[0],value=cmet[1],color=cmet[2],inverse=cmet[3],beforesample=cmet[4])
	show()
	return newexp,ax


def onplotkeyclick(event):
	if not hasattr(event,'inaxes'):
		print "boogs"
		return
	cax=event.inaxes
	if cax is None:
		print "basdoogs"
		return
	cylim=cax.get_ylim()
	cxlim=cax.get_xlim()
	cexp=cax.expdat
	if event.key=='q':
		cax.set_ylim(cylim[0], cylim[0]+(cylim[1]-cylim[0])/2)
		tight_layout()
		cax.ofig.canvas.draw()
	if event.key=='a':
		cax.set_ylim(cylim[0], cylim[0]+(cylim[1]-cylim[0])*2)
		tight_layout()
		cax.ofig.canvas.draw()
	if event.key=='down':
		cax.set_ylim(cylim[0]-(cylim[1]-cylim[0]), cylim[0])
		tight_layout()
		cax.ofig.canvas.draw()
	if event.key=='up':
		cax.set_ylim(cylim[1],cylim[1]+(cylim[1]-cylim[0]))
		tight_layout()
		cax.ofig.canvas.draw()
	if event.key=='left':
		cax.set_xlim(cxlim[0]-(cxlim[1]-cxlim[0]), cxlim[0])
		tight_layout()
		cax.ofig.canvas.draw()
	if event.key=='right':
		cax.set_xlim(cxlim[1],cxlim[1]+(cxlim[1]-cxlim[0]))
		tight_layout()
		cax.ofig.canvas.draw()
	if event.key==',':
		# select next bacteria
		cax.guiwin.clearselection()
		cax.lastselect+=1
		cax.guiwin.selectbact([cax.lastselect])
		cax.guiwin.updateinfo(cax.guiwin.csamp,cax.lastselect)
		if cexp.cdb:
			info = hs.cooldb.getseqinfo(cexp.cdb,cexp.seqs[cax.lastselect])
			if cax.guiwin:
				cax.guiwin.updatecdb(info)
			else:
				for cinfo in info:
					print (cinfo)
				sys.stdout.flush()
	if event.key=='.':
		# select prev bacteria
		cax.guiwin.clearselection()
		cax.lastselect-=1
		cax.guiwin.selectbact([cax.lastselect])
		cax.guiwin.updateinfo(cax.guiwin.csamp,cax.lastselect)
		if cexp.cdb:
			info = hs.cooldb.getseqinfo(cexp.cdb,cexp.seqs[cax.lastselect])
			if cax.guiwin:
				cax.guiwin.updatecdb(info)
			else:
				for cinfo in info:
					print (cinfo)
				sys.stdout.flush()

	if event.key=='<':
		cx=cax.guiwin.csamp
		cx=cx-1
		cy=cax.guiwin.cseq
		if cax.sampline in cax.lines:
			cax.lines.remove(cax.sampline)
		cax.sampline=cax.plot([cx,cx],[-0.5,len(cexp.sids)-0.5],':w')[0]
		cax.guiwin.updateinfo(cx,cy)
		cax.ofig.canvas.draw()
	if event.key=='>':
		cx=cax.guiwin.csamp
		cx=cx+1
		cy=cax.guiwin.cseq
		if cax.sampline in cax.lines:
			cax.lines.remove(cax.sampline)
		cax.sampline=cax.plot([cx,cx],[-0.5,len(cexp.sids)-0.5],':w')[0]
		cax.guiwin.updateinfo(cx,cy)
		cax.ofig.canvas.draw()
	# show taxonomies
	if event.key=='h':
		showtaxonomies(cexp,cax)
	# nice taxonomies (genus+species)
	if event.key=='n':
		labs=[]
		for ctax in cexp.tax:
			cstr=ctax.split(';')
			labs.append(cstr[-2]+';'+cstr[-1])
		cax.set_yticks(np.array(range(len(cexp.seqs))))
		cax.set_yticklabels(labs)
		cax.set_ylim(cylim[0], cylim[1])
		tight_layout()
		cax.ofig.canvas.draw()


def getlabelnames(cexp,showdb=True,showcontam=True):
	"""
	Get the sequence label list for an experiment

	input:
	"""
	pass


def showtaxonomies(cexp,cax,show=True,showdb=True,showcontam=True,maxtax=250):
	"""
	show the y-lables (taxonomies) for the plot window

	input:
	cexp : Experiment
	cax : axis (matplotlib)
		the plot window axis
	show : bool
		True (default) to show the labels, False to remove them
	showdb : bool
		True (default) to add '*' to sequences in the cooldb database, False to not add '*'
	showcontam : bool
		True (default) to show suspected contamination bacteria in red, False to not color separately
	maxtax : int
		Maximal number of taxonomies to show (to prevent slow repsonse when looking at big experiment) or 0 to show all
	"""

	cylim=cax.get_ylim()
	cxlim=cax.get_xlim()

	# check if we show too many bacteria - don't show taxonomy labels
	if maxtax>0:
		if cylim[1]-cylim[0]>maxtax:
			cax.set_yticklabels([])
			hs.Debug(7,'Too many bacteria - zoom in to show labels')
			return

	contamlist=[]
	pathogenlist=[]
	labs=hs.clipstrings(cexp.tax,25,reverse=True)
	if showdb or contam:
		if cexp.cdb:
			for idx,cseq in enumerate(cexp.seqs):
				info=hs.cooldb.getseqinfo(cexp.cdb,cseq)
				if len(info)>0:
					for cinfo in info:
						if "ontamination" in cinfo:
							contamlist.append(idx)
						if "patric" in cinfo:
							pathogenlist.append(idx)
					if showdb:
						labs[idx]+='*'
	cax.set_yticks(np.array(range(len(cexp.seqs))))
	cax.tick_params(axis='y', which='major', labelsize=8)
	cax.set_yticklabels(labs)
	if showcontam:
		for idx,clab in enumerate(cax.get_yticklabels()):
			if idx in contamlist:
				clab.set_color("red")
			if idx in pathogenlist:
				clab.set_color("blue")
	cax.set_ylim(cylim[0], cylim[1])
	cax.set_xlim(cxlim[0], cxlim[1])
	tight_layout()
	cax.ofig.canvas.draw()


def onplotmouseclick(event):
	hs.Debug(0,event.xdata,event.ydata)
	hs.Debug(0,event.key)
	cexp=event.inaxes.expdat
	ax=event.inaxes
	rx=int(round(event.xdata))
	ry=int(round(event.ydata))
	cseq=cexp.seqs[ry]
	if ax.guiwin:
		# plot the sample line (vertical)
		if ax.sampline in ax.lines:
			ax.lines.remove(ax.sampline)
		ax.sampline=ax.plot([rx,rx],[-0.5,len(cexp.sids)-0.5],':w')[0]
		if event.key:
			if 'super' not in event.key:
				ax.guiwin.clearselection()
			if 'shift' in event.key:
				p1=min(ry,ax.lastselect)
				p2=max(ry,ax.lastselect)
				ax.guiwin.selectbact(range(p1,p2+1),flip=False)
			else:
				ax.guiwin.selectbact([ry])
				ax.lastselect=ry
		else:
			ax.guiwin.clearselection()
			if ax.lastselect!=ry:
				ax.guiwin.selectbact([ry])
				ax.lastselect=ry
			else:
				ax.lastselect=-1
		ax.guiwin.updateinfo(rx,ry)
#			bactdb.GetSeqInfo(cexp.seqdb,cseq)
	else:
		print('*********************************')
		print('Bacteria: %s - %s' % (cexp.sids[ry],cexp.tax[ry]))
		print(cseq)
		print('Sample: %s' % cexp.samples[rx])
		sys.stdout.flush()
	if cexp.cdb:
		info = hs.cooldb.getseqinfo(cexp.cdb,cexp.seqs[ry])
		if ax.guiwin:
			ax.guiwin.updatecdb(info)
		else:
			for cinfo in info:
				print (cinfo)
			sys.stdout.flush()




def addplotmetadata(expdat,field,value=False,inverse=False,color='g',ax=False,beforesample=True,partial=False):
	"""
	plot lines on an experiment plot from plotexp. NOTE: need to use with the output of plotexp!!!
	input:
	expdat : Experiment
	field : string
		name of the metadata field to use
	value : string
		the value that when present we will plot, or False to plot whenever not empty
	inverse : bool
		inverse the logic of value
	color : string
		the color to plot (matplotlib plot style - i.e. 'r:')
	beforesample : bool
		True if it happened between the prev. and current sample, False if in the middle of the sample
	partial : bool
		True to allow substring match for field/value, False for exact match
	"""

	if not ax:
		ax=plt.gca()
	if beforesample:
		offset=-0.5
	else:
		offset=0
	for idx,csamp in enumerate(expdat.samples):
		plotit=False
		if value:
			if partial:
				if value in expdat.smap[csamp][field]:
					plotit=True
			else:
				if expdat.smap[csamp][field]==value:
					plotit=True
		else:
			if expdat.smap[csamp][field]!='':
				plotit=True
		if inverse:
			plotit=not plotit
		if plotit:
			hs.Debug(1,'Plot line %d',idx)
			plt.plot([idx+offset,idx+offset],[-0.5,len(expdat.sids)-0.5],color)
