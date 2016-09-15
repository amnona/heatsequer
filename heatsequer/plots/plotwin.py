#!/usr/bin/env python


"""
heatsequer heatmap plot module
"""

# amnonscript

from __future__ import absolute_import

__version__ = "0.9"


import heatsequer as hs

import sys
import numpy as np
import matplotlib as mpl
#mpl.use('Qt4Agg')
import matplotlib.pyplot as plt
#from matplotlib.pyplot import *


def plotexp(exp,sortby=False,numeric=False,minreads=4,rangeall=False,seqdb=None,cdb=None,showline=True,ontofig=False,usegui=True,showxall=False,showcolorbar=False,ptitle=False,lowcutoff=1,uselog=True,showxlabel=True,colormap=False,colorrange=False,linewidth=2,subline='',showhline=True,newfig=True):
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
	showhline : bool
		True (default) to plot the horizontal lines listed in exp.hlines. False to not plot them
	newfig : bool
		True (default) to open figure in new window, False to use current

	output:
	newexp - the plotted experiment (sorted and filtered)
	ax - the plot axis
	"""

	hs.Debug(1,"Plot experiment %s" % exp.studyname)
	hs.Debug(1,"Commands:")
	for ccommand in exp.commands:
		hs.Debug(1,"%s" % ccommand)

	if exp.sparse:
		hs.Debug(9,'Sparse matrix - converting to dense')
		exp=hs.copyexp(exp,todense=True)

	vals=[]
	if cdb is None:
		cdb=hs.cdb
	if seqdb is None:
		seqdb=hs.bdb
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
	newexp.scdb=hs.scdb

	# if usegui:
	# 	hs.Debug(1,"Using the GUI window")
	# 	import heatsequer.plots.plotwingui
	# 	from PyQt4 import QtGui

	# 	app = QtGui.QApplication(sys.argv)
	# 	guiwin = heatsequer.plots.plotwingui.PlotGUIWindow(newexp)

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
	if newfig:
		f=plt.figure(tight_layout=True)
	else:
		f=plt.gcf()
	# set the colormap to default if not supplied
	if not colormap:
		colormap=plt.rcParams['image.cmap']
	# plot the image
	if colorrange:
		hs.Debug(1,"colormap range is 0,10")
		iax=plt.imshow(ldat,interpolation='nearest',aspect='auto',clim=colorrange,cmap=plt.get_cmap(colormap))
	elif rangeall:
		hs.Debug(1,"colormap range is all")
		iax=plt.imshow(ldat,interpolation='nearest',aspect='auto',cmap=plt.get_cmap(colormap))
	else:
		hs.Debug(1,"colormap range is 0,10")
		iax=plt.imshow(ldat,interpolation='nearest',aspect='auto',clim=[0,10],cmap=plt.get_cmap(colormap))

	if not ptitle:
		hs.Debug(1,"Showing filters in title")
		if (len(newexp.filters))>4:
			cfilters=[newexp.filters[0],'...',newexp.filters[-2],newexp.filters[-1]]
		else:
			cfilters=newexp.filters
		cfilters=hs.clipstrings(cfilters,30)
		ptitle='\n'.join(cfilters)
	plt.title(ptitle,fontsize=10)

	ax=iax.get_axes()
	ax.autoscale(False)

	# plot the sublines (smaller category lines)
	if subline:
		slval=hs.getfieldvals(newexp,subline)
		prevval=slval[0]
		for idx,cval in enumerate(slval):
			if cval!=prevval:
				xpos=idx-0.5
				plt.plot([xpos,xpos],[-0.5,np.size(ldat,0)-0.5],'w:')
				prevval=cval

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
			plt.plot([cx,cx],[-0.5,np.size(ldat,0)-0.5],'k',linewidth=linewidth)
			plt.plot([cx,cx],[-0.5,np.size(ldat,0)-0.5],'w:',linewidth=linewidth)
	else:
		hs.Debug(1,"Not showing lines")
		if showxall or len(newexp.samples)<=10:
			hs.Debug(1,"less than 10 samples, showing all sample names")
			ax.set_xticklabels(svals,rotation=90)
			ax.set_xticks(range(len(newexp.samples)))
	# f.tight_layout()
	ax.set_ylim(-0.5,np.size(ldat,0)+0.5)

	if showcolorbar:
		hs.Debug(1,"Showing colorbar")
		cb=plt.colorbar(ticks=list(np.log2([2,10,100,500,1000])))
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
		import heatsequer.plots.plotwingui
#		from PyQt4 import QtGui

#		app = QtGui.QApplication(sys.argv)
		guiwin = heatsequer.plots.plotwingui.PlotGUIWindow(newexp)
		from heatsequer.plots import plotwingui
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
	if showhline:
		if newexp.hlines:
			for cpos in newexp.hlines:
				plt.plot([0,np.shape(newexp.data)[1]],[cpos-0.5,cpos-0.5],'g')
	plt.show()

#	if usegui:
#		app.exec_()

	return newexp,ax


def onplotkeyclick(event):
	if not hasattr(event,'inaxes'):
		print("boogs")
		return
	cax=event.inaxes
	if cax is None:
		print("basdoogs")
		return
	cylim=cax.get_ylim()
	cxlim=cax.get_xlim()
	cexp=cax.expdat
	if event.key=='C':
		plt.figure()
		plt.plot([0,1],[0,1],'k')
	if event.key=='q':
		cax.set_ylim(cylim[0], cylim[0]+(cylim[1]-cylim[0])/2)
		plt.tight_layout()
		cax.ofig.canvas.draw()
	if event.key=='a':
		cax.set_ylim(cylim[0], cylim[0]+(cylim[1]-cylim[0])*2)
		plt.tight_layout()
		cax.ofig.canvas.draw()
	if event.key=='Q':
		cax.set_xlim(cxlim[0], cxlim[0]+(cxlim[1]-cxlim[0])/2)
		plt.tight_layout()
		cax.ofig.canvas.draw()
	if event.key=='A':
		cax.set_xlim(cxlim[0], cxlim[0]+(cxlim[1]-cxlim[0])*2)
		plt.tight_layout()
		cax.ofig.canvas.draw()
	if event.key=='down':
		if cylim[0]>0:
			cax.set_ylim(cylim[0]-(cylim[1]-cylim[0]), cylim[0])
			plt.tight_layout()
			cax.ofig.canvas.draw()
	if event.key=='up':
		if cylim[1]<len(cexp.seqs):
			cax.set_ylim(cylim[1],cylim[1]+(cylim[1]-cylim[0]))
			plt.tight_layout()
			cax.ofig.canvas.draw()
	if event.key=='left':
		if cxlim[0]>0:
			cax.set_xlim(cxlim[0]-(cxlim[1]-cxlim[0]), cxlim[0])
			plt.tight_layout()
			cax.ofig.canvas.draw()
	if event.key=='right':
		if cxlim[1]<len(cexp.samples)-1:
			cax.set_xlim(cxlim[1],cxlim[1]+(cxlim[1]-cxlim[0]))
			plt.tight_layout()
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
					print(cinfo)
				sys.stdout.flush()
		if cexp.scdb:
			info = hs.supercooldb.getcurationstrings(cexp.scdb,cexp.seqs[cax.lastselect])
			if cax.guiwin:
				cax.guiwin.addtocdblist(info)
			else:
				for cinfo in info:
					print(cinfo)
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
					print(cinfo)
				sys.stdout.flush()
		if cexp.scdb:
			info = hs.supercooldb.getcurationstrings(cexp.scdb,cexp.seqs[cax.lastselect])
			if cax.guiwin:
				cax.guiwin.addtocdblist(info)
			else:
				for cinfo in info:
					print(cinfo)
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
	if event.key=='H':
		showsampleids(cexp,cax)
	# nice taxonomies (genus+species)
	if event.key=='n':
		labs=[]
		for ctax in cexp.tax:
			cstr=ctax.split(';')
			labs.append(cstr[-2]+';'+cstr[-1])
		cax.set_yticks(np.array(range(len(cexp.seqs))))
		cax.set_yticklabels(labs)
		cax.set_ylim(cylim[0], cylim[1])
		plt.tight_layout()
		cax.ofig.canvas.draw()
	# hide the cursor
	if event.key=='!':
		cax.guiwin.clearselection()
		cax.lines.remove(cax.sampline)
		cax.ofig.canvas.draw()


def getlabelnames(cexp,showdb=True,showcontam=True):
	"""
	Get the sequence label list for an experiment

	input:
	"""
	pass


def showtaxonomies(cexp,cax,show=True,showdb=True,showcontam=True,maxtax=250,nicenames=True):
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
	nicenames :  bool
		True to show genus+species names, False to show last 25 chars
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
	if nicenames:
		labs=hs.getnicetaxnames(cexp.tax)
	else:
		labs=hs.clipstrings(cexp.tax,25,reverse=True)

	if showdb or showcontam:
		if cexp.cdb:
			for idx,cseq in enumerate(cexp.seqs):
				info=hs.cooldb.getseqinfo(cexp.cdb,cseq)
				if len(info)>0:
					for cinfo in info:
						if "patric" in cinfo:
							pathogenlist.append(idx)
						if "ontamination" in cinfo:
							contamlist.append(idx)
					if showdb:
						labs[idx]+='*'
	cax.set_yticks(np.array(range(len(cexp.seqs))))
	cax.tick_params(axis='y', which='major', labelsize=8)
	cax.set_yticklabels(labs)
	if showcontam:
		for idx,clab in enumerate(cax.get_yticklabels()):
			if idx in pathogenlist:
				clab.set_color("blue")
			if idx in contamlist:
				clab.set_color("red")
	cax.set_ylim(cylim[0], cylim[1])
	cax.set_xlim(cxlim[0], cxlim[1])
	plt.tight_layout()
	cax.ofig.canvas.draw()


def showsampleids(cexp,cax,show=True,maxsamples=25):
	"""
	show the y-lables (taxonomies) for the plot window

	input:
	cexp : Experiment
	cax : axis (matplotlib)
		the plot window axis
	show : bool
		True (default) to show the labels, False to remove them
	maxsamples : int
		Maximal number of samples to show (to prevent slow repsonse when looking at big experiment) or 0 to show all
	"""

	cylim=cax.get_ylim()
	cxlim=cax.get_xlim()

	# check if we show too many bacteria - don't show taxonomy labels
	if maxsamples>0:
		if cxlim[1]-cxlim[0]>maxsamples:
			cax.set_xticklabels([])
			hs.Debug(7,'Too many samples - zoom in (Q/A) to show labels')
			return

	labs=hs.clipstrings(cexp.samples,25,reverse=True)

	cax.set_xticks(np.array(range(len(labs))))
	cax.tick_params(axis='x', which='major', labelsize=8)
	cax.set_xticklabels(labs,rotation=45)

# 	nax = cax.twiny()
# 	nax.set_xticks(np.array(range(len(labs))))
# 	nax.tick_params(axis='x', which='major', labelsize=8)
# 	nax.set_xticklabels(labs,rotation=45)
# #	nax.set_ylim(cylim[0], cylim[1])
# 	nax.set_xlim(cxlim[0], cxlim[1])
# 	nax.expdat=cax.expdat
# 	nax.ofig=cax.ofig
# 	nax.guiwin=cax.guiwin

	cax.set_ylim(cylim[0], cylim[1])
	cax.set_xlim(cxlim[0], cxlim[1])
	plt.sca(cax)

	plt.tight_layout()
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
				print(cinfo)
			sys.stdout.flush()
	if cexp.scdb:
		info = hs.supercooldb.getcurationstrings(cexp.scdb,cexp.seqs[ry])
		if ax.guiwin:
			ax.guiwin.addtocdblist(info)
		else:
			for cinfo in info:
				print(cinfo)
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


def getnicetaxnames(tax,separator=';',level=5,maxnum=2):
	"""
	get a nice string of the taxonomy names (genus+species if available)
	input:
	tax : list of str
		the full taxonomy names (separated by separator)
	separator : str
		the separator between the levels (default is ';')
	level : int
		the default taxonomic level to show (0=bacteria, 5=genus)

	output:
	nicetax : list of str
		the nice names taxonomy (genus+species if there, if no genus - highest level available)
	"""

	hs.Debug(1,'getting nice names for %d taxonomies' % len(tax))
	nicetax=[]
	for ctax in tax:
		cntax='NA'
		stax=ctax.split(separator)
		# remove the trailing empty taxonomies
		while stax[-1]=='':
			del stax[-1]

		if len(stax)<=level:
			cntax=stax[-1]
		else:
			endpos=min(level+maxnum,len(stax))
			cntax=separator.join(stax[level:endpos])
		nicetax.append(cntax)
	hs.Debug(1,'created %d nice names' % len(nicetax))
	return nicetax
