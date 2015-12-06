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


def plotexp(exp,sortby=False,numeric=False,minreads=4,rangeall=False,seqdb=None,cdb=None,showline=True,ontofig=False,usegui=True,showxall=False,showcolorbar=False,ptitle=False,lowcutoff=1,uselog=True):
	"""
	Plot an experiment
	input:
	exp - from load()
	sortby - name of mapping file field to sort by or Flase to not sort
	numeric - True if the field is numeric
	minreads - minimum number of reads per bacteria in order to show it
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

	output:
	newexp - the plotted experiment (sorted and filtered)
	ax - the plot axis
	"""

	vals=[]
	if sortby:
		for csamp in exp.samples:
			vals.append(exp.smap[csamp][sortby])
		if numeric:
			vals=hs.tofloat(vals)
		svals,sidx=hs.isort(vals)
		newexp=hs.reordersamples(exp,sidx)
	else:
		svals=hs.getfieldvals(exp,'#SampleID')
		newexp=hs.copyexp(exp)
	newexp=hs.filterminreads(newexp,minreads)
	newexp.seqdb=seqdb
	newexp.cdb=cdb

#	ldat=ldat[:,sidx]
	ldat=newexp.data
	if uselog:
		ldat[np.where(ldat<1)]=1
		ldat=np.log2(ldat)
	oldparams=plt.rcParams
	mpl.rc('keymap',back='c, backspace')
	mpl.rc('keymap',forward='v')
	mpl.rc('keymap',all_axes='A')
	f=figure()
	if rangeall:
		iax=imshow(ldat,interpolation='nearest',aspect='auto')
	else:
		iax=imshow(ldat,interpolation='nearest',aspect='auto',clim=[0,10])

	if not ptitle:
		if (len(newexp.filters))>4:
			cfilters=[newexp.filters[0],'...',newexp.filters[-2],newexp.filters[-1]]
		else:
			cfilters=newexp.filters
		ptitle='\n'.join(cfilters)
	title(ptitle,fontsize=10)

	ax=iax.get_axes()
	ax.autoscale(False)
	if showline:
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
		ax.set_xticks(labpos)
		ax.set_xticklabels(labs,rotation=45,ha='right')
		for cx in linepos:
			plot([cx,cx],[-0.5,np.size(ldat,0)-0.5],'k',linewidth=2)
	else:
		if showxall or len(newexp.samples)<=10:
			ax.set_xticklabels(svals,rotation=90)
			ax.set_xticks(range(len(newexp.samples)))
	tight_layout()
	ax.set_ylim(-0.5,np.size(ldat,0)+0.5)

	if showcolorbar:
		cb=colorbar(ticks=list(np.log2([2,10,100,500,1000])))
		cb.ax.set_yticklabels(['<0.02%','0.1%','1%','5%','>10%'])

	# create the plot
	ax.expdat=newexp
	ax.lastselect=-1
	ax.sampline=''
	ax.ofig=f
	f.canvas.mpl_connect('button_press_event', onplotmouseclick)
	f.canvas.mpl_connect('key_press_event', onplotkeyclick)
#	show()
	plt.rcParams=oldparams

	# if want the ontology analysis for a given category:
	if ontofig:
		newexp.ontofigname=ontofig
	else:
		newexp.ontofigname=False

	# if we want gui, open it
	if usegui:
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

	if event.key=='left':
		cx=cax.guiwin.csamp
		cx=cx-1
		cy=cax.guiwin.cseq
		if cax.sampline in cax.lines:
			cax.lines.remove(cax.sampline)
		cax.sampline=cax.plot([cx,cx],[-0.5,len(cexp.sids)-0.5],':w')[0]
		cax.guiwin.updateinfo(cx,cy)
		cax.ofig.canvas.draw()
	if event.key=='right':
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
		contamlist=[]
		cax.set_yticks(np.array(range(len(cexp.seqs))))
		labs=hs.clipstrings(cexp.tax,25,reverse=True)
		if cexp.cdb:
			for idx,cseq in enumerate(cexp.seqs):
				info=hs.cooldb.getseqinfo(cexp.cdb,cseq)
				if len(info)>0:
					for cinfo in info:
						if "ontamination" in cinfo:
							contamlist.append(idx)
					labs[idx]+='*'
		cax.tick_params(axis='y', which='major', labelsize=8)
		cax.set_yticklabels(labs)
		for idx,clab in enumerate(cax.get_yticklabels()):
			if idx in contamlist:
				clab.set_color("red")
		cax.set_ylim(cylim[0], cylim[1])
		cax.set_xlim(cxlim[0], cxlim[1])
		tight_layout()
		cax.ofig.canvas.draw()
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
			if not 'super' in event.key:
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
