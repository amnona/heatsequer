#!/usr/bin/env python

"""
amnonscript
analysis.py
Analyze biom tables

# to setup databases:
import bactdb
import cooldb

db=bactdb.initdb()
cdb=cooldb.loaddb()

exp=analysis.load(tablename='/Users/amnon/Projects/shorttime/short.clean.single.table.txt',mapname='/Users/amnon/Projects/shorttime/map.txt')
analysis.plotexp(exp)
"""

__version__ = "0.2"

import amnonutils as au

import os.path
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Qt4Agg')
import matplotlib.pyplot as plt
import biom
from matplotlib.pyplot import *
#from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
#from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
#from PyQt4 import QtGui, QtCore, uic
import csv
from scipy import cluster
from scipy import spatial
from scipy import stats
from sklearn.preprocessing import scale
from collections import defaultdict
import copy
# import Tkinter
# for debugging - use XXX()
from pdb import set_trace as XXX
import bactdb
import cooldb
#import plotgui

class experiment:
	def __init__(self):
		# the data matrix (non sparse)
		self.data=[]

		# the sample dictionary (double hash - sampleid and then mapping file field)
		self.smap=[]

		# name of all the fields in the mapping data
		self.fields=[]

		# list of sampleids ordered according to the data matrix
		self.samples=[]

		# the sequences in the table
		self.seqs=[]

		# dictionary holding all the sequences and their position (for fast lookup)
		self.seqdict={}

		# taxonomies
		self.tax=[]

		# the hashed ids for the sequences
		self.sids=[]

		# the original name for each otu (from the biom table)
		self.origotunames=[]

		# original table name
		self.tablefilename=''

		# original mapping file
		self.mapfilename=''

		# name of the study (or the table file name without path)
		self.studyname=''

		# number of reads for each sample in the biom table
		self.origreads=[]

		# the history of actions performed
		self.filters=[]

		# and the command list
		self.commands=[]

		# the complete sequence database
		self.seqdb=None

		# the cool sequences (manually curated) database
		self.cdb=None

		# the list of annotations to add to plot (for addplotmetadata)
		self.plotmetadata=[]

def load(tablename, mapname='map.txt', taxfile='', nameisseq=True,addsname='',studyname=False):
	"""
	Load an experiment - a biom table and a mapping file
	input:
	tablename - the name of the biom table file
	mapname - name of the mapping file
	taxfile - empty ('') to load taxonomy from biom table, non-empty to load
	from rdp output file (web)
	nameisseq - False to keep otu name as sid without hashing it, True to treat otuid as sequence
	addsname - a string to add to each sample name (or empty to not add)
	studyname - Flase to assign from table file name, otherwise string to store as study name
	output:
	an experiment class for the current experiment
	"""

	au.Debug(6,'Loading biom table')
	# load the biom table
	table = biom.load_table(tablename)

	au.Debug(6,'Loading mapping file')
	# load the mapping file
	mapf = open(mapname, 'rU')
	reader = csv.DictReader(mapf, delimiter='\t')
	fields = reader.fieldnames
	smap = {}
	mapsamples = []
	for cline in reader:
		cid = cline['#SampleID']
		smap[cid] = cline
		mapsamples.append(cid)
	mapf.close()
	au.Debug(6,'number of samples in map is %d' % len(mapsamples))

	if addsname!='':
		idtable={}
		ids=table.ids(axis='sample')
		for cid in ids:
			idtable[cid]=addsname+cid
		table=table.update_ids(idtable,axis='sample')

	tablesamples = table.ids(axis='sample')
	au.Debug(6,'number of samples in table is %d' % len(tablesamples))
	removelist=[]
	for cid in tablesamples:
		if cid not in mapsamples:
			removelist.append(cid)
			au.Debug(6,'Table sample %s not found in mapping file' % cid)
	au.Debug(6,'removing %s samples' % len(removelist))
	table=table.filter(removelist,axis='sample',invert=True)

	# get the sampleIDs not in mapping file
	datsamples = {}
	tablesamples = table.ids(axis='sample')

	au.Debug(6,'deleted. number of samples in table is now %d' % len(tablesamples))
	for cpos, cid in enumerate(tablesamples):
		if cid in mapsamples:
			datsamples[cid] = cpos
		else:
			au.Debug(10,'Errro! should have been removed table sample %s not found in mapping file' % cid)

	removemap=[]
	# get sampleids not in table
	for idx,cmap in enumerate(mapsamples):
		if cmap not in tablesamples:
			removemap.append(idx)
			try:
				del smap[cmap]
			except:
				au.Debug(8,'Duplicate SampleID %s in mapping file' % cmap)
	if len(removemap)>0:
		au.Debug(7,'need to remove %d samples from mapping file' % len(removemap))
	mapsamples=au.delete(mapsamples,removemap)


	allsamples = set(mapsamples)
	allsamples = allsamples.intersection(tablesamples)
	tableseqs = table.ids(axis='observation')

	sids = []
	tax = []
	osnames=[]
	for cid in tableseqs:
		osnames.append(cid)
		if nameisseq:
			sids.append(au.mlhash(cid, emod=10000000))
		else:
			sids.append(cid)
		md = table.metadata(cid, axis='observation')
		if md:
			if 'taxonomy' in md:
				ctax = md['taxonomy']
				if not isinstance(ctax,str):
					ctax=[x[3:] if x[2]=='_' else x for x in ctax]
					ctax = ';'.join(ctax)
				tax.append(ctax)
			else:
				tax.append('unknown')
		else:
			tax.append('unknown')

	if not studyname:
		studyname=os.path.basename(tablename)

	exp=experiment()
	exp.data=table.matrix_data.todense().A
	exp.smap=smap
	exp.samples=tablesamples
	exp.seqs=tableseqs
	for idx,cseq in enumerate(exp.seqs):
		exp.seqdict[cseq]=idx
	exp.sids=sids
	exp.origotunames=osnames
	exp.tax=tax
	exp.tablefilename=tablename
	exp.studyname=studyname
	exp.mapfilename=tablename
	exp.filters = [tablename]
	exp.fields = fields
	colsum=np.sum(exp.data,axis=0,keepdims=False)
	exp.origreads=list(colsum)
	# add the original number of reads as a field to the experiment
	exp.fields.append('origReads')
	for idx,csamp in enumerate(exp.samples):
		exp.smap[csamp]['origReads']=str(exp.origreads[idx])

	# normalize samples to 10k reads per samples
	colsum=np.sum(exp.data,axis=0,keepdims=True)
	okreads=np.where(colsum>0)
	if np.size(colsum)-np.size(okreads[1])>0:
		print("Samples with 0 reads: %d" % (np.size(colsum)-np.size(okreads[1])))
		exp=reordersamples(exp,okreads[1])
		colsum=np.sum(exp.data,axis=0,keepdims=True)

	exp.data=10000*exp.data/colsum

	exp.commands.append("load('"+tablename+"',mapname='"+mapname+")")
	exp.filters.append('loaded '+tablename)
	exp=sortbacteria(exp)
	return(exp)


def reordersamples(exp,newpos,inplace=False):
	"""
	reorder the samples of the experiment
	input:
	exp - the experiment
	newpos - array - the new positions (can skip positions to delete them)
	output:
	newexp - the new experiment
	"""

	if inplace:
		newexp=exp
	else:
		newexp=copy.deepcopy(exp)
	newexp.data=newexp.data[:,newpos]
	newexp.samples=au.reorder(newexp.samples,newpos)
	newexp.origreads=au.reorder(newexp.origreads,newpos)
	return newexp




def plotexp(exp,sortby=False,numeric=False,minreads=4,rangeall=False,seqdb=None,cdb=None,showline=True,ontofig=False,usegui=True,showxall=False,showcolorbar=False,ptitle=False):
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

	output:
	newexp - the plotted experiment (sorted and filtered)
	ax - the plot axis
	"""

	vals=[]
	if sortby:
		for csamp in exp.samples:
			vals.append(exp.smap[csamp][sortby])
		if numeric:
			vals=au.tofloat(vals)
		svals,sidx=au.isort(vals)
		newexp=reordersamples(exp,sidx)
	else:
		svals=getfieldvals(exp,'#SampleID')
		newexp=copy.deepcopy(exp)
	newexp=filterminreads(newexp,minreads)
	newexp.seqdb=seqdb
	newexp.cdb=cdb

#	ldat=ldat[:,sidx]
	ldat=newexp.data
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

	print('lala')
	# if want gui, open it
	if usegui:
		import plotgui
		guiwin = plotgui.PlotGUIWindow(newexp)
		ax.guiwin=guiwin
		guiwin.plotfig=f
		guiwin.plotax=ax
		guiwin.show()
	else:
		ax.guiwin=False
		au.Debug(7,'Not using gui')

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
		cax.set_yticks(np.array(range(len(cexp.seqs))))
		labs=au.clipstrings(cexp.tax,25,reverse=True)
		if cexp.cdb:
			for idx,cseq in enumerate(cexp.seqs):
				info=cooldb.getseqinfo(cexp.cdb,cseq)
				if len(info)>0:
					labs[idx]+='*'
		cax.tick_params(axis='y', which='major', labelsize=8)
		cax.set_yticklabels(labs)
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
	au.Debug(0,event.xdata,event.ydata)
	au.Debug(0,event.key)
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
		info = cooldb.getseqinfo(cexp.cdb,cexp.seqs[ry])
		if ax.guiwin:
			ax.guiwin.updatecdb(info)
		else:
			for cinfo in info:
				print (cinfo)
			sys.stdout.flush()


def filterminreads(exp,minreads):
	"""
	filter away all bacteria that contain less than minreads in all samples together (out of 10k/samples)
	input:
	exp - the experiment
	minreads - the minimum number of reads total for all samples (and out of 10k/sample) for a bacteria to be kept
	output:
	newexp - the filtered experiment
	"""
	numreads=np.sum(exp.data,axis=1)
	keep=np.where(numreads>=minreads)
	newexp=reorderbacteria(exp,keep[0])
	newexp.filters.append('filter min reads %d' % minreads)
	au.Debug(6,'%d Bacteria left' % len(newexp.sids))
	return newexp


def filterpresence(expdat,frac):
	"""
	filter away bacteria present in less than frac of the samples
	input:
	expdat
	frac - the minimal fraction of samples to appear in for a beacteria to be kept

	output:
	newexp - the filtered experiment
	"""
	fracreads=(np.sum(expdat.data>0,axis=1)+0.0)/len(expdat.samples)
	keep=np.where(fracreads>=frac)
	newexp=reorderbacteria(expdat,keep[0])
	newexp.filters.append('filter presence %f' % frac)
	au.Debug(6,'%d Bacteria left' % len(newexp.sids))
	return newexp

def filtermean(expdat,meanval):
	"""
	filter keeping bacteria with a mean frequency >= meanval
	input:
	exp - the experiment
	meanval - the minimum mean reads of per sample (and out of 10k/sample) for a bacteria to be kept
	output:
	newexp - the filtered experiment
	"""
	meanreads=np.mean(expdat.data,axis=1)
	keep=np.where(meanreads>=meanval)
	newexp=reorderbacteria(expdat,keep[0])
	newexp.filters.append('filter mean reads %f' % meanval)
	au.Debug(6,'%d Bacteria left' % len(newexp.sids))
	return newexp


def filterorigreads(exp,minreads):
	"""
	filter away all samples that contained originally less than minreads
	input:
	exp - the experiment
	minreads - the minimum number of reads of the sample in the biom table to filter if less
	output:
	newexp - the filtered experiment
	"""
	numreads=np.array(exp.origreads)
	keep=np.where(numreads>=minreads)
	newexp=reordersamples(exp,keep[0])
	au.Debug(6,'%d Samples left' % len(newexp.samples))
	return newexp


def reorderbacteria(exp,order,inplace=False):
	"""
	reorder the bacteria in an experiment (can delete if bacteria not in new order)
	input:
	exp - the experiment
	order - the new order
	output:
	newexp
	"""
	if inplace:
		newexp=exp
	else:
		newexp=copy.deepcopy(exp)
	newexp.data=newexp.data[order,:]
	newexp.seqs=au.reorder(newexp.seqs,order)
	newexp.seqdict={}
	for idx,cseq in enumerate(newexp.seqs):
		newexp.seqdict[cseq]=idx
	newexp.tax=au.reorder(newexp.tax,order)
	newexp.sids=au.reorder(newexp.sids,order)
	return newexp


def filtersamples(exp,field,filt,exact=True,exclude=False,numexpression=False):
	"""
	filter samples in experiment according to field
	input:
	exp
	field - name of the field
	filt - the string to filter
	exact - True for exact match, False for substring
	exclude - False to keep only matching samples, True to exclude matching samples
	numexpression - True if val is a python expression, False if just a value. For an expression assume value is the beggining of the line (i.e. '<=5')
	"""
	keep=[]
	for cidx,csamp in enumerate(exp.samples):
		keepit=False
		if numexpression:
			cval=exp.smap[csamp][field]
			if eval(cval+filt):
				keepit=True
		elif exact:
			if exp.smap[csamp][field]==filt:
				keepit=True
		else:
			if filt in exp.smap[csamp][field]:
				keepit=True
		# if exclude reverse the decision
		if exclude:
			keepit=not keepit
		if keepit:
			keep.append(cidx)
	newexp=reordersamples(exp,keep)
	fstr="filter data %s in %s " % (filt,field)
	if exact:
		fstr=fstr+"(exact)"
	else:
		fstr=fstr+"(substr)"
	if exclude:
		fstr+=" (exclude)"
	newexp.filters.append(fstr)
	newexp.commands.append("filterdata(exp,'%s','%s',exact=%s,exclude=%s)" % (field,filt,exact,exclude))
	au.Debug(6,'%d Samples left' % len(newexp.samples))
	return newexp


def getfieldvals(exp,field):
	"""
	get a list of the field values in all samples
	"""
	vals=[]
	for cid in exp.samples:
		vals.append(exp.smap[cid][field])
	return vals

def clusterbacteria(exp,minreads=0,uselog=True):
	"""
	cluster bacteria in an experiment according to similar behavior
	input:
	exp
	minreads - the minimal number of reads to keep before clustering (to make faster)
	uselog - True to log transform reads before normalizing, false to use full reads
	"""

	newexp=filterminreads(exp,minreads)
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

	newexp=reorderbacteria(newexp,order)
	newexp.commands.append("clusterbacteria(exp,minreads=%d)" % minreads)
	newexp.filters.append("cluster bacteria minreads=%d" % minreads)
	return newexp


def clustersamples(exp,minreads=0):
	"""
	cluster samples in an experiment according to similar behavior
	input:
	exp
	minreads - the minimal original number of reads per sample to keep it
	"""

	newexp=filterorigreads(exp,minreads)
	# normalize each row (bacteria) to sum 1
	dat=copy.copy(newexp.data)
	dat=np.transpose(dat)
	dat[dat<=2]=2
	dat=np.log2(dat)
	# cluster
	dm=spatial.distance.pdist(dat,metric='braycurtis')
	ll=cluster.hierarchy.single(dm)
	order=cluster.hierarchy.leaves_list(ll)

	newexp=reordersamples(newexp,order)
	newexp.commands.append("clustersamples(exp,minreads=%d)" % minreads)
	newexp.filters.append("cluster samples minreads=%d" % minreads)
	return newexp


def sortbacteria(exp):
	"""
	sort bacteria according to taxonomy
	"""
	tax=exp.tax
	svals,sidx=au.isort(tax)
	newexp=reorderbacteria(exp,sidx)
	newexp.filters.append('sorted bacteria by taxonomy')
	return newexp


def sortsamples(exp,field,numeric=False):
	"""
	sort samples according to field
	input:
	exp
	field - name of the field to sort by
	numeric - True for numeric values in field, false for text
	output:
	newexp - the sorted experiment
	"""
	fvals=getfieldvals(exp,field)
	if numeric:
		fvals=au.tofloat(fvals)
	svals,sidx=au.isort(fvals)
	newexp=reordersamples(exp,sidx)

	newexp.filters.append('sorted samples by field %s' % field)
	return newexp


def filterid(expdat,sids,exclude=False):
	"""
	filter bacteria keeping only ones in sids
	input:
	expdat
	sids - a list of ids
	exclude - False to keep these bacteria, True to filter away
	output:
	newexp - the filtered experiment
	"""

	if not type(sids) is list:
		sids=[sids]
	keep=[]
	au.Debug(1,'filter ids',sids)
	for cid in sids:
		for idx,tid in enumerate(expdat.sids):
			if tid==cid:
				keep.append(idx)
	if exclude:
		keep=set(range(len(expdat.sids))).difference(keep)
	keep=list(set(keep))
	au.Debug(1,'keep pos',keep)
	newexp=reorderbacteria(expdat,keep)
	if exclude:
		newexp.filters.append('Filter %d ids (exclude)' % len(sids))
	newexp.filters.append('Filter %d ids' % len(sids))
	return newexp


def filtertaxonomy(exp,tax,exact=False,inverse=False):
	"""
	filter bacteria matching a given taxonomy name
	input:
	exp
	tax - the taxonomy name to filter by
	exact - True for exact matches to tax string, false for substring
	inverse - True to throw away matching taxonomy, False to keep matching
	"""

	match=[]
	for cidx,ctax in enumerate(exp.tax):
		keep=False
		if exact:
			if ctax==tax:
				keep=True
		else:
			if tax in ctax:
				keep=True
		if inverse:
			keep=not keep
		if keep:
			match.append(cidx)
	newexp=reorderbacteria(exp,match)
	filt='filter taxonomy '
	if exact:
		filt+='exact match '
	if inverse:
		filt+='inverse '
	filt+=tax
	newexp.filters.append(filt)
	au.Debug(6,'%d bacteria left' % len(newexp.sids))
	return newexp


def calcdist(sample1,sample2,distmetric='bc'):
	"""
	calculate the distance between 2 samples using a given distance metric
	input:
	sample1,sample2 - the column arrays for the 2 samples
	distmetric - the distance metric to use:
		'bc' - bray curtis
		'bj' - binary jaccard
		'fbj' - filtered binary jaccard (using a low read threshold)
	output:
	dist - the distance
	"""
	if distmetric=='bc':
		dist=float(np.sum(np.abs(sample1-sample2)))/(np.sum(sample1)+np.sum(sample2))
	elif distmetric=='bj':
		dist=1-float(np.sum((sample1>0)*(sample2>0)))/np.sum((sample1+sample2)>0)
	elif distmetric=='fbj':
		thresh=1.9
		dist=1-float(np.sum((sample1>thresh)*(sample2>0)+(sample2>thresh)*(sample1>0)))/np.sum((sample1>thresh)+(sample2>thresh))
	else:
		au.Debug(10,'Distance meteric %s not supported' % distmetric)
		raise
	return dist


def calcdistmat(expdat,distmetric='bc'):
	"""
	calculate the distance matrix between all samples of the experiment
	input:
	expdat
	distmetric - the name of the distance metric (see calcdist for options)
	output:
	dist - the distance matrix (numpy array)
	dsamp - a dictionary with position in matrix for each sample name (key)
	"""
	msize=len(expdat.samples)
	dist=np.zeros([msize,msize])
	dsamp={}
	for cs1 in range(msize):
		for cs2 in range(msize):
			dist[cs1,cs2]=calcdist(expdat.data[:,cs1],expdat.data[:,cs2],distmetric)
	for idx,csamp in enumerate(expdat.samples):
		dsamp[csamp]=idx
	return dist,dsamp


def getgroupdist(expdat,field,distmat,dsamp,plotit=True):
	"""
	calculate the distance matrix based on groups of samples according to field
	using a distance matrix and mapping
	input:
	expdat
	field - name of the field to group by
	distmat - the distance matrix (from calcdistmat or loaddistmat)
	dsamp - the mapping of each sample id to the distance matrix position
	plotit - True to plot heatmap, False to no plot
	output:
	gdist - the group distance matrix
	uvals - a list of group names in the matrix (ordered)
	"""

	vals=getfieldvals(expdat,field)
	uvals=list(set(vals))
	gdist=np.empty([len(uvals),len(uvals)])
	gdist.fill(np.NaN)
	gmap=defaultdict(list)
	for idx,cval in enumerate(vals):
		gmap[cval].append(idx)
	for idx1,cg1 in enumerate(uvals):
		pos1=gmap[cg1]
		for idx2,cg2 in enumerate(uvals):
			pos2=gmap[cg2]
			adist=[]
			for p1 in pos1:
				for p2 in pos2:
					if p1==p2:
						continue
					adist.append(distmat[dsamp[expdat.samples[p1]],dsamp[expdat.samples[p2]]])
			gdist[idx1,idx2]=np.mean(adist)
	if plotit:
		figure()
		iax=imshow(gdist,interpolation='nearest',aspect='auto',vmin=0,vmax=1)
		ax=iax.get_axes()
		ax.set_xticks(range(len(uvals)))
		ax.set_xticklabels(uvals,rotation=90)
		ax.set_yticks(range(len(uvals)))
		ax.set_yticklabels(uvals)
		title(expdat.studyname+' '+field)
	return gdist,uvals


def joinfields(expdat,field1,field2,newfield):
	"""
	join 2 fields to create a new field for each sample
	input:
	expdat
	field1,field2 - name of the 2 fields to join
	newfield - name of new field to add
	"""

	for csamp in expdat.samples:
		expdat.smap[csamp][newfield]=expdat.smap[csamp][field1]+';'+expdat.smap[csamp][field2]
	expdat.fields.append(newfield)
	return expdat


def clearexp(expdat):
	"""
	clear experiment from missing samples and bacteria
	input:
	expdat
	output:
	newexp - the new experiment without 0 read samples or bacteria
	"""

	newexp=filterorigreads(expdat,1)
	newexp=filterminreads(expdat,1)
	return newexp


def findmislabels(expdat,field,distmetric='bc'):
	""""
	find mislabelled samples according to field
	input:
	expdat
	field - name of the field to examine (i.e. subjectid)
	distmetric - the distance meteric to use
	"""

	expdat=sortsamples(expdat,field)
	fvals=getfieldvals(expdat,field)
	ufvals=list(set(fvals))
	onames=[]
	for idx,csamp in enumerate(expdat.samples):
		onames.append(csamp+';'+fvals[idx])
	omat=np.zeros([len(fvals),len(ufvals)])
	for groupidx,groupval in enumerate(ufvals):
		cexp=filtersamples(expdat,field,groupval,exact=True)
		for aidx,aval in enumerate(expdat.samples):
			cdist=[]
			for gidx,gval in enumerate(cexp.samples):
				# don't measure distance to ourselves
				if gval==aval:
					continue
				cdist.append(calcdist(cexp.data[:,gidx],expdat.data[:,aidx],distmetric=distmetric))
			omat[aidx,groupidx]=np.mean(cdist)
	figure()
	iax=imshow(omat,interpolation='nearest',aspect='auto')
	ax=iax.get_axes()
	ax.set_xticks(range(len(ufvals)))
	ax.set_xticklabels(ufvals,rotation=90)
	ax.set_yticks(range(len(onames)))
	ax.set_yticklabels(onames)


def plotnumotus(expdat,newfig=True,threshold=0.001):
	"""
	normalized histogram of number of OTUs per sample
	input:
	expdat
	newfig - True to plot a new figure
	threshold - minimal number of reads for an otu to be present
	"""

	notus=[]
	for idx,csamp in enumerate(expdat.samples):
		cnotus=np.sum(expdat.data[:,idx]>=threshold)
		notus.append(cnotus)
	if newfig:
		figure()
	hist(notus,bins=20,range=[0,200],normed=True)
	title('num OTUs for %s' % expdat.tablefilename)


def filterseqs(expdat,seqs,exclude=False,subseq=False):
	"""
	filter sequences from the list seqs (keeping sequences appearing in the list)
	input:
	expdat
	seqs - a list of (ACGT) sequences to keep
	exclude - True to filter away instead of keep, False to keep
	subseq - if true, the sequences can be subsequence (slower). False - look only for exact match.

	output:
	newexp - the filtered experiment
	"""
	keeplist=[]
	for cseq in seqs:
		if subseq:
			for idx,coseq in enumerate(expdat.seqs):
				if len(cseq)>len(coseq):
					if coseq in cseq:
						keeplist.append(idx)
						break
				else:
					if cseq in coseq:
						keeplist.append(idx)
						break
		else:
			if cseq in expdat.seqdict:
				keeplist.append(expdat.seqdict[cseq])
			else:
				au.Debug(7,'sequence not in experiment',cseq)
	if exclude:
		keeplist=list(set(range(len(expdat.seqs))).difference(keeplist))
	newexp=reorderbacteria(expdat,keeplist)
	return newexp


def filterknownbact(expdat,cdb,exclude=False):
	"""
	filter keeping only bacteria which we know about in cooldb
	input:
	expdat
	cdb - the cooldb (fromn cooldb.loaddb)
	exclude - True to throw away known bacteria, False to keep only them
	output:
	newexp - the filtered experiment
	"""
	known=[]
	for idx,cseq in enumerate(expdat.seqs):
		if len(cooldb.getseqinfo(cdb,cseq))>0:
			known.append(idx)
	au.Debug(2,'Found %d sequences known in cooldb' % len(known))
	if exclude:
		known=set(range(len(expdat.seqs))).difference(known)
	newexp=reorderbacteria(expdat,known)
	if not exclude:
		newexp.filters.append('filter cooldb known bacteria')
	else:
		newexp.filters.append('filter exclude cooldb known bacteria')
	au.Debug(6,'%d bacteria left' % len(newexp.sids))
	return newexp



def pulsetaxplot(expdat,taxonomy):
	expdat=filtertaxonomy(expdat,taxonomy)
	expdat=filterminreads(expdat,30)
	for cid in expdat.sids:
		pulseplot(expdat,cid)

def pulseplot(expdat,cid):
	plen=[14,10,6,3,1]
	figure()
	expdat=filterid(expdat,cid)
	title('%s' % (expdat.tax[0]))
	for day in ['5','1','3','4','2']:
		cexp=filtersamples(expdat,'Day',day)
		cexp=sortsamples(cexp,'Elapsed Time',numeric=True)
		times=au.tofloat(getfieldvals(cexp,'Elapsed Time'))
		print times
		print cexp.sids
		plot(times,cexp.data[0,:])
	legend(plen)


def addplotmetadata(expdat,field,value=False,inverse=False,color='g',ax=False,beforesample=True,partial=False):
	"""
	plot lines on an experiment plot from plotexp. NOTE: need to use with the output of plotexp!!!
	input:
	expdat
	field - name of the metadata field to use
	value - the value that when present we will plot, or False to plot whenever not empty
	inverse - inverse the logic of value
	color - the color to plot
	beforesample - True if it happened between the prev. and current sample, False if in the middle of the sample
	partial - True to allow substring match, False for exact match
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
			au.Debug(1,'Plot line %d',idx)
			plt.plot([idx+offset,idx+offset],[-0.5,len(expdat.sids)-0.5],color)



def norm2(expdat,numremove):
	mv=np.mean(expdat.data,axis=1)
	sv,si=au.isort(mv)
	keeppos=(si[:-numremove-1])
	nval=np.sum(expdat.data[keeppos,:],axis=0)
	newexp=copy.deepcopy(expdat)
	for cs in range(len(newexp.samples)):
		newexp.data[:,cs]=newexp.data[:,cs]*10000/nval[cs]
	return newexp


def normalizeprctile(expdat,percent=80):
	"""
	normalize reads per experiment so percentile (rather than mean) will be normalized
	used to reduce effect of outliers
	note normalization is done on the same set of bacteria for all samples
	input:
	expdat
	percent - the percentile to normalize (0-100)

	output:
	newexp - the new normalized experiment
	"""

	# select the bacteria to use
	newexp=filterminreads(expdat,1*len(expdat.samples))

	percvals=np.percentile(newexp.data,percent,axis=0)
	plt.figure()
	plt.plot(percvals)
	percvals=percvals/np.mean(percvals)
	newexp=copy.deepcopy(expdat)
	for idx,samp in enumerate(expdat.samples):
		newexp.data[:,idx]=newexp.data[:,idx]*percvals[idx]
	return newexp


def sortbyfreq(expdat,field=False,value=False,exact=False):
	"""
	sort bacteria in experiment according to frequency
	sorting is performed based on a subset of samples (field/val/exact) and then
	all the experiment is sorted according to them
	input:
	expdat
	field - name of the field to filter samples for freq. sorting or False for all samples
	value - value of samples to use for the freq. sorting
	exact - is the value exact or partial string

	output:
	newexp - the experiment with bacteria sorted according to subgroup freq.
	"""

	if field:
		texp=filtersamples(expdat,field,value,exact=exact)
	else:
		texp=copy.deepcopy(expdat)

	meanvals=np.mean(texp.data,axis=1)
	svals,sidx=au.isort(meanvals)

	newexp=reorderbacteria(expdat,sidx)
	return newexp


def savecoolseqs(expdat,cdb,seqs,description):
	"""
	save sequences to the cooldb database
	input:
	expdat
	seqs - a list of sequences (ACGT) to save
	description - the name to put in coolseqs database
	"""

	numwritten=0
	for cseq in seqs:
		if cseq in expdat.seqdict:
			seqpos=expdat.seqdict[cseq]
			cooldb.saveseq(db=cdb,seq=cseq,ggid=expdat.sids[seqpos],taxonomy=expdat.tax[seqpos],filename=expdat.tablefilename,expdescription=expdat.studyname,description=description)
			numwritten+=1
		else:
			au.Debug(9,'sequence:\n%s\nnot in experiment!!!' % cseq)
	au.Debug(6,'Wrote %d sequences to cooldb' % numwritten)



def plotseqfreq(expdat,seqs,toaxis=False,xfield=False,normalizey=False):
	"""
	plot the frequency of sequences in seq as a function of the sortfield
	input:
	expdat
	seqs - a list of sequnces (acgt) to plot
	toaxis - if not empty - the axis to plot to, False plot a new figure
	xfield - if not False - space the points on the x axis according to (numeric) value of in xfield, False - just according to the sorted order
	normalizey - True: normalize all y values to 0-1, False: no normalization
	"""

	if xfield:
		xdat=au.tofloat(getfieldvals(expdat,xfield))
	else:
		xdat=range(len(expdat.samples))
	sv,si=au.isort(xdat)
	if not toaxis:
		figure()
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
	labels=au.clipstrings(labels,20,reverse=True)
	toaxis.legend(labels,prop={'size':6})
	toaxis.set_xticks(xdat)
#	toaxis.set_xticklabels(labs,rotation=45,ha='right')


def filtersamplesfromfile(expdat,filename,field='#SampleID',exclude=False):
	"""
	filter samples based on a list in a file (one line per sample)
	input:
	expdat
	filename - name of the text file (1 line per sample)
	field - field that the file contains
	exclude - true to throw away instead of keep
	output:
	newexp - the filtered experiment
	"""

	fl=open(filename,'r')
	vals=[]
	for cline in fl:
		vals.append(cline)
	fl.close()

	keep=[]
	for cidx,csamp in enumerate(expdat.samples):
		keepit=False
		if expdat.smap[csamp][field] in vals:
			keepit=True
		# if exclude reverse the decision
		if exclude:
			keepit=not keepit
		if keepit:
			keep.append(cidx)
	newexp=reordersamples(expdat,keep)

	fstr="filter data from file %s in %s " % (filename,field)
	if exclude:
		fstr+=" (exclude)"
	newexp.filters.append(fstr)
	au.Debug(6,'%d Samples left' % len(newexp.samples))
	return newexp



def joinexperiments(exp1,exp2,missingval='NA',origfieldname='origexp'):
	"""
	join 2 experiments into a new experiment. adding a new field origfieldname
	input:
	exp1,exp2 - the experiments to join
	missingval - string to put when field not in mapping file of one of the experiments
	origfieldname - name of the new field to add which contains the original experiment name
	"""

	# join the sequences of both experiments
	# ASSUMING SAME SEQ LENGTH!!!!
	allseqs=list(set(exp1.seqs) | set(exp2.seqs))
	alldict={}
	alltax=[]
	allids=[]
	for idx,cseq in enumerate(allseqs):
		alldict[cseq]=idx
	# make the new joined data for each experiment
	dat1=np.zeros((len(allseqs),np.size(exp1.data,1)))
	for idx,cseq in enumerate(allseqs):
		if cseq in exp1.seqdict:
			dat1[idx,:]=exp1.data[exp1.seqdict[cseq],:]
			alltax.append(exp1.tax[exp1.seqdict[cseq]])
			allids.append(exp1.sids[exp1.seqdict[cseq]])
		else:
			alltax.append(exp2.tax[exp2.seqdict[cseq]])
			allids.append(exp2.sids[exp2.seqdict[cseq]])

	dat2=np.zeros((len(allseqs),np.size(exp2.data,1)))
	for idx,cseq in enumerate(allseqs):
		if cseq in exp2.seqdict:
			dat2[idx,:]=exp2.data[exp2.seqdict[cseq],:]

	newexp=copy.deepcopy(exp1)
	# concatenate the reads
	newexp.data=np.concatenate((dat1,dat2), axis=1)
	newexp.seqdict=alldict
	newexp.seqs=allseqs
	newexp.tax=alltax
	newexp.sids=allids
	newexp.sids=newexp.seqs
	newexp.samples = list(exp1.samples) + list(exp2.samples)
	newexp.origreads=exp1.origreads+exp2.origreads
	newexp.fields=list(set(exp1.fields+exp2.fields))

	for cfield in newexp.fields:
		if cfield in exp1.fields:
			continue
		for csamp in exp1.samples:
			newexp.smap[csamp][cfield]=missingval

	for csamp in exp2.samples:
		newexp.smap[csamp]={}
		for cfield in newexp.fields:
			if cfield in exp2.fields:
				newexp.smap[csamp][cfield]=exp2.smap[csamp][cfield]
			else:
				newexp.smap[csamp][cfield]=missingval

	for csamp in exp1.samples:
		newexp.smap[csamp][origfieldname]=exp1.studyname
	for csamp in exp2.samples:
		newexp.smap[csamp][origfieldname]=exp2.studyname

	newexp.filters.append('joined with %s' % exp2.studyname)
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

	exp1=filtersamples(expdat,field,val1,exact=True)
	exp2=filtersamples(expdat,field,val2,exact=True)
	m1=np.mean(np.log2(exp1.data+2),axis=1)
	m2=np.mean(np.log2(exp2.data+2),axis=1)
	diff=(m1-m2)/(m1+m2+20)
	sv,si=au.isort(diff)
	newexp=reorderbacteria(expdat,si)
	print(sv[-10:])
	return newexp


def getdiffsig(expdat,field,val1,val2=False,method='mean',numperm=1000,maxfval=0.1):
	"""
	test the difference between 2 groups (val1 and val2 in field field)
	for bacteria that have a high difference.
	input:
	expdat
	field - the field for the 2 categories
	val1 - values for the first group
	val2 - value for the second group or false to compare to all other
	method - the test to compare the 2 groups:
		mean - absolute difference in mean frequency
		bin - abs diff in binary presence/absence
		rank - abs diff in rank order (to ignore outliers)
		freqpres - abs diff in frequency
	numperm - number of random permutations to run
	maxfval - the maximal f-value (FDR) for a bacteria to keep

	output:
	newexp - the experiment with only significant (FDR<=maxfval) difference, sorted according to difference
	"""

	minthresh=2
	exp1=filtersamples(expdat,field,val1,exact=True)
	if val2:
		exp2=filtersamples(expdat,field,val2,exact=True)
	else:
		exp2=filtersamples(expdat,field,val1,exact=True,exclude=True)
	cexp=joinexperiments(exp1,exp2)
	len1=len(exp1.samples)
	len2=len(exp2.samples)
	dat=cexp.data
	dat[dat<minthresh]=minthresh
	dat=np.log2(dat)
	numseqs=len(cexp.seqs)

	eps=0.000001
	if method=='mean':
		eps=20
	if method=='rank':
		for idx in range(numseqs):
			dat[idx,:]=stats.rankdata(dat[idx,:])
	if method=='bin':
		dat=(dat>np.log2(minthresh))


	m1=np.mean(dat[:,0:len1],axis=1)
	m2=np.mean(dat[:,len1:],axis=1)
	odif=(m1-m2)/(m1+m2+eps)
	print(np.max(odif))
	print(np.min(odif))
	alldif=np.zeros([len(cexp.sids),numperm])
	for x in range(numperm):
		rp=np.random.permutation(len1+len2)
		m1=np.mean(dat[:,rp[0:len1]],axis=1)
		m2=np.mean(dat[:,rp[len1:]],axis=1)
		diff=(m1-m2)/(m1+m2+eps)
		alldif[:,x]=diff
	pval=[]
	for crow in range(len(odif)):
		cpval=float(np.sum(np.abs(alldif[crow,:])>=np.abs(odif[crow])))/numperm
		# need to remember we only know as much as the number of permutations - so add 1 as upper bound for fdr
		cpval=min(cpval+(1.0/numperm),1)
		pval.append(cpval)
	# NOTE: maybe use better fdr (this one not tested...)
	fval=au.fdr(pval)
	keep=np.where(np.array(fval)<=maxfval)
	seqlist=[]
	for cidx in keep[0]:
		seqlist.append(cexp.seqs[cidx])
	newexp=filterseqs(expdat,seqlist)
	odif=odif[keep[0]]
	sv,si=au.isort(odif)
	newexp=reorderbacteria(newexp,si)
	newexp.filters.append('differential expression (%s) in %s between %s and %s' % (method,field,val1,val2))
	return newexp


def saveseqsfasta(expdat,seqs,filename):
	"""
	save bacteria from list seqs to a fasta file
	input:
	expdat - (for taxonomy) or False
	seqs - the sequences to save
	filename - the name of the output fasta file
	"""

	fl=open(filename,'w')
	for idx,cseq in enumerate(seqs):
		if expdat:
			if cseq in expdat.seqdict:
				ctax=expdat.tax[expdat.seqdict[cseq]]
				cid=expdat.sids[expdat.seqdict[cseq]]
			else:
				ctax='NA'
				cid='NA'
		fl.write('>%d-%s-%s\n' % (idx,str(cid),str(ctax)))
		fl.write(cseq+'\n')
	fl.close()


def saveexpseqs(expdat,filename):
	"""
	save all bacteria from an experiment as a fasta file
	input:
	expdat
	filename - the name of the output fasta file
	"""
	saveseqsfasta(expdat,expdat.seqs,filename)


def forgeorg(expdat,filename,num):
	fl=open(filename,'w')
	fl.write('#SampleID\tfreq%d\n' % num)
	for idx,csamp in enumerate(expdat.samples):
		fl.write(csamp)
		fl.write('\t')
		fl.write(str(np.sum(expdat.data[:,idx])))
		fl.write('\n')
	fl.close()



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
	figure()
	plot(x,y)
	plot(x,yb)
	plot(x,ynb)
	title(expdat.studyname+' (%d samples)' % len(expdat.samples))
	xlabel('number of reads')
	ylabel('fraction of samples with >= reads')
	legend(['all','blanks','non blanks'])
	show()


def filterfasta(expdat,filename,exclude=False,subseq=False):
	"""
	Filter bacteria based on a fasta file (keeping only sequences in the fasta file)
	input:
	expdat
	filename - the fasta file name
	exclude - False to keep only seqs from the file, True to throw away seqs from the file
	subseq - match subsequences (slower)

	output:
	newexp - the filtered experiment
	"""

	seqs=au.readfastaseqs(filename)
	newexp=filterseqs(expdat,seqs,exclude=exclude,subseq=subseq)
	filt='Filter sequences from file '+filename
	if exclude:
		filt+=' (Exclude)'
	if subseq:
		filt+=' (subseq)'
	newexp.filters.append(filt)
	return newexp


def filterbacteriafromfile(expdat,filename,exclude=False,subseq=False):
	"""
	filter bacteria from an experiment based on a file with sequences (one per line)
	input:
	expdat
	filename - name of the sequence file (1 per line)
	exclude - remove bacteria from the file instead of keeping them
	subseq - the sequences in the file can be subsequences of the experiment sequences (different lengths). but slower.

	output:
	newexp - the filtered experiment
	"""
	fl=open(filename,'rU')
	seqs=[]
	for cline in fl:
		seqs.append(cline.strip())
	newexp=filterseqs(expdat,seqs,exclude=exclude,subseq=False)
	filt='Filter sequences from file '+filename
	if exclude:
		filt+=' (Exclude)'
	if subseq:
		filt+=' (subseq)'
	newexp.filters.append(filt)
	return newexp


def bicluster(expdat,numiter=5,startb=False,starts=False,method='zscore'):
	"""
	cluster bacteria and samples from subgroup
	input:
	expdat
	numiter - number of iterations to run the biclustering
	startb - start list of bacteria [acgt] of False for random
	method - 'zscore' or 'ranksum'. the method for choosing which sample/bacteria to keep
	"""

	dat=copy.copy(expdat.data)
#	dat[dat<20]=20
#	dat=np.log2(dat)
	dat=(dat>1)
	print('li')
	bdat=dat
#	bdat=scale(dat,axis=1,copy=True)
	nbact=np.size(dat,0)
	nsamp=np.size(dat,1)
	allsamp=np.arange(nsamp)
	allbact=np.arange(nbact)

	bthresh=0.5
	sthresh=0
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

	x=np.setdiff1d(allsamp,usamp)
	sampo=np.concatenate((usamp,x))
	bacto=np.concatenate((ubact,np.setdiff1d(allbact,ubact)))

	newexp=reordersamples(expdat,sampo)
	newexp=reorderbacteria(newexp,bacto)
	newexp.filters.append('biclustering')
	return newexp



def savebiom(expdat,filename):
	"""
	save experiment to text biom table and mapping file
	"""

	mf=open(filename+'.map.txt','w')
	mf.write('#SampleID')
	for cfield in expdat.fields:
		if cfield=='#SampleID':
			continue
		mf.write('\t%s' % cfield)
	mf.write('\n')
	for csamp in expdat.samples:
		mf.write('%s' % csamp)
		for cfield in expdat.fields:
			if cfield=='#SampleID':
				continue
			mf.write('\t')
			mf.write(expdat.smap[csamp][cfield])
		mf.write('\n')
	mf.close()

	tf=open(filename+'.table.txt','w')
	tf.write('# Saved biom table from python analysis\n')
	tf.write('#OTUID')
	for csamp in expdat.samples:
		tf.write('\t%s' % csamp)
	tf.write('\tTaxonomy\n')
	for idxseq,cseq in enumerate(expdat.seqs):
		tf.write('%s' % cseq)
		for idx,csamp in enumerate(expdat.samples):
			tf.write('\t%d' % expdat.data[idxseq,idx])
		tf.write('\t%s\n' % expdat.tax[idxseq])
	tf.close()



def getexpdbsources(expdat):
	if not expdat.seqdb:
		au.Debug(9,'No sequence database loaded')
		return
	dat=bactdb.GetDBSource(expdat.seqdb,expdat.seqs)

	newexp=copy.deepcopy(expdat)

	THRESH=0.001
	used=np.arange(np.size(dat,0))
	done=False
	while not done:
		npsamp=np.sum(dat[used,:]>=THRESH,axis=0)
		pos=np.argmax(npsamp)
		print('sample is %d, size is %d' % (pos,npsamp[pos]))
		sid,sname=bactdb.GetSampleStudy(expdat.seqdb,pos+1)
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
	return newexp



def clipseqs(expdat,startpos):
	"""
	clip the first nucleotides in all sequences in experiment
	input:
	expdat
	startpos - the position to start from (0 indexed)
	output:
	newexp - new experiment with all sequences clipped and joined identical sequences
	"""

	newexp=copy.deepcopy(expdat)
	newseqs=[]
	newdict={}
	keeppos=[]
	for idx,cseq in enumerate(newexp.seqs):
		cseq=cseq[startpos:]
		if cseq in newdict:
			newexp.data[newdict[cseq],:] += newexp.data[idx,:]
		else:
			newdict[cseq]=idx
			newseqs.append(cseq)
			keeppos.append(idx)
	newexp=reorderbacteria(newexp,keeppos)
	newexp.seqs=newseqs
	newexp.seqdict=newdict
	newexp.filters.append("trim %d nucleotides" % startpos)
	return newexp


def normalizereads(expdat,fixorig=True,numreads=10000):
	"""
	normalize the number of reads per sample to 10k
	input:
	expdat
	fixorig - True to fix origreads with the same ratio, False to keep as before
	numreads - the number of reads to normalize to

	output:
	newexp - the normalized experiment
	"""
	newexp=copy.deepcopy(expdat)
	for idx,csamp in enumerate(newexp.samples):
		totreads=np.sum(newexp.data[:,idx])
		if totreads==0:
			continue
		ratio=float(numreads)/totreads
		newexp.data[:,idx]=newexp.data[:,idx]*ratio
		if fixorig:
			au.Debug(6,'fixing original frequencies')
			newexp.origreads[idx]=float(newexp.origreads[idx])/ratio
	newexp.filters.append("renormalized reads to sum %d" % numreads)
	return newexp



def filterandnormalize(expdat,seqs,exclude=False,subseq=False,numreads=10000):
	"""
	filter away sequences in seqs
	and then renormalize the data and recalculate origreads per sample
	input:
	expdat
	seqs - a list of sequences (ACGT) to remove
	exclude - False to remove seqs, True to only keep seqs
	subseq - False to look for exact match only, True to look for subsequence match (slower)

	output:
	newexp - the experiment without seqs and renormalied to 10k reads/sample
	"""

	newexp=filterseqs(expdat,seqs,exclude=not(exclude),subseq=subseq)
	newexp=normalizereads(newexp,fixorig=True,numreads=numreads)
	return newexp


def filterannotations(expdat,annotation,cdb,exclude=False):
	"""
	filter keeping only samples which have annotation in their cooldb description
	input:
	expdat
	annotation - substring of the annotation (case insensitive)
	cdb - the database of cool sequences (from cooldb.load())
	exclude - False to keep matching bacteria, True to remove matching bacteria

	output:
	newexp - the filtered experiment
	"""

	keeplist=[]
	for idx,cseq in enumerate(expdat.seqs):
		keep=False
		info=cooldb.getseqinfo(cdb,cseq)
		for cinfo in info:
			if annotation.lower() in str(cinfo).lower():
				keep=True
		if exclude:
			keep = not keep
		if keep:
			keeplist.append(idx)
	newexp=reorderbacteria(expdat,keeplist)
	newexp.filters.append('Filter annotations %s' % annotation)
	au.Debug(6,'%d bacteria found' % len(keeplist))
	return newexp


def savetsvtable(expdat,filename,logtransform=True):
	"""
	save an experiment as a tab separated table, with columns for samples and rows for bacteria
	for jose navas long babies paper
	input:
	expdat
	filename - name of the output tsv file
	minreads - save only bacteria with >=minreads reads
	logtransform - True to save the log2 of the reads, False to save the reads
	"""

	ldat=copy.copy(expdat.data)
	if logtransform:
		ldat[np.where(ldat<1)]=1
		ldat=np.log2(ldat)

	of=open(filename,'w')
	of.write("Taxonomy\tSequence")
	for csamp in expdat.samples:
		of.write("\t%s" % csamp)
	of.write("\n")
	for idx,cseq in enumerate(expdat.seqs):
		of.write("%s\t%s" % (expdat.tax[idx],cseq))
		for cval in ldat[idx,:]:
			of.write("\t%f" % cval)
		of.write("\n")
	of.close()


def getnucdistribution(expdat,position):
	"""
	get the distribution of nucleotides in the positions in position
	note this is unweighted
	input:
	expdat
	position - a list of positions (0 based) to test

	output:
	"""

	retv=np.zeros((len(position),6))
	for cseq in expdat.seqs:
		cseqn=au.SeqToArray(cseq)
		for idx,cpos in enumerate(position):
			retv[idx,cseqn[cpos]]+=1
	figure()
	for cnuc in range(np.size(retv,axis=1)):
		for crow in range(np.size(retv,axis=0)-1):
			bar(np.arange(np.size(retv,axis=0)),retv[crow+1,:],bottom=retv[crow,:])
	return (retv)
