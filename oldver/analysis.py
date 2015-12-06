#!/usr/bin/env python

# amnonscript

"""
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

__version__ = "0.9"

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
import sklearn.metrics
import sklearn.cross_validation
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

		# the tree structure of the sequences (from loadexptree)
		self.tree=False

		# the experiment type ('biom' or 'meta' for metabolite)
		self.datatype=''


def copyexp(expdat):
	"""
	copy an experiment (duplicating the important fields)
	input:
	expdat
	output:
	newexp
	"""

	newexp=copy.copy(expdat)
	newexp.data=copy.deepcopy(expdat.data)
	newexp.samples=copy.deepcopy(expdat.samples)
	newexp.seqs=copy.deepcopy(expdat.seqs)
	newexp.sids=copy.deepcopy(expdat.sids)
	newexp.seqdict=copy.deepcopy(expdat.seqdict)
	newexp.tax=copy.deepcopy(expdat.tax)
	newexp.plotmetadata=copy.deepcopy(expdat.plotmetadata)
	newexp.smap=copy.deepcopy(expdat.smap)
	newexp.fields=copy.deepcopy(expdat.fields)
	newexp.filters=copy.deepcopy(expdat.filters)
	newexp.origreads=copy.deepcopy(expdat.origreads)
	newexp.origotunames=copy.deepcopy(expdat.origotunames)

	return newexp


def load(tablename, mapname='map.txt', taxfile='', nameisseq=True,addsname='',studyname=False,tabletype='biom',normalize=True):
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
	tabletype:
		'biom' - a biom table
		'meta' - a metabolomics table (row per sample, col per metabolite, can contain duplicate metaboliteids)
	normalize - True to normalize to 10k reads per sample, False to not normalize (change to mean 10k reads/sample)
	output:
	an experiment class for the current experiment
	"""

	# au.Debug(6,'Loading mapping file')
	# # load the mapping file
	# mapf = open(mapname, 'rU')
	# reader = csv.DictReader(mapf, delimiter='\t')
	# fields = reader.fieldnames
	# smap = {}
	# mapsamples = []
	# for cline in reader:
	# 	cid = cline['#SampleID']
	# 	smap[cid] = cline
	# 	mapsamples.append(cid)
	# mapf.close()
	# au.Debug(6,'number of samples in map is %d' % len(mapsamples))

	if tabletype=='biom':
		au.Debug(6,'Loading biom table')
		# load the biom table
		table = biom.load_table(tablename)
	elif tabletype=='meta':
		au.Debug(6,'Loading metabolite table')
		# load the metabolite table and turn it into a biom table
		fl=open(tablename,'rU')
		head=fl.readline().rstrip('\n')
		# its a csv
		headsplit=head.split(',')
		headsplit=headsplit[1:]
		# look if we have strange entries (like 'x','y','z')
		usepos=[]
		metabolites=[]
		for idx,cmet in enumerate(headsplit):
			if cmet[0].isdigit():
				usepos.append(idx+1)
				metabolites.append("%s-%d" % (cmet,idx))
			else:
				au.Debug(7,'Metabolite %s (col %d) not numeric!' % (cmet,idx))
		# load sample names
		sampnames=[]
		for cline in fl:
			cline=cline.rstrip('\n')
			cdat=cline.split(',')
			sampnames.append(cdat[0])
		fl.close()
		# now load the table data:
		dat=np.genfromtxt(tablename,skip_header=1,usecols=usepos,delimiter=',')
		dat=np.transpose(dat)
		# and create the biom table:
		table=biom.table.Table(dat,metabolites,sampnames)
		# and add metabolite name as taxonomy:
		taxdict={}
		for cmet in metabolites:
			taxdict[cmet]={'taxonomy': cmet}
		table.add_metadata(taxdict,axis='observation')
	else:
		au.Debug(9,'Table type %s not supported' % tabletype)
		return False

	if addsname!='':
		idtable={}
		ids=table.ids(axis='sample')
		for cid in ids:
			idtable[cid]=addsname+cid
		table=table.update_ids(idtable,axis='sample')

	smap = {}
	mapsamples = []
	if mapname:
		au.Debug(6,'Loading mapping file')
		# load the mapping file
		mapf = open(mapname, 'rU')
		reader = csv.DictReader(mapf, delimiter='\t')
		fields = reader.fieldnames
		for cline in reader:
			cid = cline['#SampleID']
			smap[cid] = cline
			mapsamples.append(cid)
		mapf.close()
		au.Debug(6,'number of samples in map is %d' % len(mapsamples))
	else:
		au.Debug(6,'No mapping file supplied - using just sample names')
		tablesamples = table.ids(axis='sample')
		for cid in tablesamples:
			smap[cid]={'#SampleID':cid}
			mapsamples.append(cid)
		fields=['#SampleID']

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
	exp.datatype=tabletype
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
	if tabletype=='meta':
		normalize=False

	if normalize:
		exp.data=10000*exp.data/colsum
	else:
		exp.data=10000*exp.data/np.mean(colsum)


	exp.commands.append("load('"+tablename+"',mapname='"+mapname+")")
	exp.filters.append('loaded table=%s, map=%s' % (tablename,mapname))
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
	if event.key==',':
	# select next bacteria
		cax.guiwin.clearselection()
		cax.lastselect+=1
		cax.guiwin.selectbact([cax.lastselect])
		cax.guiwin.updateinfo(cax.guiwin.csamp,cax.lastselect)
		if cexp.cdb:
			info = cooldb.getseqinfo(cexp.cdb,cexp.seqs[cax.lastselect])
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
			info = cooldb.getseqinfo(cexp.cdb,cexp.seqs[cax.lastselect])
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
		labs=au.clipstrings(cexp.tax,25,reverse=True)
		if cexp.cdb:
			for idx,cseq in enumerate(cexp.seqs):
				info=cooldb.getseqinfo(cexp.cdb,cseq)
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


def filterorigreads(exp,minreads,inplace=False):
	"""
	filter away all samples that contained originally less than minreads
	input:
	exp - the experiment
	minreads - the minimum number of reads of the sample in the biom table to filter if less
	inplace - True to replace current experiment, False to create a new one
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


def filtersamples(exp,field,filtval,exact=True,exclude=False,numexpression=False):
	"""
	filter samples in experiment according to field
	input:
	exp
	field - name of the field
	filtval - the string to filter (if a list of strings, filter if any in the list)
	exact - True for exact match, False for substring
	exclude - False to keep only matching samples, True to exclude matching samples
	numexpression - True if val is a python expression, False if just a value. For an expression assume value is the beggining of the line (i.e. '<=5')
	"""

	if not isinstance(filtval,list):
		filtval=[filtval]

	keep=[]
	for cidx,csamp in enumerate(exp.samples):
		keepit=False
		for filt in filtval:
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


def filtertaxonomy(exp,tax,exact=False,exclude=False):
	"""
	filter bacteria matching a given taxonomy name
	input:
	exp
	tax - the taxonomy name to filter by
	exact - True for exact matches to tax string, false for substring
	exclude - True to throw away matching taxonomy, False to keep matching
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
		if exclude:
			keep=not keep
		if keep:
			match.append(cidx)
	newexp=reorderbacteria(exp,match)
	filt='filter taxonomy '
	if exact:
		filt+='exact match '
	if exclude:
		filt+='exclude '
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



def loaddistmat(expdat,dmfilename):
	"""
	load a distance matrix (from qiime) for analysis
	input:
	expdat
	dmfilename - name of the distance matrix file

	output:
	distmat - the distance matrix
	dsamp - the mapping to position in the mapping file for each distmat entry
	"""

	fl=open(dmfilename,'rU')
	# get the column ids
	head=fl.readline().strip('\n')
	ids=head.split('\t')
	ids=ids[1:]
	dist=np.array([])
	snames={}
	for idx,cline in enumerate(fl):
		cline=cline.strip("\n")
		vals=cline.split('\t')
		dist=np.vstack((dist,au.tofloat(vals[1:]))) if dist.size else np.array(au.tofloat(vals[1:]))
		if not vals[0]==ids[idx]:
			au.Debug(9,"strange! line %d row head %s but col head %s" % (idx,vals[0],ids[idx]))
		snames[vals[0]]=idx
	fl.close()
	expkeep=[]
	distorder=[]
	dsamp={}
	for idx,csamp in enumerate(expdat.samples):
		if csamp in snames:
			distorder.append(snames[csamp])
			expkeep.append(idx)
			dsamp[csamp]=snames[csamp]
	au.Debug(6,"%d samples in dist mat, %d samples in experiment" % (len(ids),len(expdat.samples)))
	au.Debug(6,"%d samples to keep from dist mat, %d samples to keep from experiment" % (len(distorder),len(expkeep)))
#	dist=dist[distorder,:]
#	dist=dist[:,distorder]
	return dist,dsamp


def getgroupdist(expdat,field,distmat,dsamp,plotit=True,plottype='heatmap',uvals=False):
	"""
	calculate the distance matrix based on groups of samples according to field
	using a distance matrix and mapping
	input:
	expdat
	field - name of the field to group by
	distmat - the distance matrix (from calcdistmat or loaddistmat)
	dsamp - the mapping of each sample id to the distance matrix position
	plotit - True to plot heatmap, False to no plot
	plottype:
		'heatmap'
		'hist'
	uvals - false to plot all values, or a list of values to plot only them (in field)
	output:
	gdist - the group distance matrix
	uvals - a list of group names in the matrix (ordered)
	"""

	vals=getfieldvals(expdat,field)
	if not uvals:
		uvals=list(set(vals))
	gdist=np.empty([len(uvals),len(uvals)])
	gdist.fill(np.NaN)
	gmap=defaultdict(list)
	distdict={}
	for idx,cval in enumerate(vals):
		gmap[cval].append(idx)
	for idx1,cg1 in enumerate(uvals):
		pos1=gmap[cg1]
		for idx2,cg2 in enumerate(uvals):
			pos2=gmap[cg2]
			adist=[]
			for p1 in pos1:
				if expdat.samples[p1] not in dsamp:
					continue
				for p2 in pos2:
					if expdat.samples[p2] not in dsamp:
						continue
					if p1==p2:
						continue
					adist.append(distmat[dsamp[expdat.samples[p1]],dsamp[expdat.samples[p2]]])
			distdict[(cg1,cg2)]=adist
			gdist[idx1,idx2]=np.mean(adist)
	if plotit:
		figure()
		if plottype=='heatmap':
			iax=imshow(gdist,interpolation='nearest',aspect='auto',vmin=0,vmax=1)
			ax=iax.get_axes()
			ax.set_xticks(range(len(uvals)))
			ax.set_xticklabels(uvals,rotation=90)
			ax.set_yticks(range(len(uvals)))
			ax.set_yticklabels(uvals)
			title(expdat.studyname+' '+field)
		elif plottype=='hist':
			pl=[]
			pairs=[]
			names=[]
			for k,v in distdict.items():
				ks=set(k)
				if ks in pairs:
					continue
				pl.append(v)
				pairs.append(ks)
				names.append(k)
			hist(pl,alpha=0.5,normed=True,bins=50,range=[0,1])
			legend(names)
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
	normalized histogram of number of OTUs per sample (alpha diversity)
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
	newexp.filters.append("sort by freq field=%s value=%s" % (field,value))
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
	newexp.fields.append(origfieldname)

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


def getdiffsigall(expdat,field,val1,val2=False,numperm=1000,maxfval=0.1):
	"""
	Get the differentially abundant bacteria (using getdiffsig) using all methods possible.
	Sort results according to combined effect order
	input:
	see getdiffsig()

	output:
	newexp - the new experiment with bacteria significantly differentiating the 2 groups by at least 1 method
	"""

	methods=['mean','binary','ranksum','freqpres']
	res=[]
	for cmethod in methods:
		res.append(getdiffsig(expdat,field=field,val1=val1,val2=val2,method=cmethod,numperm=numperm,maxfval=maxfval))

	keep=[]
	keeporder=[]
	for cidx,cseq in enumerate(expdat.seqs):
		pos=[]
		for cres in res:
			if cseq in cres.seqdict:
				pos.append(float(cres.seqdict[cseq])/len(cres.seqs))
		if len(pos)>0:
			keep.append(cidx)
			keeporder.append(np.mean(pos))
	keep=np.array(keep)
	if len(keep)>0:
		si=np.argsort(keeporder)
		newexp=reorderbacteria(expdat,keep[si])
		newexp.filters.append('differential expression (all) in %s between %s and %s' % (field,val1,val2))
		return newexp
	else:
		au.Debug(6,'No bacteria found')
		return False


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
		pass
#		eps=20
	elif method=='ranksum':
		for idx in range(numseqs):
			dat[idx,:]=stats.rankdata(dat[idx,:])
	elif method=='binary':
		dat=(dat>np.log2(minthresh))
	elif method=='freqpres':
		dat[dat<=minthresh]=np.nan
	else:
		au.Debug(9,"Method not supported!",method)
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
#		cpval=float(np.sum(np.abs(alldif[crow,:])>=np.abs(odif[crow])))/numperm
		# need to remember we only know as much as the number of permutations - so add 1 as upper bound for fdr
		cpval=min(cpval+(1.0/cnumperm),1)
#		cpval=min(cpval+(1.0/numperm),1)
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

	seqs,headers=au.readfastaseqs(filename)
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


def bicluster(expdat,numiter=5,startb=False,starts=False,method='zscore',sampkeep=0.5,bactkeep=0.25,justcount=False,numruns=1):
	"""
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
			au.Debug(6,"bactkeep %f" % bactkeep)
		if sampkeep==0:
			sampkeep=np.random.uniform(0.25,0.75)
			au.Debug(6,"sampkeep %f" % sampkeep)
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
				au.Debug(9,"biclustering method %s not supported")
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
			newexp=reordersamples(expdat,sampo)
			newexp=reorderbacteria(newexp,bacto,inplace=True)
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
	figure()
	cc,seqs,samps=bicluster(expdat,method='binary',justcount=True,numruns=numruns,numiter=numiter)
	for idx,cseqs in enumerate(seqs):
		plot(len(samps[idx]),len(cseqs),'xr')
	rp=randomizeexp(expdat,normalize=False)
	cc,seqs,samps=bicluster(rp,method='binary',justcount=True,numruns=numruns,numiter=numiter)
	for idx,cseqs in enumerate(seqs):
		plot(len(samps[idx]),len(cseqs),'xk')
	rp=randomizeexp(expdat,normalize=True)
	cc,seqs,samps=bicluster(rp,method='binary',justcount=True,numruns=numruns,numiter=numiter)
	for idx,cseqs in enumerate(seqs):
		plot(len(samps[idx]),len(cseqs),'xk')
	xlabel('# Samples in cluster')
	ylabel('# Bacteria in cluster')


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



def getexpdbsources(expdat,seqdb=False):
	if not seqdb:
		if not expdat.seqdb:
			au.Debug(9,'No sequence database loaded')
			return
		else:
			seqdb=expdat.seqdb

	dat=bactdb.GetDBSource(seqdb,expdat.seqs)

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


def clipseqs(expdat,startpos,addseq='TAC'):
	"""
	clip the first nucleotides in all sequences in experiment
	to fix offset in sequencing
	input:
	expdat
	startpos - the position to start from (0 indexed) or negative to add nucleotides
	addseq - the sequence to add (just a guess) if startpos is negative
	output:
	newexp - new experiment with all sequences clipped and joined identical sequences
	"""

	newexp=copy.deepcopy(expdat)
	newseqs=[]
	newdict={}
	keeppos=[]
	for idx,cseq in enumerate(newexp.seqs):
		if startpos>=0:
			cseq=cseq[startpos:]
		else:
			cseq=addseq[:abs(startpos)]+cseq
			cseq=cseq[:len(expdat.seqs[0])]
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

	if inplace:
		newexp=expdat
	else:
		newexp=copy.deepcopy(expdat)

	for idx,csamp in enumerate(newexp.samples):
		totreads=np.sum(newexp.data[:,idx])
		if totreads==0:
			continue
		ratio=float(numreads)/totreads
		newexp.data[:,idx]=newexp.data[:,idx]*ratio
		if fixorig:
			au.Debug(2,'fixing original frequencies')
			newexp.origreads[idx]=float(newexp.origreads[idx])/ratio
	newexp.filters.append("renormalized reads to sum %d" % numreads)
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

	newexp=copy.deepcopy(expdat)
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
	for jose clemente long babies paper
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



def findsamples(expdat,field,value,exclude=False):
	"""
	return the positions of samples in expdat matching value in field
	similar to filtersamples but returns a list of indices (for the data matrix)
	input:
	expdat
	field - name of the field to test
	value - the value to look for
	exclude - True to get positions without that value, False to get positions of the value

	output:
	pos - a list of positions matching the field/val (for use as indices in expdat.data)
	"""

	pos=[]
	for cidx,csamp in enumerate(expdat.samples):
		if expdat.smap[csamp][field]==value:
			if not exclude:
				pos.append(cidx)
		else:
			if exclude:
				pos.append(cidx)
	return pos


def BaysZeroClassifyTest(oexpdat,field,val1,val2=False,n_folds=10,istreeexpand=False,randseed=False,numiter=1):
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

	output:
	auc - the auc of each iteration
	"""
	if randseed:
		np.random.seed(randseed)

	# we want only to keep the 2 classes
	if val2:
		oexpdat=filtersamples(oexpdat,field,[val1,val2])

	# remove the validation samples and keep them in a different experiment valexp
	ind1=findsamples(oexpdat,field,val1)
	types=np.zeros(len(oexpdat.samples))
	types=types>10
	types[ind1]=True
	aucres=[]
	for citer in range(numiter):
		rs = sklearn.cross_validation.StratifiedKFold(types, n_folds=n_folds,shuffle=True)
		for trainidx,testidx in rs:
			valexp=reordersamples(oexpdat,testidx)
			expdat=reordersamples(oexpdat,trainidx)
			# classify
			lrscores=BayZeroClassify(expdat,valexp,field,val1,val2,istreeexpand)

			# prepare the correct answer list
			typeis1=[]
			for vsamp in valexp.samples:
				if valexp.smap[vsamp][field]==val1:
					vtype=1
				else:
					vtype=2
				typeis1.append(vtype==1)

			cauc=sklearn.metrics.roc_auc_score(typeis1, lrscores, average='macro', sample_weight=None)
			au.Debug(4,"auc=%f" % cauc)
			aucres.append(cauc)

	au.Debug(8,"mean is : %f" % np.mean(aucres))
	au.Debug(8,"s.err is : %f" % (np.std(aucres)/np.sqrt(len(aucres))))
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
		expdat=filtersamples(expdat,field,[val1,val2])
	ind1=findsamples(expdat,field,val1)
	types=np.zeros(len(expdat.samples))
	types=types>10
	types[ind1]=True

	# if an expanded tree, keep the best subtrees
	if istreeexpand:
		expdat=keeptreebest(expdat,field,val1,val2)
		valexp=filterseqs(valexp,expdat.seqs)

	# prepare the claissifier
	g1ind=findsamples(expdat,field,val1)
	if val2:
		g2ind=findsamples(expdat,field,val2)
	else:
		g2ind=findsamples(expdat,field,val1,exclude=True)

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
		au.Debug(2,"Classifying sample %s" % vsamp)
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
		au.Debug(2,"LRScore %f, zscore %f, nzscore %f, freqscore %f" % (totratio,crat0,cratn0,cratfreq))
	au.Debug(3,"Finished classifying")
	return lrscores


def loadexptree(expdat,treefilename):
	"""
	load a tree file associated with an experiment
	DONT USE - CANNOT DEEP COPY IT!
	input:
	expdat
	treefilename - the name of the newick tree file (from make_phylogeny.py). note that need to use the sequences as the fasta sequence ids (use -s in the CreateTable)

	output:
	expdat - with a new field - tree
	"""

	import skbio.tree

	tree=skbio.tree.TreeNode.read(treefilename)
	au.Debug(4,'Loaded tree')
	expdat.tree=tree

	return expdat



def insertbacteria(expdat,freqs=[],seq="unknown",tax="unknown"):
	"""
	insert a new bacteria to an experiment

	input:
	expdat
	freqs - the frequency of the bacteria in all samles of expdat or False to add zeros
	seq - the sequence of the new bacteria
	tax - taxonomy of the new bacteria

	output:
	pos - position of the new bacteria
	"""

	if len(freqs)==0:
		freqs=np.zeros([1,len(expdat.seqs)])

	expdat.data=np.vstack((expdat.data,freqs))
	expdat.tax.append(tax)

	if seq in expdat.seqdict:
		au.Debug(6,'Sequence already in experiment',seq)
	# get a unique sequence
		cid=0
		while seq+str(cid) in expdat.seqdict:
			cid+=1
		expdat.seqs.append()
		seq=seq+str(cid)

	expdat.seqs.append(seq)
	expdat.seqdict[seq]=len(expdat.seqs)-1
	expdat.sids.append(seq)
	return expdat,len(expdat.seqs)-1


def addsubtrees(expdat,tree,inplace=False):
	"""
	add otus for all subtrees with the frequency being the sum of all bacteria in the subtree
	input:
	expdat - the experiment
	tree - the tree for the experiment
	inplace - if true, replace current experiment

	output:
	newexp - the new experiment with twice-1 number of otus
	"""

#	if not expdat.tree:
#		au.Debug(8,"No tree loaded for experiment")
#		return False

	if inplace:
		newexp=expdat
	else:
		newexp=copy.deepcopy(expdat)

	subtrees=tree.subsets()
	for csubtree in subtrees:
		newname=""
		newtax=""
		newfreq=np.zeros([1,len(newexp.samples)])
		for cbact in csubtree:
			if not cbact in newexp.seqdict:
				au.Debug(6,'sequence not in seqdict',cbact)
				continue
			cpos=newexp.seqdict[cbact]
			newfreq+=newexp.data[cpos,:]
			newname+='%d,' % cpos
			if newtax=='':
				newtax=newexp.tax[cpos]
			else:
				newtax=au.common_start(newtax,newexp.tax[cpos])
		newexp,newpos=insertbacteria(newexp,freqs=newfreq,seq=newname,tax=newtax)
	newexp.filters.append("Add subtrees")
	return(newexp)


def keeptreebest(expdat,field,val1,val2,method="meandif"):
	"""
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

	pos1=findsamples(expdat,field,val1)
	if val2:
		pos2=findsamples(expdat,field,val2)
	else:
		pos2=findsamples(expdat,field,val1,exclude=True)
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
#				print("added %s" % str(expdat.seqdict[cseq]))
#			else:
#				print('cannot use %s' % str(expdat.seqdict[cseq]))
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
#			print("added %s" % poss)
#		else:
#			print("cannot use %s" % poss)

	newexp=reorderbacteria(expdat,keep)
	newexp.filters.append("keeptreebest field %s val1 %s val2 %s" % (field,val1,str(val2)))
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

	newexp=copy.deepcopy(expdat)
	numsamps=len(newexp.samples)
	for idx,cseq in enumerate(newexp.seqs):
		rp=np.random.permutation(numsamps)
		newexp.data[idx,:]=newexp.data[idx,rp]
	if normalize:
		newexp=normalizereads(newexp,inplace=True,fixorig=False)

	newexp.filters.append("RANDOMIZED!!! normalize = %s" % normalize)

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

	vals=getfieldvals(expdat,field)
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
		vals=getfieldvals(expdat,cfield)
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
		fval=au.fdr(justp)
		keep=np.where(np.array(fval)<=fdr)
		keep=keep[0]
	else:
		keep=np.arange(len(justp))

	if len(keep)==0:
		au.Debug(6,'No significant cateogries found')

	upv=[]
	for ckeep in keep:
		upv.append(allpv[ckeep])
		au.Debug(6,allpv[ckeep])

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
			au.Debug('method %s not supported' % method)
	si=np.argsort(effects)
	newenrich=au.reorder(enrich,si[::-1])
	return newenrich


def toorigreads(expdat,inplace=False):
	"""
	convert the number of reads to absolute using the origreads field
	input:
	expdat
	inplace - True to replace current exp, false to create a new one

	output:
	newexp - each sample has origreads reads (instead of 10k)
	"""
	if inplace:
		newexp=expdat
	else:
		newexp=copy.deepcopy(expdat)

	for idx,csamp in enumerate(newexp.samples):
		totreads=np.sum(newexp.data[:,idx])
		origreads=newexp.origreads[idx]
		if totreads==0:
			continue
		ratio=float(origreads)/totreads
		newexp.data[:,idx]=newexp.data[:,idx]*ratio
	newexp.filters.append("changed reads to origread value")
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

	newexp=filterorigreads(expdat,numreads,inplace)
	newexp=toorigreads(newexp,inplace=True)

	table=biom.table.Table(newexp.data,newexp.seqs,newexp.samples)
	table=table.subsample(numreads,axis='observation')
	tids=table.ids(axis='sample')
	for idx,cid in enumerate(tids):
		if not cid==newexp.samples[idx]:
			print('problem with sample ids!!!!')
	newpos=[]
	for cseq in table.ids(axis='observation'):
		newpos.append(newexp.seqdict[cseq])
	newexp=reorderbacteria(newexp,newpos,inplace=True)
	newexp.data=table.matrix_data.todense().A
	newexp=normalizereads(newexp,numreads=10000,inplace=True,fixorig=False)
	for cidx in range(len(newexp.samples)):
		newexp.origreads[cidx]=numreads
	newexp=updateorigreads(newexp)
	newexp.filters.append("subsample to %d" % numreads)
	return newexp


def updateorigreads(expdat):
	for idx,csamp in enumerate(expdat.samples):
		expdat.smap[csamp]['origReads']=expdat.origreads[idx]
	return expdat


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
			au.Debug(9,'testenrichment method not supported',method)
			return False
	if fdr:
		fval=au.fdr(justp)
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

def findseqsinexp(expdat,seqs):
	"""
	find sequences from seqs in expdat sequences and return the indices
	input:
	expdat
	seqs - a list of sequences

	output:
	res - a list of indices where seqs are in expdat sequences
	"""
	res=[]
	for cseq in seqs:
		res.append(expdat.seqdict[cseq])
	return res


def testbactenrichment(expdat,seqs,cdb=False,bdb=False,dbexpres=False,translatestudy=False):
	"""
	test for enrichment in bacteria database categories for the bacteria in the list seqs
	enrichment is tested against manual curation (if cdb not False) and automatic curation (if bactdb not false)

	input:
	expdat
	seqs - the sequences in the cluster
	cdb - the cooldb (manual curation) or false to skip
	bactdb - the automatic database or false to skip
	dbexpres - the assignment of values to all bacteria in the experiment (for bactdb) or false to calulate it. it is the output of bactdb.GetSeqListInfo()

	output:
	dbexpres - new if calculated
	"""

	# maybe need to keep similar freq bacteria?
	if cdb:
		cooldb.testenrichment(cdb,expdat.seqs,seqs)
	if bdb:
		if not dbexpres:
			dbexpres=bactdb.GetSeqListInfo(bdb,expdat.seqs,info='studies')
		seqpos=findseqsinexp(expdat,seqs)
		plist=testenrichment(dbexpres,seqpos,printit=False)
		for cpv in plist:
			cname=cpv['name']
			if translatestudy:
				studyname=bactdb.StudyNameFromID(bdb,cname)
			else:
				studyname=cname
			print("%s - observed %f, expected %f, pval %f" % (studyname,cpv['observed'],cpv['expected'],cpv['pval']))
	return dbexpres




def saveillitable(expdat,filename):
	"""
	create a tsv table file for use with the illi visualization tool
	input:
	expdat
	filename - name of the table file to create
	"""

	fl=open(filename,'w')
	for idx,ctax in enumerate(expdat.tax):
		atax=ctax.split(';')
		btax=[]
		for tpart in atax:
			if not tpart=='':
				btax.append(tpart)
		if len(btax)>2:
			btax=btax[-2:]
		taxstr=''
		for tpart in btax:
			taxstr+=';%s' % tpart
		cname='%s-%d' % (taxstr,idx)
		fl.write('\t%s' % cname)
	fl.write('\n')
	for sampidx,csamp in enumerate(expdat.samples):
		fl.write('%s' % csamp)
		for idx,cseq in enumerate(expdat.tax):
			fl.write('\t%f' % expdat.data[idx,sampidx])
		fl.write('\n')
	fl.close()


def sortcorrelation(expdat,method='all'):
	"""
	sort bacteria according to highest correlation/anti-correlation

	input:
	expdat
	method:
		pres - use correlation only on samples where present in one of the two sequnences
		all - use correlation on all samples (default)

	output:
	newexp - the experiment with bacteria sorted by correlation (each time next bacteria the most abs(corr) to the current bacteria)
	"""

	cdat=copy.copy(expdat.data)
	cdat[cdat<=2]=2
	cdat=np.log2(cdat)
	cdat=scale(cdat,axis=1,copy=False,with_mean=False)
	au.Debug(6,"Calculating correlation matrix")
	cmat=np.corrcoef(cdat)
	au.Debug(6,"sorting bacteria")
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
	newexp=reorderbacteria(expdat,order)
	newexp.filters.append("correlation sort")
	return newexp


def plotdiffsummary2(expdatlist1,expdatlist2,seqs1,seqs2,field,val1,val2=False,method='mean',sortit=True):
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
	diff1=plotdiffsummary(expdatlist1,seqs1,field,val1,val2,method,sortit)
	diff2=plotdiffsummary(expdatlist2,seqs2,field,val1,val2,method,sortit)
	print(np.shape(diff1))
	print(np.shape(diff2))
	diff=np.vstack([diff1,diff2])
	maxdiff=np.nanmax(np.abs(diff))
	figure()
	imshow(diff,interpolation='nearest',aspect='auto',cmap=plt.get_cmap("coolwarm"),clim=[-maxdiff,maxdiff])
	colorbar()
	title("log2 fold change between %s and %s in field %s" % (val1,val2,field))
	plot([-0.5,len(expdatlist1)-0.5],[np.shape(diff1)[0]-0.5,np.shape(diff1)[0]-0.5],'k')
	autoscale(tight=True)

def plotdiffsummary(expdatlist,seqs,field,val1,val2=False,method='mean',sortit=True):
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

	output:
	diffsum - the same as the plotted heatmap (row per otu, column per experiment)
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
	diff=np.array(getdiffsummary(expdatlist[0],seqs,field,val1[0],val2[0],method))
	odiff=copy.copy(diff)
	odiffnotnan=np.where(np.isfinite(odiff))[0]
	diffsum=[]
	for cidx,cexp in enumerate(expdatlist[1:]):
		cdiff=np.array(getdiffsummary(cexp,seqs,field,val1[cidx+1],val2[cidx+1],method))
		diff=np.vstack([diff,cdiff])
		notnan=np.where(np.isfinite(cdiff))[0]
		notnan=np.intersect1d(notnan,odiffnotnan)
		if len(notnan)>0:
			cdiffsum=float(np.sum((cdiff[notnan]>0)==(odiff[notnan]>0)))/len(notnan)
		else:
			cdiffsum=np.nan
		diffsum.append(cdiffsum)
	if sortit:
		si=np.argsort(diff[0,:])
		diff=diff[:,si]
	figure()
	maxdiff=np.nanmax(np.abs(diff))
	diff=np.transpose(diff)
	imshow(diff,interpolation='nearest',aspect='auto',cmap=plt.get_cmap("coolwarm"),clim=[-maxdiff,maxdiff])
	colorbar()
	title("log2 fold change between %s and %s in field %s" % (val1,val2,field))
	return diff

def getdiffsummary(expdat,seqs,field,val1,val2=False,method='mean'):
	"""
	plot the fold change between 2 groups in each of the sequences in seqs
	for zech chinese ibd paper
	input:
	expdat
	seqs - the sequences to examine
	field - name of the field dividing the 2 groups
	val1 - value of the field for group 1
	val2 - value of the field for group 2 or False for all the rest (not val1)
	method:
		- mean - calculate the difference in the mean of the 2 groups

	output:
	diff - a list of the difference between the 2 groups for each sequence
	"""

	pos1=findsamples(expdat,field,val1)
	if val2:
		pos2=findsamples(expdat,field,val2)
	else:
		pos2=findsamples(expdat,field,val1,exclude=True)

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
			threshold=0.1
		elif method=='binary':
			cval1=np.mean(expdat.data[seqpos,pos1]>0)
			cval2=np.mean(expdat.data[seqpos,pos2]>0)
			threshold=0.001
		else:
			au.Debug(9,"Unknown method %s for getdiff" % method)
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

def saveseqsforrdp(expdat,outfilename):
	"""
	save sequences of an experiment into a fasta file
	with the header identical to the sequence (for easy reloading of rdp taxonomy)
	input:
	expdat
	outfilename - name of the output fasta file
	"""
	fl=open(outfilename,'w')
	for cseq in expdat.seqs:
		fl.write('>'+cseq+'\n')
		fl.write(cseq+'\n')
	fl.close()


def loadrdptax(expdat,rdpfilename,fastaname=False,threshold=60):
	"""
	load rdp taxonomy (the output of download allrank in the rdp classifier website) and add to biom table
	input:
	expdat - the biom table for which the taxonomy was assigned (sequenced were saved)
	rdpfilename - name of the saved allrank rdp assignment
	fastaname - name of fasta file used for rdp assignment (if it was not from saveseqsforrdp) or False if sequences are in the header of the fasta
	threshold - the assignemt probability threshold under which to not include the assignment (for each level)
	"""

	if fastaname:
		seqs,headers=au.readfastaseqs(fastaname)
		hdict={}
		for idx,chead in enumerate(headers):
			hdict[chead]=seqs[idx]

	fl=open(rdpfilename,'r')
	for cline in fl:
		cline=cline.rstrip()
		cdat=cline.split(';')
		# skip header lines
		if len(cdat)<2:
			continue
		# check if sequence in experiment
		cseq=cdat[0]
		if fastaname:
			if cdat[0] in hdict:
				cseq=hdict[cseq]
			else:
				au.Debug(6,'sequence %s not found in fasta file' % cseq)
		if not cseq in expdat.seqdict:
			au.Debug(6,'sequence %s not found in experiment' % cseq)
			XXX()
			continue
		cpos=expdat.seqdict[cseq]
		ctax=''
		for idx in np.arange(2,len(cdat),2):
			cp=cdat[idx+1].rstrip('%')
			if float(cp)<60:
				break
			ctax+=';'
			ctax+=cdat[idx]
		expdat.tax[cpos]=ctax
	fl.close()
	expdat.filters.append("loaded rdp taxonomy from file %s" % rdpfilename)
	return(expdat)
