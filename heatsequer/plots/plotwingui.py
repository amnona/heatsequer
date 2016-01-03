#!/usr/bin/env python


"""
heatsequer plot window gui (2nd plot window) module
imported from plotwin.py when you plotexp() and set usegui=True
"""

# amnonscript

__version__ = "0.91"


import heatsequer as hs

import numpy as np
import matplotlib as mpl
mpl.use('Qt4Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from PyQt4 import QtGui, QtCore, uic
# for debugging - use XXX()
from pdb import set_trace as XXX


""""
for the GUI
"""


class MyMplCanvas(FigureCanvas):
	"""Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
	def __init__(self, parent=None, width=5, height=4, dpi=100):
		self.fig = Figure(figsize=(width, height), dpi=dpi)
		self.axes = self.fig.add_subplot(111)
		# We want the axes cleared every time plot() is called
#		self.axes.hold(False)
		FigureCanvas.__init__(self, self.fig)
		self.setParent(parent)
		FigureCanvas.setSizePolicy(self,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding)
		FigureCanvas.updateGeometry(self)


class PlotGUIWindow(QtGui.QDialog):
	cexp=[]

	def __init__(self,expdat):
		super(PlotGUIWindow, self).__init__()
		uic.loadUi('ui/plotguiwindow.py', self)
		self.bGetSequence.clicked.connect(self.getsequence)
		self.bExport.clicked.connect(self.export)
		self.bView.clicked.connect(self.view)
		self.bSave.clicked.connect(self.save)
		self.bDBSave.clicked.connect(self.dbsave)
		self.bEnrich.clicked.connect(self.enrich)
		self.connect(self.cSampleField, QtCore.SIGNAL('activated(QString)'), self.samplefield)
		self.FigureTab.connect(self.FigureTab, QtCore.SIGNAL("currentChanged(int)"),self.tabchange)
		self.cSampleField.setCurrentIndex(0)
		self.cexp=expdat
		self.selectionlines={}
		self.selection=[]

		self.setWindowTitle(self.cexp.studyname)

		for cfield in self.cexp.fields:
			self.cSampleField.addItem(cfield)
			self.cPlotXField.addItem(cfield)

		if self.cexp.seqdb:
			ontofields,ontonames=hs.bactdb.getontonames(self.cexp.seqdb)
			for conto in ontofields:
#			for conto in self.cexp.seqdb.OntoGraph.keys():
				self.cOntology.addItem(conto)


		# add the matplotlib figure
		self.frame = QtGui.QWidget(self)
		self.dc = MyMplCanvas(self.frame, width=5, height=4, dpi=100)
		# add it to an hboxlayout to make it resize with window
		layout = QtGui.QHBoxLayout(self)
		layout.insertSpacing(0,250)
#		layout.addWidget(self.dc)
#		self.setLayout(layout)
		layout2 = QtGui.QVBoxLayout()
		layout.addLayout(layout2)
		layout2.addWidget(self.dc)
		self.mpl_toolbar = NavigationToolbar(self.dc, self)
		layout2.addWidget(self.mpl_toolbar)
		self.setLayout(layout)


	def tabchange(self,newtab):
		hs.Debug(0,"new tab",newtab)
		if newtab==2:
			self.plotxgraph()
		if newtab==1:
			self.plotontology()

	def plotontology(self):
		if self.cexp.seqdb:
			self.dc.axes.clear()
			hs.Debug(2,"plotting taxonomy for seq %s onto %s" % (self.cexp.seqs[self.cseq],self.cexp.ontofigname))
			hs.bactdb.PlotOntologyGraph(self.cexp.seqdb,self.cexp.seqs[self.cseq],field=str(self.cOntology.currentText()),toax=self.dc.axes)
			self.dc.draw()


	def plotxgraph(self):
		self.dc.axes.clear()
		seqs=self.getselectedseqs()
		if self.cPlotNormalizeY.checkState()==0:
			normalizey=False
		else:
			normalizey=True
		if self.cPlotXNumeric.checkState()==0:
			xfield=False
		else:
			xfield=str(self.cPlotXField.currentText())
		hs.plotseqfreq(self.cexp,seqs=seqs,toaxis=self.dc.axes,normalizey=normalizey,xfield=xfield)
		self.dc.draw()

	def samplefield(self,qstr):
		cfield=str(qstr)
		self.lSampleFieldVal.setText(self.cexp.smap[self.cexp.samples[self.csamp]][cfield])

	def getsequence(self):
		seq=self.cexp.seqs[self.cseq]
		val,ok=QtGui.QInputDialog.getText(self,'Sequence',self.cexp.tax[self.cseq],text=seq)

	def view(self):
		slist=[]
		for cseq in self.selection:
			slist.append(self.cexp.tax[cseq]+'-'+str(self.cexp.sids[cseq]))
		val,ok=QtGui.QInputDialog.getItem(self,'Selected bacteria','',slist)

	def getselectedseqs(self):
		slist=[]
		for cseq in self.selection:
			slist.append(self.cexp.seqs[cseq])
		return slist

	def dbsave(self):
		"""
		save the selected list to the coolseq database
		"""
		val,ok=QtGui.QInputDialog.getText(self,'Save %d bacteria to coolseqDB' % len(self.selection),'Enter description')
		print(ok)
		if ok:
			seqs=[]
			for cid in self.selection:
				seqs.append(self.cexp.seqs[cid])
			hs.cooldb.savecoolseqs(self.cexp,self.cexp.cdb,seqs,val)

	def enrich(self):
		"""
		check for annotation enrichment for selected sequences (compared to other sequences in this experiment)
		"""
		if self.cexp.cooldb:
			selseqs=[]
			for cid in self.selection:
				selseqs.append(self.cexp.seqs[cid])
			bmd=hs.cooldb.testenrichment(self.cexp.cooldb,self.cexp.seqs,selseqs)
			bmd=hs.sortenrichment(bmd)
			for cbmd in bmd:
				# if cbmd['observed']<cbmd['expected']:
				# 	ccolor=QtGui.QColor(155,0,0)
				# else:
				# 	ccolor=QtGui.QColor(0,155,0)
				# item = QtGui.QListWidgetItem()
				# item.setText("%s (p:%f o:%d e:%f)" % (cbmd['description'],cbmd['pval'],cbmd['observed'],cbmd['expected']))
				# item.setTextColor(ccolor)
				# self.lBacteria.addItem(item)
				print("%s (p:%f o:%d e:%f)" % (cbmd['description'],cbmd['pval'],cbmd['observed'],cbmd['expected']))


	def save(self):
		"""
		save the selected list to a fasta file
		"""
		fname = str(QtGui.QFileDialog.getSaveFileName(self, 'Save selection fasta file name','pita'))
		slist=[]
		for cseq in self.selection:
			slist.append(self.cexp.seqs[cseq])
		hs.saveseqsfasta(self.cexp,slist,fname)
		hs.Debug(6,'Saved %d sequences to file %s' % (len(slist),fname))

	def export(self):
		"""
		export the selected bacteria list to the global variable 'selectlist'
		"""
		global selectlist

		hs.Debug(0,'exporting')
		selectlist=[]
		for cseq in self.selection:
			selectlist.append(self.cexp.seqs[cseq])

	def updateinfo(self,csamp,cseq):
		"""
		update the information about the sample/bacteria
		"""
		self.csamp=csamp
		self.cseq=cseq
		self.lSample.setText(self.cexp.samples[self.csamp])
		self.lTaxonomy.setText(self.cexp.tax[self.cseq])
		self.lID.setText(str(self.cexp.sids[self.cseq]))
		self.lReads.setText('%f' % (float(self.cexp.data[self.cseq,self.csamp])/100))
		self.lSampleFieldVal.setText(self.cexp.smap[self.cexp.samples[self.csamp]][str(self.cSampleField.currentText())])
		# update the stats about the database:
		if self.cexp.seqdb:
			self.lStudies.clear()
			totappear,numstudies,allstudies,studysamples,totdbsamples=hs.bactdb.GetSeqInfo(self.cexp.seqdb,self.cexp.seqs[self.cseq])
			if totappear>0:
				self.lNumSamples.setText(str('%d/%dK' % (totappear,int(totdbsamples/1000))))
				self.lNumStudies.setText(str(numstudies))
				res=studysamples.items()
				vlens=[]
				for cv in res:
					totsamps=hs.bactdb.SamplesInStudy(self.cexp.seqdb,cv[0])
					vlens.append(float(len(cv[1]))/len(totsamps))
				sv,si=hs.isort(vlens,reverse=True)
				for cind in si:
					studyname=hs.bactdb.StudyNameFromID(self.cexp.seqdb,res[cind][0])
					self.lStudies.addItem('%s (%f)' % (studyname,vlens[cind]))
		if self.FigureTab.currentIndex()==2:
			self.plotxgraph()
		if self.FigureTab.currentIndex()==1:
			self.plotontology()

	def updatecdb(self,info):
		"""
		update the coolseq database info for the bacteria
		by adding all lines in list to the listbox
		"""
		self.lCoolDB.clear()
		for cinfo in info:
			self.lCoolDB.addItem(cinfo)

	def selectbact(self,bactlist,flip=True):
		"""
		add bacteria from the list bactlist (position in exp) to the selection
		flip - if true, if bacteria from list is already in selection, remove it
		"""
		for cseq in bactlist:
			# if already in list and can flip, remove from list instead
			if flip:
				if cseq in self.selectionlines:
					self.clearselection([cseq])
					hs.Debug(0,'Flip')
					continue
			if cseq in self.selectionlines:
				continue
			cline=self.plotax.plot([-0.5,len(self.cexp.samples)-0.5],[cseq,cseq],':w')
			self.selectionlines[cseq]=cline
			self.selection.append(cseq)
		self.plotfig.canvas.draw()
		self.lSelection.setText('%d bacteria' % len(self.selection))

	def clearselection(self,seqlist=False):
		if not seqlist:
			seqlist=self.selectionlines.keys()
		for cseq in seqlist:
			cline=self.selectionlines[cseq]
			self.plotax.lines.remove(cline[0])
			del self.selectionlines[cseq]
			self.selection.remove(cseq)
		self.plotfig.canvas.draw()
		self.lSelection.setText('%d bacteria' % len(self.selection))
