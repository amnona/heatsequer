#!/usr/bin/env python


"""
heatsequer plot window gui (2nd plot window) module
imported from plotwin.py when you plotexp() and set usegui=True
"""

# amnonscript

__version__ = "0.91"


import heatsequer as hs

import os
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PyQt5 import QtGui, QtCore, QtWidgets, uic
from PyQt5.QtCore import Qt
#from PyQt4 import QtGui
from PyQt5.QtWidgets import QCompleter,QMessageBox,QListWidgetItem
from PyQt5.QtCore import QStringListModel
import pickle
# for debugging - use XXX()
from pdb import set_trace as XXX


""""
for the GUI
"""



class SListWindow(QtWidgets.QDialog):
	def __init__(self,listdata=[],listname=''):
		"""
		create a list window with items in the list and the listname as specified
		input:
		listdata - the data to show in the list (a list)
		listname - name to display above the list
		"""
		super(SListWindow, self).__init__()
#		uic.loadUi('./ui/listwindow.py', self)
		uic.loadUi(os.path.join(hs.heatsequerdir,'ui/listwindow.py'), self)
#		uic.loadUi(hs.get_data_path('listwindow.py','ui'), self)
		for citem in listdata:
			self.lList.addItem(citem)
		if listname:
			self.lLabel.setText(listname)



class MyMplCanvas(FigureCanvas):
	"""Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
	def __init__(self, parent=None, width=5, height=4, dpi=100):
		self.fig = Figure(figsize=(width, height), dpi=dpi)
		self.axes = self.fig.add_subplot(111)
		# We want the axes cleared every time plot() is called
#		self.axes.hold(False)
		FigureCanvas.__init__(self, self.fig)
		self.setParent(parent)
		FigureCanvas.setSizePolicy(self,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding)
		FigureCanvas.updateGeometry(self)


class PlotGUIWindow(QtWidgets.QDialog):
	cexp=[]

	def __init__(self,expdat):
		super(PlotGUIWindow, self).__init__()
		hs.Debug(1,hs.get_data_path('plotguiwindow.py','ui'))
		uic.loadUi(os.path.join(hs.heatsequerdir,'ui/plotguiwindow.py'), self)
		self.bGetSequence.clicked.connect(self.getsequence)
		self.bExport.clicked.connect(self.export)
		self.bView.clicked.connect(self.view)
		self.bSave.clicked.connect(self.save)
		self.bDBSave.clicked.connect(self.dbsave)
		self.bEnrich.clicked.connect(self.enrich)
		self.bExpInfo.clicked.connect(self.expinfo)
		self.bSampleInfo.clicked.connect(self.sampleinfo)
		self.lCoolDB.doubleClicked.connect(self.showannotation)
		self.cSampleField.activated.connect(self.samplefield)
		self.FigureTab.currentChanged.connect(self.tabchange)
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
				# for conto in self.cexp.seqdb.OntoGraph.keys():
				self.cOntology.addItem(conto)
		self.dc=None
		self.createaddplot(useqt=True)
		# right click menu
		self.lCoolDB.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
		self.lCoolDB.customContextMenuRequested.connect(self.listItemRightClicked)

	def listItemRightClicked(self, QPos):
		self.listMenu= QtWidgets.QMenu()
		menuitem = self.listMenu.addAction("Delete annotation")
		menuitem.triggered.connect(self.menuDeleteAnnotation)

		parentPosition = self.lCoolDB.mapToGlobal(QtCore.QPoint(0, 0))
		self.listMenu.move(parentPosition + QPos)
		self.listMenu.show()

	def menuDeleteAnnotation(self):
		if len(self.lCoolDB.selectedItems())>1:
			print('more than 1 item')
		for citem in self.lCoolDB.selectedItems():
			cdetails=citem.data(Qt.UserRole)
			if cdetails is None:
				print('no details')
			else:
				print('delete id %d?' % cdetails['annotationid'])


	def showannotation(self):
		citem=self.lCoolDB.currentItem()
		cdetails=citem.data(Qt.UserRole)
		print('-----')
		print(cdetails)
		showannotationdata(cdetails)


	def createaddplot(self,useqt=True):
		"""
		create the additional figure for the ontology/line plots
		input:
		useqt : boolean
			True to embed the plot in the qtgui window, false to open a new figure window (so don't need the qtagg)
		"""
		if useqt:
			# add the matplotlib figure
			self.frame = QtWidgets.QWidget(self)
			self.dc = MyMplCanvas(self.frame, width=5, height=4, dpi=100)
			# add it to an hboxlayout to make it resize with window
			layout = QtWidgets.QHBoxLayout(self)
			layout.insertSpacing(0,250)
	#		layout.addWidget(self.dc)
	#		self.setLayout(layout)
			layout2 = QtWidgets.QVBoxLayout()
			layout.addLayout(layout2)
			layout2.addWidget(self.dc)
			self.mpl_toolbar = NavigationToolbar(self.dc, self)
			layout2.addWidget(self.mpl_toolbar)
			self.setLayout(layout)
		else:
			addfig=plt.figure()
			addax=addfig.add_subplot(1,1,1)
	#		addax.hold(False)
			self.dc=addax

	def sampleinfo(self):
		if not self.csamp:
			return
		csamp=self.cexp.samples[self.csamp]
		cmap=self.cexp.smap[csamp]
		info=[]
		for k,v in cmap.items():
			info.append(k+':'+v)
		slistwin = SListWindow(info,cmap['#SampleID'])
		slistwin.exec_()


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
		if self.dc is None:
			self.createaddplot()
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
		# is this needed?
#		self.dc.draw()
		self.dc.figure.canvas.draw_idle()


	def samplefield(self,qstr):
		cfield=str(qstr)
		self.lSampleFieldVal.setText(self.cexp.smap[self.cexp.samples[self.csamp]][cfield])

	def getsequence(self):
		seq=self.cexp.seqs[self.cseq]
		val,ok=QtWidgets.QInputDialog.getText(self,'Sequence',self.cexp.tax[self.cseq],text=seq)

	def view(self):
		slist=[]
		for cseq in self.selection:
			slist.append(self.cexp.tax[cseq]+'-'+str(self.cexp.sids[cseq]))
		val,ok=QtWidgets.QInputDialog.getItem(self,'Selected bacteria','',slist)

	def getselectedseqs(self):
		slist=[]
		for cseq in self.selection:
			slist.append(self.cexp.seqs[cseq])
		return slist

	def dbsave(self):
		"""
		save the selected list to the coolseq database
		"""
		val,ok=QtWidgets.QInputDialog.getText(self,'Save %d bacteria to coolseqDB' % len(self.selection),'Enter description')
		hs.Debug(1,ok)
		if ok:
			seqs=[]
			for cid in self.selection:
				seqs.append(self.cexp.seqs[cid])
			hs.cooldb.savecoolseqs(self.cexp,self.cexp.cdb,seqs,val)

	def enrich(self):
		"""
		check for annotation enrichment for selected sequences (compared to other sequences in this experiment)
		"""
		if not self.cexp.cdb:
			hs.Debug(8,'No cooldb loaded')
			return
		selseqs=[]
		for cid in self.selection:
			selseqs.append(self.cexp.seqs[cid])
		bmd=hs.cooldb.testenrichment(self.cexp.cdb,self.cexp.seqs,selseqs)
		# bmd=hs.annotationenrichment(self.cexp,selseqs)
		hs.Debug(6,'found %d items' % len(bmd))
		if len(bmd)>0:
			slistwin = SListWindow(listname='Enrichment')
			bmd=hs.sortenrichment(bmd)
			for cbmd in bmd:
				if cbmd['observed']<cbmd['expected']:
					ccolor=QtGui.QColor(155,0,0)
				else:
					ccolor=QtGui.QColor(0,155,0)
				item = QtWidgets.QListWidgetItem()
				item.setText("%s (p:%f o:%d e:%f)" % (cbmd['description'],cbmd['pval'],cbmd['observed'],cbmd['expected']))
				item.setForeground(ccolor)
				slistwin.lList.addItem(item)
				print("%s (p:%f o:%d e:%f)" % (cbmd['description'],cbmd['pval'],cbmd['observed'],cbmd['expected']))
			slistwin.exec_()


	def save(self):
		"""
		save the selected list to a fasta file
		"""
		fname = str(QtWidgets.QFileDialog.getSaveFileName(self, 'Save selection fasta file name','pita'))
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
				res=list(studysamples.items())
				vlens=[]
				for cv in res:
					totsamps=hs.bactdb.SamplesInStudy(self.cexp.seqdb,cv[0])
					vlens.append(float(len(cv[1]))/len(totsamps))
				sv,si=hs.isort(vlens,reverse=True)
				for cind in si:
					studyname=hs.bactdb.StudyNameFromID(self.cexp.seqdb,res[cind][0])
					self.lStudies.addItem('%s (%f)' % (studyname,vlens[cind]))
			else:
				self.lNumSamples.setText(str('%d/%dK' % (0,int(totdbsamples/1000))))
				self.lNumStudies.setText("0")
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
		self.addtocdblist(info)

	def addtocdblist(self,info):
		"""
		add to cdb list without clearing
		"""
		for cinfo in info:
			# test if the supercooldb annotation
			if type(cinfo)==tuple:
				details=cinfo[0]
				newitem=QListWidgetItem(cinfo[1])
				newitem.setData(Qt.UserRole,details)
				if details['annotationtype']=='diffexp':
					ccolor=QtGui.QColor(0,0,200)
				elif details['annotationtype']=='contamination':
					ccolor=QtGui.QColor(200,0,0)
				elif details['annotationtype']=='common':
					ccolor=QtGui.QColor(0,200,0)
				elif details['annotationtype']=='highfreq':
					ccolor=QtGui.QColor(0,200,0)
				else:
					ccolor=QtGui.QColor(0,0,0)
				newitem.setForeground(ccolor)
				self.lCoolDB.addItem(newitem)
			else:
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
			seqlist=list(self.selectionlines.keys())
		for cseq in seqlist:
			cline=self.selectionlines[cseq]
			self.plotax.lines.remove(cline[0])
			del self.selectionlines[cseq]
			self.selection.remove(cseq)
		self.plotfig.canvas.draw()
		self.lSelection.setText('%d bacteria' % len(self.selection))

	def expinfo(self):
		# get the selected sequences
		sequences=[]
		for cid in self.selection:
			sequences.append(self.cexp.seqs[cid])

		self.cexp.selectedseqs=sequences
		dbs = DBAnnotateSave(self.cexp)
		res=dbs.exec_()
		if res==QtWidgets.QDialog.Accepted:
			# fl=open('/Users/amnon/Python/git/heatsequer/db/ontologyfromid.pickle','rb')
			# ontologyfromid=pickle.load(fl)
			# fl.close()
			ontologyfromid=hs.scdb.ontologyfromid
			description=str(dbs.bdescription.text())
			# TODO: need to get primer region!!!!
			primerid='V4'
			method=str(dbs.bmethod.text())
			if method=='':
				method='na'
			submittername='Amnon Amir'
			curations=[]
			# if it is differential abundance
			for citem in qtlistiteritems(dbs.blistall):
				cdat=qtlistgetdata(citem)
				cval=cdat['value']
				ctype=cdat['type']
				if cval in ontologyfromid:
					cval=ontologyfromid[cval]
				else:
					hs.Debug(1,"item %s not found in ontologyfromid" % cval)
				curations.append((ctype,cval))
			if dbs.bdiffpres.isChecked():
				curtype='DIFFEXP'
			elif dbs.bisa.isChecked():
				curtypeval=dbs.bisatype.currentText()
				if 'Common' in curtypeval:
					curtype='COMMON'
				elif 'Contam' in curtypeval:
					curtype='CONTAMINATION'
				elif 'High' in curtypeval:
					curtype='HIGHFREQ'
				else:
					curtype='OTHER'
			else:
				hs.Debug(9,"No annotation type selected")
				return
			scdb=hs.scdb
			cdata=hs.supercooldb.finddataid(scdb,datamd5=self.cexp.datamd5,mapmd5=self.cexp.mapmd5)
			# if study not in database, ask to add some metadata for it
			if cdata is None:
				okcontinue=False
				while not okcontinue:
					hs.Debug(6,'study data info not found based on datamd5, mapmd5. need to add one!!!')
					qres=QtWidgets.QMessageBox.warning(self,"No study data","No information added about study data. Add info?",QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No,QtWidgets.QMessageBox.Cancel)
					if qres==QtWidgets.QMessageBox.Cancel:
						return
					if qres==QtWidgets.QMessageBox.No:
						cdata=hs.supercooldb.addexpdata(scdb,( ('DataMD5',self.cexp.datamd5), ('MapMD5',self.cexp.mapmd5) ) )
						okcontinue=True
					if qres==QtWidgets.QMessageBox.Yes:
						okcontinue=getstudydata(self.cexp)
						cdata=hs.supercooldb.finddataid(scdb,datamd5=self.cexp.datamd5,mapmd5=self.cexp.mapmd5)
						hs.Debug(1,'new cdata is %s' % cdata)
			hs.Debug(6,'Data found. id is %s' % cdata)
			hs.supercooldb.addannotations(scdb,expid=cdata,sequences=sequences,annotationtype=curtype,annotations=curations,submittername=submittername,description=description,method=method,primerid=primerid)
			# store the history
			try:
				hs.lastcurations.append(curations)
			except:
				hs.lastcurations=[curations]
			hs.lastdatamd5=self.cexp.datamd5


class DBStudyAnnotations(QtWidgets.QDialog):
	def __init__(self,studyid):
		super(DBStudyAnnotations, self).__init__()
		uic.loadUi(os.path.join(hs.heatsequerdir,'ui/annotationlist.py'), self)
		scdb=hs.scdb
		self.scdb=scdb
		self.studyid=studyid
		info=hs.supercooldb.getexpannotations(scdb,studyid)
		for cinfo in info:
			self.blist.addItem(cinfo)
		self.bdetails.clicked.connect(self.details)

	def details(self):
		items=self.blist.selectedItems()
		if len(items)==0:
			return
		print(str(items[0].text()))


class DBStudyInfo(QtWidgets.QDialog):
	def __init__(self,expdat):
		super(DBStudyInfo, self).__init__()
		uic.loadUi(os.path.join(hs.heatsequerdir,'ui/studyinfo.py'), self)
		scdb=hs.scdb
		self.scdb=scdb
		self.dataid=0
		dataid=hs.supercooldb.finddataid(scdb,datamd5=expdat.datamd5,mapmd5=expdat.mapmd5)
		if dataid is not None:
			info=hs.supercooldb.getexperimentinfo(scdb,dataid)
			for cinfo in info:
				qtlistadd(self.blist,cinfo[0]+':'+cinfo[1],{'fromdb':True,'type':cinfo[0],'value':cinfo[1]},color='grey')
			self.dataid=dataid
		else:
			qtlistadd(self.blist,"DataMD5:%s" % expdat.datamd5,{'fromdb':False,'type':"DataMD5",'value':expdat.datamd5},color='black')
			qtlistadd(self.blist,"MapMD5:%s" % expdat.mapmd5,{'fromdb':False,'type':"MapMD5",'value':expdat.mapmd5},color='black')
		self.bplus.clicked.connect(self.plus)
		self.bvalue.returnPressed.connect(self.plus)
		self.bminus.clicked.connect(self.minus)
		self.bannotations.clicked.connect(self.annotations)
		self.cexp=expdat
		self.setWindowTitle(self.cexp.studyname)
		self.prepstudyinfo()
		self.bvalue.setFocus()


	def keyPressEvent(self, e):
		"""
		override the enter event so will not close dialog
		"""
		e.ignore()

	def addentry(self,fromdb,ctype,value,color='black'):
		if len(ctype)>0 and len(value)>0:
			newentry='%s:%s' % (ctype,value)
			for citem in getqtlistitems(self.blist):
				if citem==newentry:
					hs.Debug(2,'item already in list %s' % newentry)
					return
			qtlistadd(self.blist,newentry,{'fromdb':False,'type':ctype,'value':value},color="black")


	def plus(self):
		ctype=str(self.btype.currentText())
		cval=str(self.bvalue.text())
		self.addentry(fromdb=False,ctype=ctype,value=cval,color='black')
		self.bvalue.setText('')

	def minus(self):
		items=self.blist.selectedItems()
		for citem in items:
			cdata=qtlistgetdata(citem)
			if cdata['fromdb']:
				print('delete from db')
			self.blist.takeItem(self.blist.row(citem))

	def annotations(self):
		dbsa = DBStudyAnnotations(self.dataid)
		dbsa.exec_()

	def prepstudyinfo(self):
		"""
		add the study info from the mapping file if available
		"""
		fieldlist=[('SRA_Study_s','sra'),('project_name_s','name'),('experiment_title','name'),('experiment_design_description','name'),('BioProject_s','sra')]
		cexp=self.cexp
		for (cfield,infofield) in fieldlist:
			if cfield in cexp.fields:
				uvals=hs.getfieldvals(cexp,cfield,ounique=True)
				if len(uvals)==1:
					self.addentry(fromdb=False,ctype=infofield,value=uvals[0].lower(),color='black')


class DBAnnotateSave(QtWidgets.QDialog):
	def __init__(self,expdat):
		super(DBAnnotateSave, self).__init__()
		print("DBAnnotateSave")
		uic.loadUi(os.path.join(hs.heatsequerdir,'ui/manualdata.py'), self)
		self.bplus.clicked.connect(self.plus)
		self.bminus.clicked.connect(self.minus)
		self.bontoinput.returnPressed.connect(self.plus)
		self.bstudyinfo.clicked.connect(self.studyinfo)
		self.bisa.toggled.connect(self.radiotoggle)
		self.bdiffpres.toggled.connect(self.radiotoggle)
		self.bisatype.currentIndexChanged.connect(self.isatypechanged)
		self.bhistory.clicked.connect(self.history)
		self.cexp=expdat
		self.lnumbact.setText(str(len(expdat.selectedseqs)))
		completer = QCompleter()
		self.bontoinput.setCompleter(completer)
		scdb=hs.scdb
		self.scdb=scdb
		self.dataid=hs.supercooldb.finddataid(scdb,datamd5=self.cexp.datamd5,mapmd5=self.cexp.mapmd5)

		model = QStringListModel()
		completer.setModel(model)
#		completer.setCompletionMode(QCompleter.InlineCompletion)
		completer.maxVisibleItems=10
		completer.setCaseSensitivity(Qt.CaseInsensitive)

		# make the completer selection also erase the text edit
		completer.activated.connect(self.cleartext,type=Qt.QueuedConnection)

		# in qt5 should work with middle complete as well...
#		completer.setFilterMode(Qt.MatchContains)
		if not hs.scdb.ontologyfromid:
			hs.scdb=hs.supercooldb.loaddbonto(hs.scdb)

		self.ontology=hs.scdb.ontology
		self.ontologyfromid=hs.scdb.ontologyfromid

		nlist=list(self.ontology.keys())
#		nlist=sorted(nlist)
		nlist=sorted(nlist, key=lambda s: s.lower())
		print("sorted ontology")

		model.setStringList(nlist)
		self.setWindowTitle(self.cexp.studyname)
		try:
			tt=hs.lastdatamd5
		except:
			hs.lastdatamd5=''
		if self.cexp.datamd5==hs.lastdatamd5:
			self.fillfromcuration(hs.lastcurations[-1],onlyall=True)

		self.prefillinfo()
		self.bontoinput.setFocus()


	def history(self):
		curtext=[]
		for cur in hs.lastcurations:
			ct=''
			for dat in cur:
				ct+=dat[0]+'-'+dat[1]+','
			curtext.append(ct)
		slistwin = SListWindow(curtext,'select curation from history')
		res=slistwin.exec_()
		if res:
			items=slistwin.lList.selectedItems()
			for citem in items:
				print(citem)
				spos=slistwin.lList.row(citem)
				print(spos)
				self.fillfromcuration(hs.lastcurations[spos],onlyall=False)


	def fillfromcuration(self,curation,onlyall=True,clearit=True):
		"""
		fill gui list from curation
		input:
		curation : from hs.lastcurations
		onlyall : bool
			True to show only curations which have ALL, False to show also HIGH/LOW
		clearit : bool
			True to remove previous curations from list, False to keep
		"""
		if clearit:
			self.blistall.clear()
		for cdat in curation:
			if onlyall:
				if cdat[0]!='ALL':
					continue
			self.addtolist(cdat[0],cdat[1])


	def radiotoggle(self):
		if self.bisa.isChecked():
			self.blow.setDisabled(True)
			self.bhigh.setDisabled(True)
		if self.bdiffpres.isChecked():
			self.blow.setEnabled(True)
			self.bhigh.setEnabled(True)

	def isatypechanged(self):
		"""
		changed the selection of isatype combobox so need to activate the isa radio button
		"""
		self.bisa.setChecked(True)

	def studyinfo(self):
		getstudydata(self.cexp)

	def keyPressEvent(self, e):
		"""
		override the enter event so will not close dialog
		"""
#		print(e.key())
		e.ignore()

	def minus(self):
		"""
		delete selected item from current list
		"""
		items=self.blistall.selectedItems()
		for citem in items:
			self.blistall.takeItem(self.blistall.row(citem))

	def cleartext(self):
		self.bontoinput.setText('')

	def plus(self):
		conto=str(self.bontoinput.text())
		cgroup=self.getontogroup()
		self.addtolist(cgroup,conto)
		self.cleartext()

	def addtolist(self,cgroup,conto):
		"""
		add an ontology term to the list

		input:
		cgroup : str
			the group (i.e. 'low/high/all')
		conto : str
			the ontology term to add
		"""
		if conto=='':
			hs.Debug(2,'no string to add to list')
			return
		print('addtolist %s %s' % (cgroup,conto))
		if conto in self.ontology:
			conto=self.ontologyfromid[self.ontology[conto]]
		else:
			hs.Debug(1,'Not in ontology!!!')
			# TODO: add are you sure... not in ontology list....

		# if item already in list, don't do anything
		for citem in qtlistiteritems(self.blistall):
			cdata=qtlistgetdata(citem)
			if cdata['value']==conto:
				hs.Debug(2,'item already in list')
				return

		if cgroup=='LOW':
			ctext="LOW:%s" % conto
			qtlistadd(self.blistall,ctext, {'type':'LOW','value':conto},color='red')
		if cgroup=='HIGH':
			ctext="HIGH:%s" % conto
			qtlistadd(self.blistall,ctext, {'type':'HIGH','value':conto},color='blue')
		if cgroup=='ALL':
			ctext="ALL:%s" % conto
			qtlistadd(self.blistall,ctext, {'type':'ALL','value':conto},color='black')

	def getontogroup(self):
		if self.ball.isChecked():
			return('ALL')
		if self.blow.isChecked():
			return('LOW')
		if self.bhigh.isChecked():
			return('HIGH')

	def prefillinfo(self):
		"""
		prefill "ALL" data fields based on mapping file
		if all samples have same info
		"""
		hs.Debug(1,'prefill info')
		ontologyfromid=self.ontologyfromid
#		fl=open('/Users/amnon/Python/git/heatsequer/db/ncbitaxontofromid.pickle','rb')
		fl=open(os.path.join(hs.heatsequerdir,'db/ncbitaxontofromid.pickle'),'rb')
		ncbitax=pickle.load(fl)
		fl.close()

		cexp=self.cexp
		for cfield in cexp.fields:
			uvals=[]
			if cfield in cexp.fields:
				uvals=hs.getfieldvals(cexp,cfield,ounique=True)
			# if we have 1 value
			if len(uvals)==1:
				cval=uvals[0]
				hs.Debug(1,'found 1 value %s' % cval)
				if cfield=='HOST_TAXID' or cfield=='host_taxid':
					hs.Debug(2,'%s field has 1 value %s' % (cfield,cval))
					# if ncbi taxonomy (field used differently)
					cval='NCBITaxon:'+cval
					if cval in ncbitax:
						hs.Debug(2,'found in ncbitax %s' % cval)
						cval=ncbitax[cval]
				else:
					# get the XXX from ENVO:XXX value
					uvalspl=cval.split(':',1)
					if len(uvalspl)>1:
						cval=uvalspl[1]
						cval=uvalspl[1]+' :'+uvalspl[0]
				if cval in self.ontology:
					cval=ontologyfromid[self.ontology[cval]]
					hs.Debug(2,'term %s found in ontologyfromid' % cval)
					conto=cval
					hs.Debug(1,'add prefill %s' % conto)
					self.addtolist('ALL',conto)
				else:
					hs.Debug(3,'term %s NOT found in ontologyfromid' % uvals[0])

			else:
				hs.Debug(1,'found %d values' % len(uvals))


def getqtlistitems(qtlist):
	"""
	get a list of strings of the qtlist
	input:
	qtlist : QTListWidget

	output:
	item : list of str
	"""
	items = []
	for index in range(qtlist.count()):
		items.append(str(qtlist.item(index).text()))
	return items



def qtlistadd(qtlist,text,data,color="black"):
	"""
	Add an entry (text) to qtlist and associaxte metadata data
	input:
	qtlist : QTListWidget
	text : str
		string to add to list
	data : arbitrary python var
		the data to associate with the item (get it by qtlistgetdata)
	color : (R,G,B)
		the color of the text in the list
	"""
	item = QtWidgets.QListWidgetItem()
	item.setText(text)
	ccol=QtGui.QColor()
	ccol.setNamedColor(color)
	item.setForeground(ccol)
	item.setData(Qt.UserRole,data)
	qtlist.addItem(item)


def qtlistgetdata(item):
	"""
	Get the metadata associated with item as position pos
	input:
	qtlist : QtListWidget
	index : QtListWidgetItem
		the item to get the info about

	output:
	data : arbitrary
		the data associated with the item (using qtlistadd)
	"""
#	item=qtlist.item(index)
	if sys.version_info[0] < 3:
		# QVariant version 1 API (python2 default)
		data=item.data(Qt.UserRole).toPyObject()
	else:
		# QVariant version 2 API (python3 default)
		data=item.data(Qt.UserRole)
	return data


def qtlistiteritems(qtlist):
	"""
	iterate all items in a list
	input:
	qtlist : QtListWidget
	"""
	for i in range(qtlist.count()):
		yield qtlist.item(i)


def getstudydata(cexp):
	"""
	open the study info window and show/get new references for the study data

	input:
	cexp : Experiment
		the experiment for which to show the data (uses the datamd5 and mapmd5)

	output:
	hasdata : Bool
		True if the study has data, False if not
	"""
	dbsi = DBStudyInfo(cexp)
	res=dbsi.exec_()
	if res==QtWidgets.QDialog.Accepted:
		newstudydata=[]
		allstudydata=[]
		for citem in qtlistiteritems(dbsi.blist):
			cdata=qtlistgetdata(citem)
			allstudydata.append( (cdata['type'],cdata['value']) )
			if cdata['fromdb']==False:
				newstudydata.append( (cdata['type'],cdata['value']) )

		if len(newstudydata)==0:
			hs.Debug(6,'No new items. not saving anything')
			return True
		# look if study already in table
		cid=hs.supercooldb.finddataid(dbsi.scdb,datamd5=cexp.datamd5,mapmd5=cexp.mapmd5)
		if cid is None:
			hs.Debug(6,'no studyid found for datamd5 %s, mapmd5 %s' % (cexp.datamd5,cexp.mapmd5))
#			cdata=hs.supercooldb.addexpdata(scdb,( ('DataMD5',cexp.datamd5), ('MapMD5',cexp.mapmd5) ) )
			hs.Debug(3,'Adding to new experiment')
		dataid=hs.supercooldb.addexpdata(dbsi.scdb,newstudydata,studyid=cid)
		hs.Debug(6,'Study data saved to id %d' % dataid)
		if len(allstudydata)>2:
			return True
	return False



def showannotationdata(annotationdetails):
	"""
	show the list of annotation details and the sequences associated with it

	intput:
	annotationdetails : dict
		dict of various fields of the annotation (includeing annotationid)
		from scdb.getannotationstrings()
	cexp : experiment
		the experiment (for rhe scdb pointer)
	"""
	info=[]
	if annotationdetails is None:
		return
	for k,v in annotationdetails.items():
		if type(v)==list:
			for cv in v:
				info.append('%s:%s' % (k,cv))
		else:
			info.append('%s:%s' % (k,v))
	# get the annotation sequences:
	if 'annotationid' in annotationdetails:
		seqs=hs.supercooldb.getannotationseqs(hs.scdb,annotationdetails['annotationid'])
		info.append('sequences: %d' % len(seqs))
	slistwin = SListWindow(info,'Annotation details')
	slistwin.exec_()
