#!/usr/bin/env python


"""
heatsequer full gui
"""

# amnonscript

__version__ = "0.91"

import heatsequer as hs


import sys
from PyQt4 import QtGui, uic, QtCore
# import cPickle as pickle
import cPickle as pickle
import matplotlib
matplotlib.use('Qt4Agg')
import os.path
import numpy as np


class JoinWindow(QtGui.QDialog):
	cexp=[]

	def __init__(self,expdat):
		super(JoinWindow, self).__init__()
		uic.loadUi('./ui/joinfields.py', self)
		self.cexp=expdat
		self.cField1.addItems(expdat.fields)
		self.cField2.addItems(expdat.fields)
#		self.bOK.clicked.connect(self.OK)
#		self.bCancel.clicked.connect(self.cancel)



class MetaDataDetailsWindow(QtGui.QDialog):
	cexp=[]

	def __init__(self,expdat):
		super(MetaDataDetailsWindow, self).__init__()
		uic.loadUi('./ui/metadatadetails.py', self)
		self.cexp=expdat
		self.bFieldValues.clicked.connect(self.values)
		self.cField.addItems(expdat.fields)


	def values(self):
		cfield=str(self.cField.currentText())
		val,ok=QtGui.QInputDialog.getItem(self,'Select field value','Field=%s' % cfield,list(set(hs.getfieldvals(self.cexp,cfield))))
		if ok:
			self.tValue.setText(val)


class MetaDataWindow(QtGui.QDialog):
	cexp=[]

	def __init__(self,expdat):
		super(MetaDataWindow, self).__init__()
		uic.loadUi('./ui/plotmetadata.py', self)
		self.cexp=expdat
		self.mddict={}
		for cmeta in expdat.plotmetadata:
			cname="%s;%s;%s;%s;%s" % (cmeta[0],cmeta[1],cmeta[2],cmeta[3],cmeta[4])
			self.lMetaData.addItem(cname)
			self.mddict[cname]=cmeta
		self.bAdd.clicked.connect(self.add)
		# the main list right mouse menu
		self.lMetaData.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
		self.lMetaData.connect(self.lMetaData, QtCore.SIGNAL("customContextMenuRequested(QPoint)"),self.listItemRightClicked)

	def listItemRightClicked(self, QPos):
		self.listMenu= QtGui.QMenu()
		menuremove = self.listMenu.addAction("Remove Item")
		self.connect(menuremove, QtCore.SIGNAL("triggered()"), self.menuRemove)
		parentPosition = self.lMetaData.mapToGlobal(QtCore.QPoint(0, 0))
		self.listMenu.move(parentPosition + QPos)
		self.listMenu.show()

	def menuRemove(self):
		if len(self.lMetaData.selectedItems())>1:
			if QtGui.QMessageBox.warning(self,"Remove metadata?","Remove %d limes?" % len(self.lMetaData.selectedItems()),QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)==QtGui.QMessageBox.No:
				return
		for currentItemName in self.lMetaData.selectedItems():
			currentItemName=str(currentItemName.text())
#		currentItemName=str(self.bMainList.currentItem().text())
			self.removemetadata(currentItemName)

	def removemetadata(self,mdname):
		del self.mddict[mdname]
		items=self.lMetaData.findItems(mdname,QtCore.Qt.MatchExactly)
		for item in items:
			self.lMetaData.takeItem(self.lMetaData.row(item))

	def add(self):
		metadatadetailswin = MetaDataDetailsWindow(self.cexp)
		res=metadatadetailswin.exec_()
		if res==QtGui.QDialog.Accepted:
			cmeta=[]
			cmeta.append(str(metadatadetailswin.cField.currentText()))
			cmeta.append(str(metadatadetailswin.tValue.text()))
			cmeta.append(str(metadatadetailswin.tColor.text()))
			cmeta.append(metadatadetailswin.cInverse.checkState())
			cmeta.append(metadatadetailswin.cBefore.checkState())
			cname="%s;%s;%s;%s;%s" % (cmeta[0],cmeta[1],cmeta[2],cmeta[3],cmeta[4])
			self.lMetaData.addItem(cname)
			self.mddict[cname]=cmeta


class BiClusterWindow(QtGui.QDialog):
	cexp=[]

	def __init__(self,expdat,cdb=False,bdb=False):
		super(BiClusterWindow, self).__init__()
		uic.loadUi('./ui/bicluster.py', self)
		self.cexp=expdat
		self.cooldb=cdb
		self.bactdb=bdb
		self.lStudy.setText(self.cexp.studyname)
		self.bBiCluster.clicked.connect(self.bicluster)
		self.bView.clicked.connect(self.view)

	def bicluster(self):
		newexp,seqs,samples=hs.bicluster(self.cexp,method='binary',sampkeep=0,bactkeep=0,justcount=True)
		self.samples=samples[0]
		self.seqs=seqs[0]
		self.lNumSamples.setText(str(len(self.samples)))
		self.lNumBacteria.setText(str(len(self.seqs)))
		self.lSamples.clear()
		self.lBacteria.clear()

		# get the new bacteria and sample order
		allsamp=np.arange(len(self.cexp.samples))
		allbact=np.arange(len(self.cexp.seqs))

		x=np.setdiff1d(allsamp,self.samples)
		sampo=np.concatenate((self.samples,x))
		ubact=[]
		for cseq in self.seqs:
			ubact.append(self.cexp.seqdict[cseq])
		bacto=np.concatenate((ubact,np.setdiff1d(allbact,ubact)))
		self.cexp.bactorder=bacto
		self.cexp.samporder=sampo

		# if we have enough bacteria and samples, update the info list
		if len(self.seqs)>=5 and len(self.samples)>=5:
			sampmd=hs.testmdenrichmentall(self.cexp,self.samples)
			for cmd in sampmd:
				if cmd['observed']<cmd['expected']:
					ccolor=QtGui.QColor(155,0,0)
				else:
					ccolor=QtGui.QColor(0,155,0)
				item = QtGui.QListWidgetItem()
				item.setText("%s - %s (p:%f o:%d e:%f)" % (cmd['field'],cmd['val'],cmd['pval'],cmd['observed'],cmd['expected']))
				item.setTextColor(ccolor)
				self.lSamples.addItem(item)
			if self.cooldb:
				bmd=hs.cooldb.testenrichment(self.cooldb,self.cexp.seqs,self.seqs)
				bmd=hs.sortenrichment(bmd)
				for cbmd in bmd:
					if cbmd['observed']<cbmd['expected']:
						ccolor=QtGui.QColor(155,0,0)
					else:
						ccolor=QtGui.QColor(0,155,0)
					item = QtGui.QListWidgetItem()
					item.setText("%s (p:%f o:%d e:%f)" % (cbmd['description'],cbmd['pval'],cbmd['observed'],cbmd['expected']))
					item.setTextColor(ccolor)
					self.lBacteria.addItem(item)
	#					self.lBacteria.addItem("%s (p:%f o:%d e:%f)" % (cbmd['description'],cbmd['pval'],cbmd['observed'],cbmd['expected']))


	def view(self):
		cexp=self.cexp
		allsamp=np.arange(len(cexp.samples))
		allbact=np.arange(len(cexp.seqs))

		x=np.setdiff1d(allsamp,self.samples)
		sampo=np.concatenate((self.samples,x))
		ubact=[]
		for cseq in self.seqs:
			ubact.append(cexp.seqdict[cseq])
		bacto=np.concatenate((ubact,np.setdiff1d(allbact,ubact)))

		newexp=hs.reorderbacteria(cexp,bacto)
		newexp=hs.reordersamples(newexp,sampo,inplace=True)
		hs.plotexp(newexp,seqdb=self.bactdb,sortby=False,numeric=False,usegui=True,cdb=self.cooldb,showline=False)


class AdvPlotWindow(QtGui.QDialog):
	cexp=[]

	def __init__(self,expdat):
		super(AdvPlotWindow, self).__init__()
		uic.loadUi('./ui/advplot.py', self)
		self.cexp=expdat
		self.cField.addItems(expdat.fields)
		self.bOK.clicked.connect(self.OK)
		self.bFieldValues.clicked.connect(self.fieldvalues)
		self.bMetaData.clicked.connect(self.metadata)
		self.cSort.stateChanged.connect(self.usesort)
		self.bCancel.clicked.connect(self.cancel)

	def metadata(self):
		metadatawin = MetaDataWindow(self.cexp)
		res=metadatawin.exec_()
		if res==QtGui.QDialog.Accepted:
			self.cexp.plotmetadata=metadatawin.mddict.values()

	def cancel(self):
		self.reject()

	def OK(self):
		self.accept()

	def fieldvalues(self):
		cfield=str(self.cField.currentText())
		val,ok=QtGui.QInputDialog.getItem(self,'Select field value','Field=%s' % cfield,list(set(hs.getfieldvals(self.cexp,cfield))))

	def usesort(self):
		# unchecked
		if self.cSort.checkState()==0:
			self.cField.setEnabled(False)
		else:
			self.cField.setEnabled(True)



class SortSamplesWindow(QtGui.QDialog):
	cexp=[]

	def __init__(self,expdat):
		super(SortSamplesWindow, self).__init__()
		uic.loadUi('./ui/sortsamples.py', self)
		self.cexp=expdat
		self.cField.addItems(expdat.fields)
		self.bOK.clicked.connect(self.OK)
		self.bFieldValues.clicked.connect(self.fieldvalues)
		self.cOverwrite.stateChanged.connect(self.overwrite)
		self.bCancel.clicked.connect(self.cancel)
		self.tNewName.setText(self.cexp.studyname+'_s')

	def cancel(self):
		self.reject()

	def OK(self):
		self.accept()

	def fieldvalues(self):
		cfield=str(self.cField.currentText())
		val,ok=QtGui.QInputDialog.getItem(self,'Select field value','Field=%s' % cfield,list(set(hs.getfieldvals(self.cexp,cfield))))

	def overwrite(self):
		# unchecked
		if self.cOverwrite.checkState()==0:
			self.tNewName.setText(self.cexp.studyname+'_s')
			self.tNewName.setReadOnly(False)
		else:
			self.tNewName.setText(self.cexp.studyname)
			self.tNewName.setReadOnly(True)


class DiffExpWindow(QtGui.QDialog):
	cexp=[]

	def __init__(self,expdat):
		super(DiffExpWindow, self).__init__()
		uic.loadUi('./ui/diffexp.py', self)
		self.cexp=expdat
		self.cField.addItems(expdat.fields)
		self.bFieldValues1.clicked.connect(self.fieldvalues1)
		self.bFieldValues2.clicked.connect(self.fieldvalues2)
		self.cAll.stateChanged.connect(self.allvalues)
#		self.cField.stateChanged.connect(self.fieldchanged)
		self.tNewName.setText(self.cexp.studyname+'_de')

#	def fieldchanged(self):
#		self.tNewName.setText(self.cexp.studyname+'_'+str(self.cField.currentText())+'_de')

	def fieldvalues1(self):
		cfield=str(self.cField.currentText())
		val,ok=QtGui.QInputDialog.getItem(self,'Select field value','Field=%s' % cfield,list(set(hs.getfieldvals(self.cexp,cfield))))
		if ok:
			self.tValue1.setText(val)

	def fieldvalues2(self):
		cfield=str(self.cField.currentText())
		val,ok=QtGui.QInputDialog.getItem(self,'Select field value','Field=%s' % cfield,list(set(hs.getfieldvals(self.cexp,cfield))))
		if ok:
			self.tValue2.setText(val)

	def allvalues(self):
		# unchecked
		if self.cAll.checkState()==0:
			self.tValue2.setReadOnly(False)
		else:
			self.tValue2.setReadOnly(True)


class ClassifyWindow(QtGui.QDialog):
	cexp=[]

	def __init__(self,expdat):
		super(ClassifyWindow, self).__init__()
		uic.loadUi('./ui/classifier.py', self)
		self.cexp=expdat
		self.cField.addItems(expdat.fields)
		self.bFieldValues1.clicked.connect(self.fieldvalues1)
		self.bFieldValues2.clicked.connect(self.fieldvalues2)
		self.cAll.stateChanged.connect(self.allvalues)
#		self.cField.stateChanged.connect(self.fieldchanged)
		self.tNewName.setText(self.cexp.studyname+'_de')

#	def fieldchanged(self):
#		self.tNewName.setText(self.cexp.studyname+'_'+str(self.cField.currentText())+'_de')

	def fieldvalues1(self):
		cfield=str(self.cField.currentText())
		val,ok=QtGui.QInputDialog.getItem(self,'Select field value','Field=%s' % cfield,list(set(hs.getfieldvals(self.cexp,cfield))))
		if ok:
			self.tValue1.setText(val)

	def fieldvalues2(self):
		cfield=str(self.cField.currentText())
		val,ok=QtGui.QInputDialog.getItem(self,'Select field value','Field=%s' % cfield,list(set(hs.getfieldvals(self.cexp,cfield))))
		if ok:
			self.tValue2.setText(val)

	def allvalues(self):
		# unchecked
		if self.cAll.checkState()==0:
			self.tValue2.setReadOnly(False)
		else:
			self.tValue2.setReadOnly(True)


class FilterSamplesWindow(QtGui.QDialog):
	cexp=[]

	def __init__(self,expdat):
		super(FilterSamplesWindow, self).__init__()
		uic.loadUi('./ui/filtersamples.py', self)
		self.cexp=expdat
		self.cField.addItems(expdat.fields)
		self.bOK.clicked.connect(self.OK)
		self.bFieldValues.clicked.connect(self.fieldvalues)
		self.cOverwrite.stateChanged.connect(self.overwrite)
		self.bCancel.clicked.connect(self.cancel)
		self.tNewName.setText(self.cexp.studyname+'_f')

	def cancel(self):
		self.reject()

	def OK(self):
		self.accept()

	def fieldvalues(self):
		cfield=str(self.cField.currentText())
		val,ok=QtGui.QInputDialog.getItem(self,'Select field value','Field=%s' % cfield,list(set(hs.getfieldvals(self.cexp,cfield))))
		if ok:
			self.tValue.setText(val)

	def overwrite(self):
		# unchecked
		if self.cOverwrite.checkState()==0:
			self.tNewName.setText(self.cexp.studyname+'_f')
			self.tNewName.setReadOnly(False)
		else:
			self.tNewName.setText(self.cexp.studyname)
			self.tNewName.setReadOnly(True)



class FilterFastaWindow(QtGui.QDialog):
	def __init__(self):
		super(FilterFastaWindow, self).__init__()
		uic.loadUi('./ui/filterfasta.py', self)
		self.bBrowse.clicked.connect(self.browse)

	def browse(self):
		fname = str(QtGui.QFileDialog.getOpenFileName(self, 'Open fasta file',filter='*.fa;;*.fasta;;*.fna'))
		self.tFileName.setText(fname)


class LoadWindow(QtGui.QDialog):
	def __init__(self):
		super(LoadWindow, self).__init__()
		uic.loadUi('./ui/load.py', self)
		self.bLoadBrowseTable.clicked.connect(self.browsetable)
		self.bLoadBrowseMap.clicked.connect(self.browsemap)
		self.bLoadLoad.clicked.connect(self.load)
		self.bLoadCancel.clicked.connect(self.cancel)

	def browsemap(self):
		fname = str(QtGui.QFileDialog.getOpenFileName(self, 'Open map file',filter='*.txt'))
		self.tLoadMap.setText(fname)


	def browsetable(self):
		fname = str(QtGui.QFileDialog.getOpenFileName(self, 'Open table file',''))
		self.tLoadTable.setText(fname)
		pname=os.path.basename(fname)
		self.tLoadName.setText(pname)

	def load(self):
		self.accept()


	def cancel(self):
		self.reject()



class ListWindow(QtGui.QDialog):
	def __init__(self,listdata=[],listname=''):
		"""
		create a list window with items in the list and the listname as specified
		input:
		listdata - the data to show in the list (a list)
		listname - name to display above the list
		"""
		super(ListWindow, self).__init__()
		uic.loadUi('./ui/listwindow.py', self)
		for citem in listdata:
			self.lList.addItem(citem)
		if listname:
			self.lLabel.setText(listname)



class AppWindow(QtGui.QMainWindow):
	# the experiments loaded for analysis
	explist={}

	def __init__(self):
		super(AppWindow, self).__init__()
		uic.loadUi('./ui/appwindow.py', self)
		self.bMainLoadNew.clicked.connect(self.load)
		self.bPickleLoad.clicked.connect(self.pickleload)
		self.bMainPlot.clicked.connect(self.plot)
		self.bMainAdvancedPlot.clicked.connect(self.advplot)
		self.bMainFilterSamples.clicked.connect(self.filtersamples)
		self.bDiffExp.clicked.connect(self.diffexp)
		self.bClassifier.clicked.connect(self.classify)
		self.bFilterOrigReads.clicked.connect(self.filterorigreads)
		self.bMainSortSamples.clicked.connect(self.sortsamples)
		self.bMainClusterBacteria.clicked.connect(self.clusterbacteria)
		self.bMainFilterMinReads.clicked.connect(self.filterminreads)
		self.bMainFilterTaxonomy.clicked.connect(self.filtertaxonomy)
		self.bFilterPresence.clicked.connect(self.filterpresence)
		self.bFilterMean.clicked.connect(self.filtermean)
		self.bSortAbundance.clicked.connect(self.sortabundance)
		self.bJoinFields.clicked.connect(self.joinfields)
		self.bJoinExps.clicked.connect(self.joinexps)
		self.bBicluster.clicked.connect(self.bicluster)
		self.bEnrichment.clicked.connect(self.enrichment)
		self.bFilterFasta.clicked.connect(self.filterfasta)
		self.bRenormalize.clicked.connect(self.renormalize)
		self.bSubsample.clicked.connect(self.subsample)


		# the main list right mouse menu
		self.bMainList.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
		self.bMainList.connect(self.bMainList, QtCore.SIGNAL("customContextMenuRequested(QPoint)"),self.listItemRightClicked)
		self.show()

		self.bactdb=False
		self.cooldb=False
		try:
			hs.Debug(6,'Loading sequence database')
			self.bactdb=pickle.load(open('dbtest2.pickle'))
			self.bactdb=hs.bactdb.dbconnect(self.bactdb)
		except:
			hs.Debug(9,'BactDB file not found')
		try:
			hs.Debug(6,'Loading coolseq database')
			self.cooldb=hs.cooldb.loaddb()
		except:
			hs.Debug(9,'CoolDB file not found')

#		self.cooldb=pickle.load(open('cooldb.pickle'))


		# put some experiments:
		hs.Debug(6,'Loading sample experiment')
		try:
			expdat=hs.load('./test_data/bears.clean.new.withtax.biom','./test_data/map.txt')
			self.addexp(expdat)
		except:
			hs.Debug(6,'Sample experiment not found. sorry')


	def listItemRightClicked(self, QPos):
		self.listMenu= QtGui.QMenu()
		menuremove = self.listMenu.addAction("Rename")
		self.connect(menuremove, QtCore.SIGNAL("triggered()"), self.menuRename)
		menuremove = self.listMenu.addAction("Delete")
		self.connect(menuremove, QtCore.SIGNAL("triggered()"), self.menuRemove)
		menusave = self.listMenu.addAction("Save (pickle) Item")
		self.connect(menusave, QtCore.SIGNAL("triggered()"), self.menuSave)
		menuexport = self.listMenu.addAction("Save (biom) Item")
		self.connect(menuexport, QtCore.SIGNAL("triggered()"), self.menuExport)
		menuinfo = self.listMenu.addAction("Info")
		self.connect(menuinfo, QtCore.SIGNAL("triggered()"), self.expinfo)
		parentPosition = self.bMainList.mapToGlobal(QtCore.QPoint(0, 0))
		self.listMenu.move(parentPosition + QPos)
		self.listMenu.show()

	def expinfo(self):
		items=self.bMainList.selectedItems()
		if len(items)!=1:
			print("Need 1 item")
			return
		for citem in items:
			cname=str(citem.text())
			cexp=self.explist[cname]
			listwin = ListWindow(cexp.filters,cexp.studyname)
			res=listwin.exec_()


	def menuRename(self):
		if len(self.bMainList.selectedItems())>1:
			return
		for citem in self.bMainList.selectedItems():
			cname=str(citem.text())
			cexp=self.explist[cname]
			val,ok=QtGui.QInputDialog.getText(self,'Rename experiment','old name=%s' % cname)
			if ok:
				self.removeexp(cname)
				cexp.studyname=val
				self.addexp(cexp)


	def menuRemove(self):
		if len(self.bMainList.selectedItems())>1:
			if QtGui.QMessageBox.warning(self,"Remove samples?","Remove %d samples?" % len(self.bMainList.selectedItems()),QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)==QtGui.QMessageBox.No:
				return
		for currentItemName in self.bMainList.selectedItems():
			currentItemName=str(currentItemName.text())
#		currentItemName=str(self.bMainList.currentItem().text())
			self.removeexp(currentItemName)

	def menuSave(self):
		cname=str(self.bMainList.currentItem().text())
		fname = str(QtGui.QFileDialog.getSaveFileName(self, 'Save experiment as pickle',''))
		fl=open(fname,'w')
		pickle.dump(self.explist[cname],fl,-1)
		fl.close()
		QtGui.QMessageBox.information(self,'Analysis','experiment %s saved as pickle' % cname)
#		picklewrapper.save('test',currentItemName)

	def menuExport(self):
		cname=str(self.bMainList.currentItem().text())
		fname = str(QtGui.QFileDialog.getSaveFileName(self, 'Save experiment as biom',''))
		hs.savebiom(self.explist[cname],fname)
		QtGui.QMessageBox.information(self,'Analysis','experiment %s saved as biom table and mapping file' % cname)
#		cname=str(self.bMainList.currentItem().text())
#		cm='global %s;%s=self.explist[cname]' % (cname,cname)
#		exec(cm)
#		hs.Debug(7,'exported',cname)


	def addexp(self,expdat):
		"""
		add a new experiment to the list
		"""

		# make sure the experiment is not already in the list
		# if so, give a new unique name
		expname=expdat.studyname
		cnum=1
		while expname in self.explist:
			expname=expdat.studyname+str(cnum)
			cnum+=1
		expdat.studyname=expname
		expdname='%s (%s-S, %s-B)' % (expname,hs.nicenum(len(expdat.samples)),hs.nicenum(len(expdat.sids)))
		self.explist[expdname]=expdat
		self.bMainList.addItem(expdname)
		self.bMainList.clearSelection()
		self.bMainList.setCurrentRow(self.bMainList.count()-1)


	def replaceexp(self,expdat):
		"""
		replace an existing experiment with new values
		"""
		expname=expdat.studyname
		self.explist[expname]=expdat
		items=self.bMainList.findItems(expname,QtCore.Qt.MatchExactly)
		for item in items:
			self.bMainList.takeItem(self.bMainList.row(item))
			self.bMainList.addItem(expname)


	def removeexp(self,expname):
		"""
		remove an experiment from the list (and clear)
		"""
		del self.explist[expname]
		items=self.bMainList.findItems(expname,QtCore.Qt.MatchExactly)
		for item in items:
			self.bMainList.takeItem(self.bMainList.row(item))

	def joinfields(self):
		items=self.bMainList.selectedItems()
		if len(items)!=1:
			print("Need 1 item")
			return
		for citem in items:
			cname=str(citem.text())
			cexp=self.explist[cname]
			joinwin = JoinWindow(cexp)
			res=joinwin.exec_()
			if res==QtGui.QDialog.Accepted:
				fieldname1=str(joinwin.cField1.currentText())
				fieldname2=str(joinwin.cField2.currentText())
				newfieldname=str(joinwin.lineEdit.text())
				expdat=hs.joinfields(cexp,fieldname1,fieldname2,newfieldname)
				self.addexp(expdat)


	def load(self):
		loadwin = LoadWindow()
		res=loadwin.exec_()
		if res==QtGui.QDialog.Accepted:
			tablefname=str(loadwin.tLoadTable.text())
			mapfname=str(loadwin.tLoadMap.text())
			expname=str(loadwin.tLoadName.text())
			metabolite=loadwin.cMetabolite.checkState()
			if metabolite:
				tabletype='meta'
			else:
				tabletype='biom'
			expdat=hs.load(tablefname,mapfname,tabletype=tabletype)
			expdat.studyname=expname
			self.addexp(expdat)
			# for biom table show the number of reads`
			if tabletype=='biom':
#				hs.analyzenumreads(expdat)
				pass

	def pickleload(self):
		fname = str(QtGui.QFileDialog.getOpenFileName(self, 'Open pickle file'))
		if fname:
			fl=open(fname,'r')
			expdat=pickle.load(fl)
			if isinstance(expdat,hs.experiment):
				self.addexp(expdat)
			else:
				hs.Debug(9,'Not an experiment pickle!')
				print(type(expdat))


	def plot(self):
		items=self.bMainList.selectedItems()
		for citem in items:
			cname=str(citem.text())
			hs.plotexp(self.explist[cname],usegui=True,sortby=False,seqdb=self.bactdb,cdb=self.cooldb)


	def advplot(self):
		items=self.bMainList.selectedItems()
		if len(items)!=1:
			print("Need 1 item")
			return
		for citem in items:
			cname=str(citem.text())
			cexp=self.explist[cname]
			dwin = AdvPlotWindow(cexp)
			res=dwin.exec_()
			if res==QtGui.QDialog.Accepted:
				minreads=dwin.sMinReads.value()
				cnumeric=dwin.cNumeric.checkState()
				if cnumeric==0:
					numeric=False
				else:
					numeric=True
				cshowline=dwin.cLines.checkState()
				if cshowline==0:
					showline=False
				else:
					showline=True
				csort=dwin.cSort.checkState()
				if csort==0:
					sortfield=False
				else:
					sortfield=str(dwin.cField.currentText())
				hs.plotexp(cexp,sortby=sortfield,numeric=numeric,showline=showline,minreads=minreads,usegui=True,cdb=self.cooldb,seqdb=self.bactdb)


	def filterfasta(self):
		items=self.bMainList.selectedItems()
		if len(items)!=1:
			print("Need 1 item")
			return
		for citem in items:
			cname=str(citem.text())
			cexp=self.explist[cname]
			filterfastawin = FilterFastaWindow()
			res=filterfastawin.exec_()
			if res==QtGui.QDialog.Accepted:
				filename=str(filterfastawin.tFileName.text())
				exclude=filterfastawin.cExclude.checkState()
				pmatch=filterfastawin.cPartialMatch.checkState()
				newname=str(filterfastawin.tNewName.text())
				if not newname:
					newname=cexp.studyname+'-ffa'
				overwrite=filterfastawin.cOverwrite.checkState()
				newexp=hs.filterfasta(cexp,filename,exclude=exclude,subseq=pmatch)
				if overwrite==0:
					newexp.studyname=newname
					self.addexp(newexp)
				else:
					self.replaceexp(newexp)


	def diffexp(self):
		items=self.bMainList.selectedItems()
		if len(items)!=1:
			print("Need 1 item")
			return
		for citem in items:
			cname=str(citem.text())
			cexp=self.explist[cname]
			diffexpwin = DiffExpWindow(cexp)
			res=diffexpwin.exec_()
			if res==QtGui.QDialog.Accepted:
				field=str(diffexpwin.cField.currentText())
				value1=str(diffexpwin.tValue1.text())
				value2=str(diffexpwin.tValue2.text())
				newname=str(diffexpwin.tNewName.text())
				method=str(diffexpwin.cMethod.currentText())
				compareall=diffexpwin.cAll.checkState()
				if compareall:
					value2=False
				if method=='all':
					newexp=hs.getdiffsigall(cexp,field,value1,value2)
				else:
					newexp=hs.getdiffsig(cexp,field,value1,value2,method=method)
				if newexp:
					newexp.studyname=newname+'_'+field
					self.addexp(newexp)


	def classify(self):
		items=self.bMainList.selectedItems()
		if len(items)!=1:
			print("Need 1 item")
			return
		for citem in items:
			cname=str(citem.text())
			cexp=self.explist[cname]
			classifywin = ClassifyWindow(cexp)
			res=classifywin.exec_()
			if res==QtGui.QDialog.Accepted:
				field=str(classifywin.cField.currentText())
				value1=str(classifywin.tValue1.text())
				value2=str(classifywin.tValue2.text())
				method=str(classifywin.cMethod.currentText())
				compareall=classifywin.cAll.checkState()
				if compareall:
					value2=False
				auc=hs.BaysZeroClassifyTest(cexp,field,value1,value2,numiter=5)
				hs.Debug(6,"AUC is %f" % auc)

	def filtersamples(self):
		items=self.bMainList.selectedItems()
		if len(items)!=1:
			print("Need 1 item")
			return
		for citem in items:
			cname=str(citem.text())
			cexp=self.explist[cname]
			filtersampleswin = FilterSamplesWindow(cexp)
			res=filtersampleswin.exec_()
			if res==QtGui.QDialog.Accepted:
				field=str(filtersampleswin.cField.currentText())
				value=str(filtersampleswin.tValue.text())
				newname=str(filtersampleswin.tNewName.text())
				overwrite=filtersampleswin.cOverwrite.checkState()
				exclude=filtersampleswin.cExclude.checkState()
				exact=filtersampleswin.cExact.checkState()
				newexp=hs.filtersamples(cexp,field,value,exclude=exclude,exact=exact)
				if overwrite==0:
					newexp.studyname=newname
					self.addexp(newexp)
				else:
					self.replaceexp(newexp)

	def sortabundance(self):
		items=self.bMainList.selectedItems()
		if len(items)!=1:
			print("Need 1 item")
			return
		for citem in items:
			cname=str(citem.text())
			cexp=self.explist[cname]
			filtersampleswin = FilterSamplesWindow(cexp)
			res=filtersampleswin.exec_()
			if res==QtGui.QDialog.Accepted:
				field=str(filtersampleswin.cField.currentText())
				value=str(filtersampleswin.tValue.text())
				newname=str(filtersampleswin.tNewName.text())
				overwrite=filtersampleswin.cOverwrite.checkState()
				exact=filtersampleswin.cExact.checkState()
				newexp=hs.sortbyfreq(cexp,field,value,exact=exact)
				if overwrite==0:
					newexp.studyname=newname
					self.addexp(newexp)
				else:
					self.replaceexp(newexp)


	def sortsamples(self):
		items=self.bMainList.selectedItems()
		if len(items)!=1:
			print("Need 1 item")
			return
		for citem in items:
			cname=str(citem.text())
			cexp=self.explist[cname]
			sortsampleswin = SortSamplesWindow(cexp)
			res=sortsampleswin.exec_()
			if res==QtGui.QDialog.Accepted:
				field=str(sortsampleswin.cField.currentText())
				newname=str(sortsampleswin.tNewName.text())
				cnumeric=sortsampleswin.cNumeric.checkState()
				if cnumeric==0:
					numeric=False
				else:
					numeric=True
				overwrite=sortsampleswin.cOverwrite.checkState()
				newexp=hs.sortsamples(cexp,field=field,numeric=numeric)
				if overwrite==0:
					newexp.studyname=newname
					self.addexp(newexp)
				else:
					self.replaceexp(newexp)


	def clusterbacteria(self):
		items=self.bMainList.selectedItems()
		if len(items)!=1:
			print("Need 1 item")
			return
		for citem in items:
			cname=str(citem.text())
			cexp=self.explist[cname]
			val,ok=QtGui.QInputDialog.getInt(self,'Cluster Bacteria','Minimal number of reads per bacteria',10,0,10000)
			if ok:
				newexp=hs.clusterbacteria(cexp,minreads=val)
				newexp.studyname=newexp.studyname+'_cb'
				self.addexp(newexp)


	def bicluster(self):
		items=self.bMainList.selectedItems()
		if len(items)!=1:
			print("Need 1 item")
			return
		for citem in items:
			cname=str(citem.text())
			cexp=self.explist[cname]

			biclusterwin = BiClusterWindow(cexp,cdb=self.cooldb,bdb=self.bactdb)
			res=biclusterwin.exec_()
			if res==QtGui.QDialog.Accepted:
				newexp=hs.reorderbacteria(cexp,cexp.bactorder)
				newexp=hs.reorderbacteria(newexp,cexp.samporder)
				newexp.studyname=newexp.studyname+'_bicluster'
				newexp.filters.append("bicluster")
				self.addexp(newexp)

#			newexp=hs.bicluster(cexp,method='binary')
#			newexp.studyname=newexp.studyname+'_bicluster'
#			self.addexp(newexp)

#			if len(newexp.samples)>5 and len(newexp.seqs)>5:
#				pv=hs.testbactenrichment(cexp,newexp.seqs,cdb=False,bdb=False,dbexpres=False,translatestudy=False)

	def enrichment(self):
		items=self.bMainList.selectedItems()
		if len(items)!=1:
			print("Need 1 item")
			return
		for citem in items:
			cname=str(citem.text())
			cexp=self.explist[cname]
			if self.cooldb:
				evals=hs.cooldb.testenrichment(self.cooldb,[],cexp.seqs,freqs=np.log(1+np.mean(cexp.data,axis=1)))
#				sevals=hs.sortenrichment(evals,method='val')
				sevals=hs.sortenrichment(evals,method='single',epsilon=1)
				for cval in sevals[::-1]:
					print('%s - obs:%f catfreq:%f' % (cval['description'],cval['observed'],cval['expected']))

	def joinexps(self):
		items=self.bMainList.selectedItems()
		if len(items)<2:
			print("Need at least 2 experiments!")
			return
		newfieldname,ok=QtGui.QInputDialog.getText(self,'Join experiments','Name of study name field')
		if ok:
			cname=str(items[0].text())
			baseexp=self.explist[cname]
			allexpname=baseexp.studyname
			for citem in items[1:]:
				cname=str(citem.text())
				cexp=self.explist[cname]
				baseexp=hs.joinexperiments(baseexp,cexp,missingval='NA',origfieldname=newfieldname)
				allexpname+='+'+cexp.studyname
			baseexp.studyname=allexpname
			self.addexp(baseexp)



	def filterorigreads(self):
		items=self.bMainList.selectedItems()
		if len(items)!=1:
			print("Need 1 item")
			return
		for citem in items:
			cname=str(citem.text())
			cexp=self.explist[cname]
			val,ok=QtGui.QInputDialog.getInt(self,'Filter Original Reads','Minimal number of reads per sample',5000,0,100000)
			if ok:
				newexp=hs.filterorigreads(cexp,minreads=val)
				newexp.studyname=newexp.studyname+'_for'
				self.addexp(newexp)


	def filterminreads(self):
		items=self.bMainList.selectedItems()
		if len(items)!=1:
			print("Need 1 item")
			return
		for citem in items:
			cname=str(citem.text())
			cexp=self.explist[cname]
			val,ok=QtGui.QInputDialog.getInt(self,'Filter min reads','Minimal number of reads per bacteria',10,0,10000)
			if ok:
				newexp=hs.filterminreads(cexp,minreads=val)
				newexp.studyname=newexp.studyname+'_fmr'
				self.addexp(newexp)


	def filterpresence(self):
		items=self.bMainList.selectedItems()
		if len(items)!=1:
			print("Need 1 item")
			return
		for citem in items:
			cname=str(citem.text())
			cexp=self.explist[cname]
			val,ok=QtGui.QInputDialog.getDouble(self,'Filter presence','Minimal fraction of samples per bacteria',value=0.5,min=0,max=1,decimals=2)
			if ok:
				newexp=hs.filterpresence(cexp,frac=val)
				newexp.studyname=newexp.studyname+'_fp'
				self.addexp(newexp)


	def filtermean(self):
		items=self.bMainList.selectedItems()
		if len(items)!=1:
			print("Need 1 item")
			return
		for citem in items:
			cname=str(citem.text())
			cexp=self.explist[cname]
			val,ok=QtGui.QInputDialog.getDouble(self,'Filter Mean','Minimal mean fraction of reads per sample bacteria ',value=0.01,min=0,max=1,decimals=5)
			if ok:
				newexp=hs.filtermean(cexp,meanval=val*10000)
				newexp.studyname=newexp.studyname+'_fmean'
				self.addexp(newexp)

	def subsample(self):
		items=self.bMainList.selectedItems()
		if len(items)!=1:
			print("Need 1 item")
			return
		for citem in items:
			cname=str(citem.text())
			cexp=self.explist[cname]
			hs.analyzenumreads(cexp)
			val,ok=QtGui.QInputDialog.getInt(self,'Subsample','Number of reads per sample',value=10000,min=0)
			if ok:
				newexp=hs.subsample(cexp,numreads=val)
				newexp.studyname=newexp.studyname+'_sub'
				self.addexp(newexp)


	def filtertaxonomy(self):
		items=self.bMainList.selectedItems()
		if len(items)!=1:
			print("Need 1 item")
			return
		for citem in items:
			cname=str(citem.text())
			cexp=self.explist[cname]
			val,ok=QtGui.QInputDialog.getText(self,'Filter taxonomy','Taxonomy to filter')
			if ok:
				newexp=hs.filtertaxonomy(cexp,tax=str(val),exact=False)
				newexp.studyname=newexp.studyname+'_ftax'
				self.addexp(newexp)

	def renormalize(self):
		items=self.bMainList.selectedItems()
		if len(items)!=1:
			print("Need 1 item")
			return
		for citem in items:
			cname=str(citem.text())
			cexp=self.explist[cname]
			newexp=hs.normalizereads(cexp)
			newexp.studyname=newexp.studyname+'_norm'
			self.addexp(newexp)


def main():
	print("starting heatsequer")
	app = QtGui.QApplication(sys.argv)
	print("almost ready")
	window = AppWindow()
	window.show()
	sys.exit(app.exec_())


if __name__ == '__main__':
#	bdb=pickle.load(open('bactdb.pickle'))
#	print(bdb.OntoGraph.keys())
	main()
