#!/usr/bin/env python

"""
amnonscript
full gui for the heatmap analysis
"""

import sys
import amnonutils as au
from PyQt4 import QtGui, uic, QtCore
# import cPickle as pickle
import cPickle as pickle
import matplotlib
matplotlib.use('Qt4Agg')
import analysis
import bactdb
import cooldb
import os.path


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
		val,ok=QtGui.QInputDialog.getItem(self,'Select field value','Field=%s' % cfield,list(set(analysis.getfieldvals(self.cexp,cfield))))
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
		val,ok=QtGui.QInputDialog.getItem(self,'Select field value','Field=%s' % cfield,list(set(analysis.getfieldvals(self.cexp,cfield))))

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
		val,ok=QtGui.QInputDialog.getItem(self,'Select field value','Field=%s' % cfield,list(set(analysis.getfieldvals(self.cexp,cfield))))

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
		val,ok=QtGui.QInputDialog.getItem(self,'Select field value','Field=%s' % cfield,list(set(analysis.getfieldvals(self.cexp,cfield))))
		if ok:
			self.tValue1.setText(val)

	def fieldvalues2(self):
		cfield=str(self.cField.currentText())
		val,ok=QtGui.QInputDialog.getItem(self,'Select field value','Field=%s' % cfield,list(set(analysis.getfieldvals(self.cexp,cfield))))
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
		val,ok=QtGui.QInputDialog.getItem(self,'Select field value','Field=%s' % cfield,list(set(analysis.getfieldvals(self.cexp,cfield))))
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
		self.bFilterOrigReads.clicked.connect(self.filterorigreads)
		self.bMainSortSamples.clicked.connect(self.sortsamples)
		self.bMainClusterBacteria.clicked.connect(self.clusterbacteria)
		self.bMainFilterMinReads.clicked.connect(self.filterminreads)
		self.bMainFilterTaxonomy.clicked.connect(self.filtertaxonomy)
		self.bFilterPresence.clicked.connect(self.filterpresence)
		self.bFilterMean.clicked.connect(self.filtermean)
		self.bJoinFields.clicked.connect(self.joinfields)
		self.bJoinExps.clicked.connect(self.joinexps)
		self.bBicluster.clicked.connect(self.bicluster)
		self.bFilterFasta.clicked.connect(self.filterfasta)
		self.bRenormalize.clicked.connect(self.renormalize)


		# the main list right mouse menu
		self.bMainList.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
		self.bMainList.connect(self.bMainList, QtCore.SIGNAL("customContextMenuRequested(QPoint)"),self.listItemRightClicked)
		self.show()

		self.bactdb=False
		self.cooldb=False
		try:
			au.Debug(6,'Loading sequence database')
			self.bactdb=pickle.load(open('dbtest2.pickle'))
			self.bactdb=bactdb.dbconnect(self.bactdb)
		except:
			au.Debug(9,'BactDB file not found')
		try:
			au.Debug(6,'Loading coolseq database')
			self.cooldb=cooldb.loaddb()
		except:
			au.Debug(9,'CoolDB file not found')

#		self.cooldb=pickle.load(open('cooldb.pickle'))


		# put some experiments:
		au.Debug(6,'Loading sample experiment')
		try:
			expdat=analysis.load('./test_data/bears.clean.new.withtax.biom','./test_data/map.txt')
			self.addexp(expdat)
		except:
			au.Debug(6,'Sample experiment not found. sorry')


	def listItemRightClicked(self, QPos):
		self.listMenu= QtGui.QMenu()
		menuremove = self.listMenu.addAction("Remove Item")
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
		analysis.savebiom(self.explist[cname],fname)
		QtGui.QMessageBox.information(self,'Analysis','experiment %s saved as biom table and mapping file' % cname)
#		cname=str(self.bMainList.currentItem().text())
#		cm='global %s;%s=self.explist[cname]' % (cname,cname)
#		exec(cm)
#		au.Debug(7,'exported',cname)


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
		expdname='%s (%s-S, %s-B)' % (expname,au.nicenum(len(expdat.samples)),au.nicenum(len(expdat.sids)))
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
				expdat=analysis.joinfields(cexp,fieldname1,fieldname2,newfieldname)
				self.addexp(expdat)


	def load(self):
		loadwin = LoadWindow()
		res=loadwin.exec_()
		if res==QtGui.QDialog.Accepted:
			tablefname=str(loadwin.tLoadTable.text())
			mapfname=str(loadwin.tLoadMap.text())
			expname=str(loadwin.tLoadName.text())
			expdat=analysis.load(tablefname,mapfname)
			expdat.studyname=expname
			self.addexp(expdat)
			analysis.analyzenumreads(expdat)

	def pickleload(self):
		fname = str(QtGui.QFileDialog.getOpenFileName(self, 'Open pickle file'))
		if fname:
			fl=open(fname,'r')
			expdat=pickle.load(fl)
			if isinstance(expdat,analysis.experiment):
				self.addexp(expdat)
			else:
				au.Debug(9,'Not an experiment pickle!')
				print(type(expdat))


	def plot(self):
		items=self.bMainList.selectedItems()
		for citem in items:
			cname=str(citem.text())
			analysis.plotexp(self.explist[cname],usegui=True,sortby=False,seqdb=self.bactdb,cdb=self.cooldb)


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
				analysis.plotexp(cexp,sortby=sortfield,numeric=numeric,showline=showline,minreads=minreads,usegui=True,cdb=self.cooldb,seqdb=self.bactdb)


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
				newexp=analysis.filterfasta(cexp,filename,exclude=exclude,subseq=pmatch)
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
					newexp=analysis.getdiffsigall(cexp,field,value1,value2)
				else:
					newexp=analysis.getdiffsig(cexp,field,value1,value2,method=method)
				if newexp:
					newexp.studyname=newname+'_'+field
					self.addexp(newexp)


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
				newexp=analysis.filtersamples(cexp,field,value,exclude=exclude)
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
				newexp=analysis.sortsamples(cexp,field=field,numeric=numeric)
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
				newexp=analysis.clusterbacteria(cexp,minreads=val)
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
			newexp=analysis.bicluster(cexp)
			newexp.studyname=newexp.studyname+'_bicluster'
			self.addexp(newexp)


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
				baseexp=analysis.joinexperiments(baseexp,cexp,missingval='NA',origfieldname=newfieldname)
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
				newexp=analysis.filterorigreads(cexp,minreads=val)
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
				newexp=analysis.filterminreads(cexp,minreads=val)
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
				newexp=analysis.filterpresence(cexp,frac=val)
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
			val,ok=QtGui.QInputDialog.getDouble(self,'Filter Mean','Minimal mean fraction of reads per sample bacteria ',value=0.01,min=0,max=1,decimals=2)
			if ok:
				newexp=analysis.filtermean(cexp,meanval=val*10000)
				newexp.studyname=newexp.studyname+'_fmean'
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
				newexp=analysis.filtertaxonomy(cexp,tax=str(val),exact=False)
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
			newexp=analysis.normalizereads(cexp)
			newexp.studyname=newexp.studyname+'_norm'
			self.addexp(newexp)


def main():
	app = QtGui.QApplication(sys.argv)
	window = AppWindow()
	sys.exit(app.exec_())


if __name__ == '__main__':
#	bdb=pickle.load(open('bactdb.pickle'))
#	print(bdb.OntoGraph.keys())
	main()
