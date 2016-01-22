#!/usr/bin/env python


"""
heatsequer experiment input/output functions
"""

# amnonscript

__version__ = "0.9"

import heatsequer as hs

import numpy as np
import csv
import biom
import os
from pdb import set_trace as XXX


def load(tablename, mapname='map.txt', taxfile='', nameisseq=True,studyname=False,tabletype='biom',normalize=True,addsname=''):
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

	params=locals()

	# load the table
	if tabletype=='biom':
		hs.Debug(6,'Loading biom table')
		table=biom.load_table(tablename)
	elif tabletype=='meta':
		hs.Debug(6,'Loading metabolite table')
		table=loadmetabuckettable(tablename)
	else:
		hs.Debug(9,'Table type %s not supported' % tabletype)
		return False

	# if need to add constant string to sample names in table
	if addsname!='':
		idtable={}
		ids=table.ids(axis='sample')
		for cid in ids:
			idtable[cid]=addsname+cid
		table=table.update_ids(idtable,axis='sample')


	smap = {}
	mapsamples = []
	if mapname:
		# if mapping file supplied, load it
		mapsamples,smap,fields=loadmap(mapname)
	else:
		# no mapping file, so just create the #SampleID field
		hs.Debug(6,'No mapping file supplied - using just sample names')
		tablesamples = table.ids(axis='sample')
		for cid in tablesamples:
			smap[cid]={'#SampleID':cid}
			mapsamples.append(cid)
		fields=['#SampleID']

	# remove table samples not in mapping file
	tablesamples = table.ids(axis='sample')
	hs.Debug(6,'number of samples in table is %d' % len(tablesamples))
	removelist=[]
	for cid in tablesamples:
		if cid not in mapsamples:
			removelist.append(cid)
			hs.Debug(6,'Table sample %s not found in mapping file' % cid)
	hs.Debug(6,'removing %s samples' % len(removelist))
	table=table.filter(removelist,axis='sample',invert=True)

	tablesamples = table.ids(axis='sample')
	hs.Debug(6,'deleted. number of samples in table is now %d' % len(tablesamples))

	# remove samples not in table from mapping file
	removemap=[]
	for idx,cmap in enumerate(mapsamples):
		if cmap not in tablesamples:
			removemap.append(idx)
			try:
				del smap[cmap]
			except:
				hs.Debug(8,'Duplicate SampleID %s in mapping file' % cmap)
	if len(removemap)>0:
		hs.Debug(7,'removing %d samples from mapping file' % len(removemap))
		mapsamples=hs.delete(mapsamples,removemap)
	hs.Debug(7,'number of samples in mapping file is now %d' % len(mapsamples))


	# get info about the sequences
	tableseqs = table.ids(axis='observation')
	sids = []
	tax = []
	osnames=[]
	for cid in tableseqs:
		# get the original sample name
		osnames.append(cid)
		# get the sid (hash )
		if nameisseq:
			sids.append(hs.hashseq(cid))
		else:
			sids.append(cid)
		# get the taxonomy string
		ctax=gettaxfromtable(table,cid)
		tax.append(ctax)

	if not studyname:
		studyname=os.path.basename(tablename)

	exp=hs.Experiment()
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
		exp=hs.reordersamples(exp,okreads[1])
		colsum=np.sum(exp.data,axis=0,keepdims=True)
	if tabletype=='meta':
		normalize=False

	if normalize:
		exp.data=10000*exp.data/colsum
	else:
		exp.data=10000*exp.data/np.mean(colsum)

	exp.uniqueid=exp.getexperimentid()
	exp=hs.sortbacteria(exp,logit=False)
	hs.addcommand(exp,"load",params=params)
	exp.filters.append('loaded table=%s, map=%s' % (tablename,mapname))
	return(exp)


def loadmetabuckettable(filename):
	'''
	load a metabolomics csv bucket table and convert to a biom table in memory
	input:
	filename : str
		the csv bucket table file name
	output:
	table: biom.Table
		a biom table initialized by the bucket table
	'''

	hs.Debug(6,'Loading metabolite table')
	# load the metabolite table and turn it into a biom table
	fl=open(filename,'rU')
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
			hs.Debug(7,'Metabolite %s (col %d) not numeric!' % (cmet,idx))
	# load sample names
	sampnames=[]
	for cline in fl:
		cline=cline.rstrip('\n')
		cdat=cline.split(',')
		sampnames.append(cdat[0])
	fl.close()
	# now load the table data:
	dat=np.genfromtxt(filename,skip_header=1,usecols=usepos,delimiter=',')
	dat=np.transpose(dat)
	# and create the biom table:
	table=biom.table.Table(dat,metabolites,sampnames)
	# and add metabolite name as taxonomy:
	taxdict={}
	for cmet in metabolites:
		taxdict[cmet]={'taxonomy': cmet}
	table.add_metadata(taxdict,axis='observation')
	return table


def gettaxfromtable(table,seq):
	"""
	get the taxonomy for a given sequence in the biom table
	and convert it to a nicer string
	if no taxonomy - use "unknown"

	input:
	table : biom.table
		the biom table containing the sequence and (possibly) the taxonomy
	seq : str
		the sequence to get taxonomy for

	output:
	tax : str
		the nice taxonomy string
	"""

	tax='unknown'
	md=table.metadata(seq,axis='observation')
	if md:
		if 'taxonomy' in md:
			tax = md['taxonomy']
			if not isinstance(tax,str):
				newtax=''
				for x in tax:
					if len(x)>2:
						if x[2]=='_':
							newtax+=x[3:]
						else:
							newtax+=x
					else:
						newtax+=x
#					tax=[x[3:] if x[2]=='_' else x for x in tax]
#					tax = ';'.join(tax)
				tax=newtax
	return tax


def loadrdptax(expdat,rdpfilename,fastaname=False,threshold=60):
	"""
	load rdp taxonomy (the output of download allrank in the rdp classifier website) and add to biom table
	input:
	expdat - the biom table for which the taxonomy was assigned (sequenced were `d)
	rdpfilename - name of the saved allrank rdp assignment
	fastaname - name of fasta file used for rdp assignment (if it was not from saveseqsforrdp) or False if sequences are in the header of the fasta
	threshold - the assignemt probability threshold under which to not include the assignment (for each level)
	"""
	params=locals()

	if fastaname:
		seqs,headers=hs.readfastaseqs(fastaname)
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
				hs.Debug(6,'sequence %s not found in fasta file' % cseq)
		if not cseq in expdat.seqdict:
			hs.Debug(6,'sequence %s not found in experiment' % cseq)
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
	hs.addcommand(expdat,"loadrdptax",params=params,replaceparams={'expdat':expdat})
	return(expdat)


def saveseqsforrdp(expdat,outfilename):
	"""
	save sequences of an experiment into a fasta file
	with the header identical to the sequence (for easy reloading of rdp taxonomy)
	input:
	expdat : Experiment
	outfilename : string
		name of the output fasta file
	"""
	fl=open(outfilename,'w')
	for cseq in expdat.seqs:
		fl.write('>'+cseq+'\n')
		fl.write(cseq+'\n')
	fl.close()
	hs.Debug(6,'wrote %d sequences to file %s' % (len(expdat.seqs),outfilename))


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
	hs.saveseqsfasta(expdat,expdat.seqs,filename)


def savebiom(expdat,filename):
	"""
	save experiment to text biom table and mapping file
	"""
	savemap(expdat,filename+'.map.txt')
	tablefilename=filename+'.table.txt'
	tf=open(tablefilename,'w')
	tf.write('# Saved biom table from python analysis\n')
	tf.write('#OTUID')
	for csamp in expdat.samples:
		tf.write('\t%s' % csamp)
	tf.write('\ttaxonomy\n')
	for idxseq,cseq in enumerate(expdat.seqs):
		tf.write('%s' % cseq)
		for idx,csamp in enumerate(expdat.samples):
			tf.write('\t%d' % expdat.data[idxseq,idx])
		tf.write('\t%s\n' % expdat.tax[idxseq])
	tf.close()
	hs.Debug(6,'Saved experiment to biom table %s and mapping %s' % (tablefilename,mapfilename))


def savemap(expdat,filename):
	hs.Debug(1,"Saving mapping file %s" % filename)
	mf=open(filename,'w')
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
			mf.write(str(expdat.smap[csamp][cfield]))
		mf.write('\n')
	mf.close()


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

	ldat=hs.copyexp(expdat.data)
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



def loadexptree(expdat,treefilename):
	"""
	EXPERIMENTAL
	load a tree file associated with an experiment
	input:
	expdat
	treefilename - the name of the newick tree file (from make_phylogeny.py). note that need to use the sequences as the fasta sequence ids (use -s in the CreateTable)

	output:
	expdat - with a new field - tree
	"""
	params=locals()

	import skbio.tree

	tree=skbio.tree.TreeNode.read(treefilename)
	hs.Debug(4,'Loaded tree')
	expdat.tree=tree
	expdat.filters.append("Load tree %s" % treefilename)
	hs.addcommand(expdat,"loadexptree",params=params,replaceparams={'expdat':expdat})
	return expdat


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



def savecommands(expdat,filename):
	"""
	save the commands used to generate an experiment to a python text file

	input:
	expdat : Experiment
	filename : string
		the filename to save commands to
	"""

	hs.savelisttofile(expdat.commands,filename,delimiter='\n')
	hs.Debug(6,"%d Commands saved to file %s" % (len(expdat.commands),filename))


def saveseqstolinefile(seqs,filename):
	"""
	save experiment sequeces to a file with 1 sequence per line
	input:
	seqs : list of strings
		The sequences to save
	filename : string
		File to save to
	"""

	fl=open(filename,'w')
	for cseq in seqs:
		fl.write('%s\n' % cseq)
	fl.close()


def loadmap(mapfilename,sampidfield='#SampleID'):
	"""
	load a tab separated mapping file and store the values and samples (from #SampleID field)
	input:
	mapfilename : string
		name of the mapping file
	sampidfield : string
		name of the field containing the sample id (default is #SampleID)

	output:
	mapsamples : list of strings
		the sampleids (in order) - the values in sampidfield
	smap : dict of dict of strings
		indexed by sampleid, then by field name (i.e. 'ENV_MATTER' etc.), values are the actual value of the field for the sample
	fields : list of strings
		names of the mapping file fields (i.e 'ENV_MATTER' etc)
	"""

	smap = {}
	mapsamples = []
	hs.Debug(6,'Loading mapping file %s' % mapfilename)
	mapf = open(mapfilename, 'rU')
	reader = csv.DictReader(mapf, delimiter='\t')
	fields = reader.fieldnames
	for cline in reader:
		cid = cline[sampidfield]
		smap[cid] = cline
		mapsamples.append(cid)
	mapf.close()
	hs.Debug(6,'number of samples in map is %d' % len(mapsamples))
	return mapsamples,smap,fields



def createbiomtablefromexp(expdat):
	"""
	Create a biom table from an experiment
	input:
	expdat : Experiment

	output:
	table - the biom table (with taxonomy)
	"""

	# init the table
	table=biom.table.Table(expdat.data,expdat.seqs,expdat.samples)
	# and add metabolite name as taxonomy:
	taxdict={}
	for idx,cseq in enumerate(expdat.seqs):
		taxdict[cseq]={'taxonomy': expdat.tax[idx]}
	table.add_metadata(taxdict,axis='observation')
	return table


def savetobiom(expdat,filename,format='hdf5'):
	"""
	Save an experiment to a biom table
	input:
	expdat : Experiment
	filename : string
		Name of the file to save to
	format : string
		Format of the file ('hdf5','json','txt')
	"""
	savemap(expdat,filename+'.map.txt')
	hs.Debug(1,'Saving biom table %s' % filename)
	tab=createbiomtablefromexp(expdat)
	if format=='hdf5':
		with biom.util.biom_open(filename, 'w') as f:
			tab.to_hdf5(f, "heatsequer")
	elif format=='json':
		with open(filename,'w') as f:
			tab.to_json("heatsequer",f)
	elif format=='txt':
		s=tab.to_tsv()
		with open(filename,'w') as f:
			f.write(s)
	else:
		hs.Debug(9,'file format not supported')
		return
	hs.Debug(6,'table saved to file %s' % filename)
	return
