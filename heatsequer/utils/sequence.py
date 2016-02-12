#!/usr/bin/env python


"""
heatsequer sequence module
for analysis of sequence data
"""

# amnonscript

__version__ = "0.9"

import heatsequer as hs

import numpy as np
from pdb import set_trace as XXX


def SeqToArray(seq):
    """ convert a string sequence to a numpy array"""
    seqa=np.zeros(len(seq),dtype=np.int8)
    for ind,base in enumerate(seq):
        if base=='A':
            seqa[ind]=0
        elif base=='a':
            seqa[ind]=0
        elif base=='C':
            seqa[ind]=1
        elif base=='c':
            seqa[ind]=1
        elif base=='G':
            seqa[ind]=2
        elif base=='g':
            seqa[ind]=2
        elif base=='T':
            seqa[ind]=3
        elif base=='t':
            seqa[ind]=3
        elif base=='-':
            seqa[ind]=4
        else:
            seqa[ind]=5
    return(seqa)


def findsimseqs(expdat,oseqs,maxmm=3):
	"""
	Find all sequences in expdat that are similar up to maxmm to any sequence in oseqs
	Note : only mismatches, no indel

	input:
	expdat : Experiment
		The experiment in which to search for similar sequences
	oseqs : string (ACGR sequence) or list of sequences
		The sequences to look for similar sequences in expdat
	maxmm : int
		The maximal number of mismatches in order to keep a similar sequence

	output:
	simseqs : dict of sequences keyed by original sequence
		for each oseq key, values are all the similar sequences from expdat
	"""

	if not isinstance(oseqs,list):
		oseqs=[oseqs]

	# convert all sequences to numpy arrays for fast comparison
	noseqs=[]
	for coseq in oseqs:
		noseqs.append(SeqToArray(coseq))
	neseqs=[]
	for ceseq in expdat.seqs:
		neseqs.append(SeqToArray(ceseq))
	simseqs={}
	for coseq in oseqs:
		simseqs[coseq]=[]
	numsimseqs=[]

	for idxo,coseq in enumerate(noseqs):
		cnumsim=0
		for idxe,ceseq in enumerate(neseqs):
			hdist=np.count_nonzero(np.not_equal(coseq,ceseq))
			if hdist<=maxmm:
				simseqs[oseqs[idxo]].append([expdat.seqs[idxe],hdist,expdat.data[idxe,0]])
				cnumsim+=1
		numsimseqs.append(cnumsim)
		hs.Debug(6,'sequence num %d has %d similar sequences' % (idxo,cnumsim))

	return simseqs


def seqdist(seq1,seq2,trim=True):
	"""
	find the hamming distance between two sequences
	input:
	seq1,seq2 : sequence (string of ACGT)
		The sequences to compare
	trim : bool
		True to compare until min(len(seq1),len(seq2)), False to compare all (count as mismatch)
	output:
	dist : int
		the hamming distance between the 2 sequences (no indel testing)
	"""
	mlen=min(len(seq1),len(seq2))
	dist=0
	for cpos in range(mlen):
		if seq1[cpos]!=seq2[cpos]:
			dist+=1
	if not trim:
		dist+=max(len(seq1),len(seq2))-mlen
	return dist
