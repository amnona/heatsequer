#!/usr/bin/env python

"""
amnonscript
analysisgui.py
a gui interface for microbiome experiment analysis
"""

__version__ = "0.1"

import amnonutils as au
import analysis
import argparse
import sys
import numpy as np
import easygui
import Tkinter as tk
from matplotlib.pyplot import *
from os import listdir
from os.path import isfile,join,basename
# for debugging - use XXX()
from pdb import set_trace as XXX


class ag:
	lastdir=''

	def bpload(self):
		print("lala")
		wload=tk.Tk()
		btable=tk.Button(wload,text='Table',command=self.bptable)
		btable.pack()

	def bptable(self):
		self.tablefile=easygui.fileopenbox()

	def analysisgui(self):
		root=tk.Tk()
		sexplist=tk.Scrollbar(root)
		sexplist.pack(side='right',fill='y')
		lexplist=tk.Listbox(root,yscrollcommand=sexplist.set)
		sexplist.config(command=lexplist.yview)
		sexplist.pack(side='left', fill='both',expand=True)
		bload=tk.Button(root,text='Load',command=self.bpload)
		bload.pack(side='right')
		tk.mainloop()

def startgui():
	tg=ag()
	tg.analysisgui()

def main(argv):
	parser=argparse.ArgumentParser(description='GUI for microbiome analysis version '+__version__)
#	parser.add_argument('-l','--lab',help='lab mapping file name')

	args=parser.parse_args(argv)
	startgui()

if __name__ == "__main__":
	main(sys.argv[1:])
