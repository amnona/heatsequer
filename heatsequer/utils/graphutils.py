#!/usr/bin/env python


"""
heatsequer graph utils module
"""

# amnonscript

__version__ = "0.4"

import networkx

def simplifygraph(g):
	"""
	remove all intermediate nodes in a graph (i.e. 1 in and 1 out)
	input:
	g : networkx directed graph
		the directed graph to simplify
	output:
	deletenodes : list of integers
		nodes deleted
	"""
	deletenodes=[]
	nodelist=g.nodes()
	for cnode in nodelist:
		suc=g.successors(cnode)
		pred=g.predecessors(cnode)
		if len(suc)==1:
			if len(pred)==1:
				print("deleting %d" % cnode)
				g.remove_node(cnode)
				g.add_edge(pred[0],suc[0])
				deletenodes.append(cnode)
	return deletenodes
