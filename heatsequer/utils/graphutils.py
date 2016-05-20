#!/usr/bin/env python


"""
heatsequer graph utils module
"""

# amnonscript

__version__ = "0.4"

import heatsequer as hs


def simplifygraph(g,keepnodes=''):
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
	hs.Debug(6,'nodes in graph %d. simplifying' % len(nodelist))
	for cnode in nodelist:
		if cnode in deletenodes:
			continue
		suc=get_nodes(g,cnode,'out')
		pred=get_nodes(g,cnode,'in')
		if len(suc)!=1:
			continue
		if len(pred)!=1:
			continue
		if keepnodes:
			if g.node[suc[0]][keepnodes]==1:
				continue
		g.remove_node(cnode)
		g.add_edge(pred[0],suc[0])
		deletenodes.append(cnode)
		print('deleted')

	nodelist=g.nodes()
	hs.Debug(6,'deleted %d nodes. nodes in graph %d' % (len(deletenodes),len(g.nodes())))
	return deletenodes


def get_nodes(g,node,direction='in'):
	"""
	get a list of nodes that have edges into node

	input:
	g : networkx DiGraph
	node : graph node
	direction : str
		'in' or 'out' for edges going into or out of node respectively

	output:
	nodelist: list of networkx nodes
	"""
	nodelist=[]
	if direction=='in':
		edges=g.in_edges(node)
		for cedge in edges:
			nodelist.append(cedge[0])
	elif direction=='out':
		edges=g.out_edges(node)
		for cedge in edges:
			nodelist.append(cedge[1])
	else:
		print("direction %s not supported" % direction)
	return nodelist
