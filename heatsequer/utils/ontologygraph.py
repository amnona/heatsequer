#!/usr/bin/env python


"""
heatsequer
parse ontology as a networkx graph
"""

# amnonscript

from ..utils.oboparse import Parser
import networkx as nx


def ontologytotree(filename):
	"""
	add ontology entries to a tree (based on 'is_a' field)

	input:
	filename : str
		name of the obo ontology file

	output:
	g : networkx graph
	names : list of str
		names of nodes of g
	"""

	g=nx.DiGraph()

	parser=Parser(open(filename))
	for citem in parser:
		cid=citem.tags["id"][0]
		if "is_a" in citem.tags:
			isa=citem.tags["is_a"][0]
			g.add_edge(isa,cid)
	return g


def ontotreetonames(g,termdict):
	"""
	convert an ontology tree with IDs as node ids to a tree with the name as the node ids

	input:
	g : networkx DiGraph
		from ontologytotree()
	termdict : dict of str:str
		a dict mapping IDs to names (such as scdb.ontologyfromid)

	output:
	ng : networkx DiGraph
		with names instead of ids of the nodes
	"""
	ng=nx.relabel_nodes(g,termdict,copy=True)
	return ng


def ontologysubtreeids(g,node):
	"""
	get the ids of a subtree of an ontology graph (from ontologytotree() )

	input:
	g : networkx graph ((from ontologytotree() )
	node : str
		name of the subtree root

	output:
	names : dict of ids
		a dict with ids as keys for the ontology (and True as value)
	"""

	g2=nx.algorithms.traversal.dfs_tree(g,node)
	names={}
	for cid in g2.nodes():
		names[cid]=True
	return names


def getnodeparents(g,node):
	"""
	get a list of all parents of a node

	input:
	g : networkx DiGraph
	node : node in the graph

	output:
	parents : list of nodes
		a list of parents of the node
	"""
	if node not in g:
		return([])
	parents=[]
	tnodes=[node]
	while tnodes:
		cnode=tnodes.pop()
		parents.append(cnode)
		tnodes+=g.predecessors(cnode)
	return parents
