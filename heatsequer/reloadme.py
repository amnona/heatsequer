import heatsequer as hs

def reloadme():
	'''
	reload heatsequer and all of it's submodules
	'''

	reload(hs.experiment.filtering)
	reload(hs.experiment.expclass)
	reload(hs.experiment.io)
	reload(hs.experiment.sorting)
	reload(hs.experiment.normalization)
	reload(hs.experiment)
	reload(hs.utils.amnonutils)
	reload(hs.utils.sequence)
	reload(hs.utils)
#	reload(hs.plots.plotwingui)
	reload(hs.plots.plotwin)
	reload(hs.plots.otherplot)
	reload(hs.plots)
	reload(hs.analysis.metrics)
	reload(hs.analysis.mislabels)
	reload(hs.analysis.analyse)
	reload(hs.analysis)
	reload(hs.bactdb.bdb)
	reload(hs.bactdb)
	reload(hs.cooldb.cdb)
	reload(hs.cooldb)
	reload(hs)
