import heatsequer as hs
import imp

def reloadme():
	'''
	reload heatsequer and all of it's submodules
	'''

	imp.reload(hs.experiment.filtering)
	imp.reload(hs.experiment.expclass)
	imp.reload(hs.experiment.io)
	imp.reload(hs.experiment.sorting)
	imp.reload(hs.experiment.normalization)
	imp.reload(hs.experiment)
	imp.reload(hs.utils.amnonutils)
	imp.reload(hs.utils.sequence)
	imp.reload(hs.utils)
#	reload(hs.plots.plotwingui)
	imp.reload(hs.plots.plotwin)
	imp.reload(hs.plots.otherplot)
	imp.reload(hs.plots)
	imp.reload(hs.analysis.metrics)
	imp.reload(hs.analysis.mislabels)
	imp.reload(hs.analysis.analyse)
	imp.reload(hs.analysis)
	imp.reload(hs.bactdb.bdb)
	imp.reload(hs.bactdb)
	imp.reload(hs.cooldb.cdb)
	imp.reload(hs.cooldb)
	imp.reload(hs)
