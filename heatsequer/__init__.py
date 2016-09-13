from __future__ import absolute_import

import os

# we need the QT backend
import matplotlib as mpl
mpl.use('Qt4Agg')


from heatsequer.experiment import *
from heatsequer.plots import *
from heatsequer.utils import *
from heatsequer.reloadme import *
from heatsequer.analysis import *
import heatsequer.cooldb
import heatsequer.bactdb
import heatsequer.supercooldb


heatsequerdir=os.path.dirname(os.path.abspath(os.path.dirname(__file__)))

# start the logger (into log.hs.log)
start_log()

# load the cooldb
cdb=False
try:
	Debug(6,'loading cooldb')
	cdb=cooldb.loaddb(get_data_path('coolseqs.txt','db'))
	Debug(6,'cooldb loaded')
except:
	Debug(9,'cooldb not found')

# load the automatic bactdb
bdb=False
try:
	Debug(6,'loading bactdb')
	bdb=bactdb.dbstart(get_data_path('SRBactDB.db','db'))
	Debug(6,'bactdb loaded')
except:
	Debug(9,'bactdb not found')

scdb=False
try:
	Debug(6,'loading supercooldb')
	scdb=supercooldb.dbstart(get_data_path('supercooldb.db','db'))
	Debug(6,'supercooldb loaded')
except:
	Debug(9,'supercool not found')

#__all__ = ['experiment','heatmap','utils']
