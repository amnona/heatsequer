from __future__ import absolute_import

import os

from heatsequer.experiment import *
from heatsequer.plots import *
from heatsequer.utils import *
from heatsequer.reloadme import *
from heatsequer.analysis import *
import heatsequer.cooldb
import heatsequer.bactdb

heatsequerdir=os.path.dirname(os.path.abspath(os.path.dirname(__file__)))
cdb=False
try:
	Debug(6,'loading cooldb')
	cdb=cooldb.loaddb(os.path.join(heatsequerdir,'db/coolseqs.txt'))
	Debug(6,'cooldb loaded')
except:
	Debug(9,'cooldb not found')


#__all__ = ['experiment','heatmap','utils']
