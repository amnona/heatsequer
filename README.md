HeatSequer
----------
look at the secret life of bacteria

HeatSequer is a program for visualization (using heatmaps) and analysis of 16S amplicon data. It is intended for large studies (100s of samples, 1000s of OTUs). While being adapted mostly for deblurred 16S sequences, it can work with closed/open reference OTUs and metabolomics data.
HeatSequer can be used as an interactive program or as a set of interactive python commands.


Installation
------------
The easiest way to install is:

1. Install miniconda (http://conda.pydata.org/miniconda.html) or anaconda ( https://www.continuum.io/downloads ). HeatSequer is tested to work with Python 2.7

2. Create a virtual environment for heatsequer:

in the command prompt, run:

conda create --name heatsequer scipy numpy networkx pyqt qt scikit-learn matplotlib hdf5 h5py

3. Activate the environment:

activate heatsequer

4. install additional required packages:

pip install biom-format

pip install h5py

5. download the heatsequer files from github (https://github.com/amnona/heatsequer)


To run the program:
-------------------
from the direcroty where heatsequer is downloaded to type:

on windows:

activate heatsequer

on mac/unix:

source activate heatsequer

and then (on win/mac/unix):

python hsgui.py


Basic Usage:
------------
1. Loading a biom table:

1a. select "Load New"

1b. select the name of the biom table (you can use the Browse button to select). The biom table can be any format (text/json/hdf5) and can originate from deblurring, closed reference or open reference OTU picking.

1c. select the name of the mapping file associated with the biom table

1d. you can change the name that will be assigned to this table in "Experiment Name". The default is the biom table file name

1e. make sure Metabolite checkbox is not selected.

1f. select "Load"


2. Loading a metabolite csv table:

1a. select "Load New"

1b. select the name of the metabolite csv bucket table (you can use the Browse button to select).

1c. select the name of the mapping file associated with the metabolite table

1d. you can change the name that will be assigned to this table in "Experiment Name". The default is the biom table file name

1e. make sure Metabolite checkbox is selected (checked).

1f. select "Load"


3. Preprocessing the data:

Once a table is loaded, it is displayed in the "Experiments" list, together with the number of samples (XXX-S) and OTUs (XXX-B). Right clicking on the table enables deletion and exporting.

Before plotting the data, it is recommended to do bit preprocessing (removing uninteresting OTUs, ordering the samples and clustering the OTUs):

3a. removing unintersting OTUs:

under the "Bacteria" tab, select "Filter MinReads". Typically we can use 10 as the total minimal number of reads (total out of 10K per sample) so we throw away all OTUs that appear less than 10 reads in all samples combined.

3b. clustering the OTUs:

under the "Bacteria" tab, select "Cluster Bacteria", then typically use 10 as the minimal number of reads (similar to 3a). This is in order to make clustering faster (usually works fast for <5K OTUs).

4. Ordering the samples:

4a. Under the "Samples" tab, select Sort Samples. Then select the field to sort by. If values in the field are numeric, check the "Numeric" checkbox


4. Plotting the data:

4a. Select the table you want to plot (from the "Experiment" list)

4b. Select "Advanced Plot"

4c. Select the field by which to order the samples. If it has numeric values, check the "Numeric" checkbox (so "10" will be after "2")

4d. If the number of unique values of the field is not too big, select the "Draw Lines" checkbox.

4e. To plot the heatmap, select "Plot"


5. The Heatmap Display:

When plotting a heatmap, 2 windows will open: the heatmap window and the info window. The info window shows information about the selected sample and OTU

The heatmap window shows a log scaled heatmap, with each sample in a columns and each OTU in a row

5a. Keyboard shortcuts:

NOTE: all keyboard shortcuts work only when the mouse is over the heatmap itself.

q - zoom in (OTUs - Y axis)

a - zoom out (OTUs - Y axis)

f - toggle full screen

h - show OTU labels (use only in zooms when less than 200 otus are visible, otherwise slow)

up and down arrows - scroll up and down the visible OTUs

, and . - move the selected otu cursor 1 up or down

left and right arrows - move the selected sample cursor 1 left or right
