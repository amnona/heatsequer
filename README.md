HeatSequer
----------
look at the secret life of bacteria

HeatSequer is a program for visualization (using heatmaps) and analysis of 16S amplicon data. It is intended for large studies (100s of samples, 1000s of OTUs). While being adapted mostly for deblurred 16S sequences, it can work with closed/open reference OTUs and metabolomics data.
HeatSequer can be used as an interactive program or as a set of interactive python commands for data exploration (also in Jupyter notebook).


Installation
------------
The easiest way to install is:

Install miniconda (http://conda.pydata.org/miniconda.html) or anaconda ( https://www.continuum.io/downloads ). HeatSequer is tested to work with Python 2 and Python 3

Create a virtual environment for heatsequer:

in the command prompt, run:

```
conda create --name heatsequer scipy numpy networkx pyqt qt scikit-learn matplotlib hdf5 h5py
```

Activate the environment:

on mac/unix:

```
source activate heatsequer
```

on windows:

```
activate heatsequer
```

Install additional required packages:

```
pip install biom-format

pip install h5py

pip install scikit-bio
```

Download the heatsequer files from github (https://github.com/amnona/heatsequer):

cd to where you want the new heatsequer directory to be created

```
git clone git@github.com:amnona/heatsequer.git
```

To run the program:
-------------------
from the direcroty where heatsequer is downloaded to type:

on windows:

```
activate heatsequer
```

on mac/unix:

```
source activate heatsequer
```

and then (on win/mac/unix):

```
python hsgui.py
```


Basic Usage:
------------
1. Loading a biom table:

1a. select "Load New"

1b. select the name of the biom table (you can use the Browse button to select). The biom table can be any format (text/json/hdf5) and can originate from deblurring, closed reference or open reference OTU picking.

1c. select the name of the mapping file associated with the biom table

1d. you can change the name that will be assigned to this table in "Experiment Name". The default is the biom table file name

1e. make sure Metabolite checkbox is not selected.

1f. select "Load"


2. Loading a metabolite csv table (for metabolomics bucket table):

2a. select "Load New"

2b. select the name of the metabolite csv bucket table (you can use the Browse button to select).

2c. select the name of the mapping file associated with the metabolite table

2d. you can change the name that will be assigned to this table in "Experiment Name". The default is the biom table file name

2e. make sure Metabolite checkbox is selected (checked).

2f. select "Load"


3. Preprocessing the data:

Once a table is loaded, it is displayed in the "Experiments" list, together with the number of samples (XXX-S) and OTUs (XXX-B). Right clicking on the table enables deletion and saving.

Before plotting the data, it is recommended to do bit preprocessing (removing uninteresting OTUs, ordering the samples and clustering the OTUs):

3a. removing unintersting OTUs:

under the "Bacteria" tab, select "Filter MinReads". Typically we can use 10 as the total minimal number of reads (total out of 10K per sample) so we throw away all OTUs that appear less than 10 reads in all samples combined.

3b. clustering the OTUs:

under the "Bacteria" tab, select "Cluster Bacteria", then typically use 10 as the minimal number of reads (similar to 3a). This is in order to make clustering faster (usually works fast for <5K OTUs).

4. Ordering the samples:

4a. Under the "Samples" tab, select Sort Samples. Then select the field to sort by. If values in the field are numeric, check the "Numeric" checkbox


5. Analyzing differential expression / correlation

In order to identify otus that differentiate between 2 sample groups in the experiment:

5.1 From the Analysis tab, select Diff. Expr.

5.2 Select the field by which to separate your 2 groups

5.3 Select the value for each group (for group 2, you can instead select the "All" checkbox, which means all samples except group 1)

5.4 Select the difference detection method (all tests are permuation based and FDR corrected):

5.4a. ranksum - find otus with different in the mean of the rank the 2 groups (less sensitive to outliers)

5.4b. binary - find otus with different absence/presence between the 2 groups (binary)

5.4c. mean - compare the mean presence of the otu between the 2 groups

5.4d. freqpres - compare the mean presence only in samples where the otu is present>0 (ignore the 0 freq. samples in each group)

5.4e. All - select all OTUs that are different between the 2 groups in at least one of the criteria



6. Plotting the data:

6a. Select the table you want to plot (from the "Experiment" list)

6b. Select "Advanced Plot"

6c. Select the field by which to order the samples. If it has numeric values, check the "Numeric" checkbox (so "10" will be after "2")

6d. If the number of unique values of the field is not too big, select the "Draw Lines" checkbox.

6e. To plot the heatmap, select "Plot"


7. The Heatmap Display:

When plotting a heatmap, 2 windows will open: the heatmap window and the info window. The info window shows information about the selected sample and OTU

The heatmap window shows a log scaled heatmap, with each sample in a columns and each OTU in a row

7a. Keyboard shortcuts:

NOTE: all keyboard shortcuts work only when the mouse is over the heatmap itself.

q - zoom in (OTUs - Y axis)

a - zoom out (OTUs - Y axis)

Q (shift+q) - zoon in (samples - X axis)

A (shift+a) - zoon out (samples - X axis)

f - toggle full screen

h - show OTU labels (use only in zooms when less than 200 otus are visible, otherwise slow)

up and down arrows - scroll up and down the visible OTUs

left and right arrows - scroll left and right the visible samples

, and . - move the selected otu cursor 1 up or down

left and right arrows - move the selected sample cursor 1 left or right

multiple otu selection can be done with shift/command and clicking


8. The information window:

Whenever you do a plot, an additional window opens with info about the selected sample/otu.

The display can be divided into 4 regions:

8.1. top part - selected sample information

8.1a. You can select the sample field to show for the selected sample from the drop down box

8.2. A few buttons:

8.2a. Get Sequence - get the sequence for the selected otu (for blast, etc)

8.2b. ExpInfo - information about the experiment being viewed (command history - how we got this experiment)

8.2c. Sample info - information about the selected sample (all mapping file metadata fields)

8.3. List of manual curated database info about current selected otu

8.4. A few buttons for the selected otus:

8.4a. Export - ignore

8.4b. Save - save the selected otus into a fasta file with all their sequences (For later filtering, blast, etc)

8.4c. View - see a list of all selected otus

8.4d. Enrich - do an enrichment analysis for the selected bacteria - what manual curation database entries are enriched in the selection compared to other bacteria in the displayed experiment

8.4e. DBSave - add an annotation in the manual curation database for the selected bacteria

8.5. Automatic database info
