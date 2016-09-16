#!/bin/bash
# run heatsequer

# change directory to where the script is located
dir=${0%/*}
if [ "$dir" = "$0" ]; then
  dir="."
fi
cd "$dir"

# run the conda env
source activate heatsequer

# and run heatsequer
python hsgui.py

