#!/bin/bash

# Make directories if they don't already exist
mkdir -p data
mkdir -p plots

mv *.dat data
cp *.in data
mv *.out data
mv *.log data
mv *.npy data
mv *.eps plots
mv *.pdf plots

# gprof ./etf | gprof2dot | dot -Tpng -o profile.png