#!/bin/bash

# Name of directory containing full run
full_dir=$1

# Change to location of all runs
cd full_inner_crust_runs

# Change to directory of full run
cd $full_dir

# Change to directory containing reruns
cd single_cell_reruns

for single_cell in rho_b*/; do
    # Name of directory for particular density, where replacements will be made
    dens_dir=${single_cell%_Z_[0-9]*}
    
    # Extract Z from directory
    z=$(echo $single_cell | sed -r 's/.*Z_([0-9]+)\//\1/g')
    
    # Extract new data from files in "single_cell_reruns"
    new_cell_shape=$(sed '2q;d' cell_shape_${dens_dir}_${z}_${z}.dat)
    new_energies=$(sed '2q;d' all_energies_${dens_dir}_${z}_${z}.dat)
    
    # Copy files from single_cell directory to original directory
    cp ${single_cell}/* ../${dens_dir}/${single_cell}
    
    # Change to directory where replacements will be made in summary files
    cd ../${dens_dir}
    
    # Insert new data
    sed -i "s/.* $z .*/$new_cell_shape/" cell_shape_${dens_dir}*.dat
    sed -i "s/.* $z .*/$new_energies/" all_energies_${dens_dir}*.dat
    
    # Change back to directory containing reruns
    cd ../single_cell_reruns
done