#!/bin/bash

# Minimum Z
z_min=20

# Maximum Z
z_max=60

# Density string
dens_string=0.02000

# Location of all Z results
all_Z_path=ale_profiles/SLY4_slices_compare_ETFSI/Slice_rho${dens_string}/

# Summary file
summary=${all_Z_path}WS_fixed_Rho${dens_string}.dat

# Create store file for N, R_WS, and profile parameters, for each Z
n_r_profiles=${all_Z_path}cell_shapes.dat
> $n_r_profiles

# Create store file for E/A for each Z
all_energies=${all_Z_path}etf.dat
> $all_energies

# Make sure etf runs in mode 3
sed -i '/run_mode = /s/[1-3]/3/g' input.in

# Loop over values of Z in steps of 2
for z in $(seq $z_min 2 $z_max); do
    # Get density from summary file
    rho_b=$(cat ${summary} | grep " ${z} " | awk '{print $1}')
    
    # Get box size from summary file
    box=$(cat ${summary} | grep " ${z} " | awk '{print $4}')
    
    # Path for individual Z
    Z_path=${all_Z_path}Z${z}/
    
    # Run fit code
    python new_densities_fit.py --num_p $z --box $box --path $Z_path
    
    # Run "etf"
    make run_etf
    
    # Extract number of neutrons and cell radius
    num_n=$(cat output.out | grep "NUM_N=" | awk '{print $2}')
    r_ws=$(cat output.out | grep "R_WS=" | awk '{print $2}')
    
    # Extract profile parameters
    n_profile=$(cat output.out | grep "N_PROFILE=" | awk -F "=|," '{print $2 $3 $4 $5 $6}')
    p_profile=$(cat output.out | grep "P_PROFILE=" | awk -F "=|," '{print $2 $3 $4 $5 $6}')
    
    # Extract all ETF energy components
    energies[0]=$(cat output.out | grep "Neutron kinetic" | awk '{print $4}')
    energies[1]=$(cat output.out | grep "Proton kinetic" | awk '{print $4}')
    energies[2]=$(cat output.out | grep "Total kinetic" | awk '{print $4}')
    energies[3]=$(cat output.out | grep "Field" | awk '{print $3}')
    energies[4]=$(cat output.out | grep "Spin-orbit" | awk '{print $3}')
    energies[5]=$(cat output.out | grep "Direct Coulomb" | awk '{print $4}')
    energies[6]=$(cat output.out | grep "Exchange Coulomb" | awk '{print $4}')
    energies[7]=$(cat output.out | grep "Total Coulomb" | awk '{print $4}')
    energies[8]=$(cat output.out | grep "Total Skyrme" | awk '{print $4}')
    energies[9]=$(cat output.out | grep "Proton shell-correction" | awk '{print $4}')
    energies[10]=$(cat output.out | grep "Neutron pairing-effect" | awk '{print $4}')
    energies[11]=$(cat output.out | grep "Proton pairing energy" | awk '{print $5}')
    energies[12]=$(cat output.out | grep "Electron kinetic" | awk '{print $4}')
    energies[13]=$(cat output.out | grep "Electron-electron Coulomb" | awk '{print $4}')
    energies[14]=$(cat output.out | grep "Proton-electron Coulomb" | awk '{print $4}')
    energies[15]=$(cat output.out | grep "Neutron chemical potential" | awk '{print $5}')
    energies[16]=$(cat output.out | grep "Proton chemical potential" | awk '{print $5}')
    energies[17]=$(cat output.out | grep "Electron chemical potential" | awk '{print $5}')
    energies[18]=$(cat output.out | grep "Overall chemical potential" | awk '{print $5}')
    energies[19]=$(cat output.out | grep "Corrected total energy" | awk '{print $5}')
    energies[20]=$(cat output.out | grep "Total E/A" | awk '{print $4}')
    energies[21]=$(cat output.out | grep "Nuclear pressure" | awk '{print $4}')
    energies[22]=$(cat output.out | grep "Electron pressure" | awk '{print $4}')
    energies[23]=$(cat output.out | grep "Exchange Coulomb pressure" | awk '{print $5}')
    energies[24]=$(cat output.out | grep "Total pressure" | awk '{print $4}')
    
    # Store rho_b, Z, N, R_WS, and profile parameters
    echo "${rho_b} ${z} ${num_n} ${r_ws} ${n_profile} ${p_profile}" >> $n_r_profiles
    
    # Store rho_b, Z and energies
    echo "${rho_b} ${z} ${num_n} ${energies[@]}" >> $all_energies
    
    # Move "profile.in" and output files to individual Z location
    mv profile.in densities_*.dat fields_*.dat sp_states_p.dat output.out ${Z_path}
done