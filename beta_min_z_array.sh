#!/bin/bash

# Baryonic density
rho_b=$1

# Minimum Z
z_min=$2

# Maximum Z
z_max=$3

# Create store file for N, R_WS, and profile parameters, for each Z
n_r_profiles=cell_shape_rho_b_${rho_b}_${z_min}_${z_max}.dat
> $n_r_profiles

# Create store file for all energies, for each Z
all_energies=all_energies_rho_b_${rho_b}_${z_min}_${z_max}.dat
> $all_energies

# Column headings
echo "rho_B Z N R_WS rho_gas_n rho_liq_n r_n a_n gamma_n rho_gas_p rho_liq_p r_p a_p gamma_p" >> $n_r_profiles
echo "rho_B Z N tau_n tau_p tau_t field spin-orbit coul_dir coul_exc coul_tot skyrme strut_p pair_n e_pair_bcs_p k_e coul_ee coul_pe mu_n mu_p mu_e mu_t e_t e/a_t p_nucl p_e p_ex p_t" >> $all_energies

# Loop over values of Z in steps of 2
for z in $(seq $z_min 2 $z_max); do
    # Store all input/output files
    store_dir=rho_b_${rho_b}_Z_${z}
    mkdir -p $store_dir
    
    # Perform minimisation
    python minimise_profiles.py --num_p $z --rho_b $rho_b
    
    # Skip to next Z if minimisation failed
    if grep -q "Maximum number of function evaluations has been exceeded." *.log
    then
        # Keep log file to see what went wrong
        mv *.log $store_dir
        continue
        
    else
        # Make sure etf runs in mode 3
        sed -i '/run_mode = /s/[1-4]/3/g' input.in
        
        # Run etf to get full output with all energies
        ./etf |& tee output.out
        
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
        energies[6]=$(cat output.out | grep "Exchange Coulomb |" | awk '{print $4}')
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
        
        # Cleanup files
        mv densities_*.dat fields_*.dat sp_states_p.dat start_profile.in profile.in output.out *.log $store_dir
        cp input.in $store_dir
        
    fi
done

# Put in nice columns
cp $n_r_profiles temp.dat
cat temp.dat | column -t > $n_r_profiles
cp $all_energies temp.dat
cat temp.dat | column -t > $all_energies
rm temp.dat

# Add "# " to first line so that plotting ignores it
sed -i '1s/^/# /' $n_r_profiles
sed -i '1s/^/# /' $all_energies

# Add 2 spaces to all other lines to keep alignment
sed -i -e '2,$ s/^/  /' $n_r_profiles
sed -i -e '2,$ s/^/  /' $all_energies