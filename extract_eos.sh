#!/bin/bash

# Minimum Z
z_beg=16

# Maximum Z
z_end=60

# Storage for all summary files
mkdir -p summary

# Storages for all profile parameters, all energies, and those of the minima
> eos_pars.dat
> mins.dat
> mins_pars.dat
> eos.dat

# Column headings
echo "rho_B Z N R_WS rho_gas_n rho_liq_n r_n a_n gamma_n rho_gas_p rho_liq_p r_p a_p gamma_p" >> eos_pars.dat
echo "rho_B Z N R_WS rho_gas_n rho_liq_n r_n a_n gamma_n rho_gas_p rho_liq_p r_p a_p gamma_p" >> mins_pars.dat
echo "rho_B Z N tau_n tau_p tau_t field spin-orbit coul_dir coul_exc coul_tot skyrme strut_p pair_n e_pair_bcs_p k_e coul_ee coul_pe mu_n mu_p mu_e mu_t e_t e/a_t p_nucl p_e p_ex p_t" >> eos.dat
echo "rho_B Z N tau_n tau_p tau_t field spin-orbit coul_dir coul_exc coul_tot skyrme strut_p pair_n e_pair_bcs_p k_e coul_ee coul_pe mu_n mu_p mu_e mu_t e_t e/a_t p_nucl p_e p_ex p_t" >> mins.dat

# Loop over directory for every density
for dir in rho_b_*/ ; do
    # Extract density from directory name, as decimal number
    rho_b=$(echo $dir | grep -Eo '0\.[0-9]*')
    
    # Enter directory for a single density
    cd $dir
    
    # Store all profile parameters
    cat cell_shape_rho_b_${rho_b}_${z_beg}_${z_end}.dat | awk '{if(NR>1)print $0}' >> ../eos_pars.dat
    
    # Store all energies
    cat all_energies_rho_b_${rho_b}_${z_beg}_${z_end}.dat | awk '{if(NR>1)print $0}' >> ../eos.dat
    
    # Store copies of summary files in "summary" directory
    cp cell_shape_rho_b_${rho_b}_${z_beg}_${z_end}.dat ../summary
    cp all_energies_rho_b_${rho_b}_${z_beg}_${z_end}.dat ../summary
    
    # Store profile parameters of minimum
        # Get Z value from energy value
    z_min=$(cat all_energies_rho_b_${rho_b}_${z_beg}_${z_end}.dat | awk '{if(NR>1)print $0}' | sort -gk24 | head -1 | awk '{print $2}')
    
        # Get row containing Z value from profile parameters file
    grep " ${z_min} " cell_shape_rho_b_${rho_b}_${z_beg}_${z_end}.dat >> ../mins_pars.dat
    
    # Store all energies of minimum
    cat all_energies_rho_b_${rho_b}_${z_beg}_${z_end}.dat | awk '{if(NR>1)print $0}' | sort -gk24 | head -1 >> ../mins.dat
    
    # Return to main directory
    cd ../
done

# Put in aligned columns
cp eos_pars.dat temp.dat
cat temp.dat | column -t > eos_pars.dat
cp mins_pars.dat temp.dat
cat temp.dat | column -t > mins_pars.dat
cp eos.dat temp.dat
cat temp.dat | column -t > eos.dat
cp mins.dat temp.dat
cat temp.dat | column -t > mins.dat
rm temp.dat

# Add "# " to first line so that plotting ignores it
sed -i '1s/^/# /' eos_pars.dat
sed -i '1s/^/# /' mins_pars.dat
sed -i '1s/^/# /' eos.dat
sed -i '1s/^/# /' mins.dat

# Add 2 spaces to all other lines to keep alignment
sed -i -e '2,$ s/^/  /' eos_pars.dat
sed -i -e '2,$ s/^/  /' mins_pars.dat
sed -i -e '2,$ s/^/  /' eos.dat
sed -i -e '2,$ s/^/  /' mins.dat

# Add long summary files to "summary" directory
cp eos_pars.dat mins_pars.dat eos.dat mins.dat summary