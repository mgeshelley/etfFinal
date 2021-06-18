#!/usr/bin/python3

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.integrate import quad

import argparse
import subprocess

from minimisation_functions import info_log, run_etf, determine_profiles



## COMMAND LINE ARGUMENTS ##
parser = argparse.ArgumentParser()
parser.add_argument("--num_n", help="number of neutrons", type=float)
parser.add_argument("--num_p", help="number of protons", type=int)
parser.add_argument("--rho_b", help="baryonic density", type=float)
args = parser.parse_args()


## OPTIONS ##
# Print larger amounts on one line
np.set_printoptions(linewidth=np.inf)

# Make overflow warnings into errors, that can be caught during minimisation
np.seterr(over='raise')

# Number of neutrons (if no beta minimisation)
num_n = 78.4
# Get from command line if specified
if args.num_n:
    num_n = args.num_n

# Number of protons
num_p = 50
# Get from command line if specified
if args.num_p:
    num_p = args.num_p

# Baryonic density
rho_b = 0.0025
# Get from command line if specified
if args.rho_b:
    rho_b = args.rho_b

# Gamma = 1?
gamOne = True

# Strutinsky?
strut = True

# Perturbative pairing?
pair_n_perturb = False
pair_p_perturb = False

# Double minimisation to get correct num_n?
determine_n = True

# Use beta-equilibrium to determine num_n, instead of minimising E/A - (Z_e * Q_{n,beta})
beta = False

# Wider range of neutron number
wide_n = True


## CREATE LOG FILE ##
if determine_n:
    log_file = 'rho_b_' + str(rho_b) + '_Z_' + str(num_p) + '.log'
else:
    log_file = 'rho_b_' + str(rho_b) + '_N_' + str(num_n) + '_Z_' + str(num_p) + '.log'

f = open(log_file, 'w')
f.close()


## BOUNDS FOR NEUTRON NUMBER
if wide_n:
    # Wider (needed for LNS, SIV, SKa, SkM*, SQMC, other "soft" PNM EoS)
    bounds = (max([70.*(np.log(3700.*rho_b + 1.)), num_p*2.]), 800.*(np.log(3700.*rho_b + 1.)))
else:
    # Narrower (fine for SLy4, BSk, KDE, SII, other "stiff" PNM EoS)
    bounds = (max([80.*(np.log(3700.*rho_b + 1.)), num_p*2.]), 450.*(np.log(3700.*rho_b + 1.)))


## MINIMISATION ##
if determine_n:
    info_log("Begin determination of num_n:", log_file)
    
    if beta:
        # Minimisation (with beta equilibrium) with respect to neutron number
        n_final, rr = brentq(determine_profiles,
                        bounds[0],
                        bounds[1],
                        #num_p*2.,
                        #num_p*10.,
                        #900.,
                        #1500.,
                        args=(num_p, rho_b, gamOne, strut, pair_n_perturb, pair_p_perturb, beta, log_file),
                        rtol=1.e-6,
                        maxiter=100,
                        full_output=True,
                        disp=True)
        
    else:
        rr = minimize_scalar(determine_profiles,
                        args=(num_p, rho_b, gamOne, strut, pair_n_perturb, pair_p_perturb, beta, log_file),
                        bounds=bounds,
                        #bounds=(num_p*2., num_p*10.),
                        #bounds=(580., 660.),
                        method='bounded',
                        options={'maxiter': 100, 'disp': 2, 'xatol': 1.e-1})
        
        # System exit if profile not converged
        if not rr['success']:
            info_log(rr['message'], log_file)
            sys.exit()
        
        # Optimal value of num_n
        n_final = rr['x']
    
    # Print convergence information
    info_log("Convergence report:", log_file)
    info_log(rr, log_file)
    
    # Final profile with optimum neutron number
    info_log("\nCreate profile at optimum neutron number:", log_file)
    determine_profiles(n_final, num_p, rho_b, gamOne, strut, pair_n_perturb, pair_p_perturb, beta, log_file)
    
    # Write warning message if n_final very close to one of the bounds on neutron number
    if any(abs(bounds-n_final) < 1.):
        info_log("WARNING: n_final very close to one bound", log_file)
        info_log("bounds:  "+str(bounds), log_file)
        info_log("n_final: "+str(n_final)+"\n\n", log_file)
        
        # Store WS cell info
        nWarningsLog = 'nWarnings_rho_b_' + str(rho_b) + '.dat'
        with open(nWarningsLog, 'a') as f2:
            f2.write(' '.join([str(i) for i in [rho_b, num_p, bounds[0], bounds[1], n_final]])+'\n')
    
else:
    info_log("Begin minimisation of energy with respect to profile parameters:", log_file)
    
    # Perform minimisation
    determine_profiles(num_n, num_p, rho_b, gamOne, strut, pair_n_perturb, pair_p_perturb, False, log_file)


### FINAL RESULT ##
## Change "input.in" so that "etf" runs in run_mode 3
#subprocess.run(["sed -i '/run_mode = /s/[1-3]/3/g' input.in"], shell=True)

## Run "etf" to get full output with all energies
#subprocess.run(["make run_etf"], shell=True)