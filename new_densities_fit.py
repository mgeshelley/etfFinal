#!/usr/bin/python3

import numpy as np

from scipy.optimize import curve_fit, minimize, leastsq, fmin_slsqp
from scipy.integrate import quad, simps

import argparse

from minimisation_functions import f_profile, write_profile



def penalty_function_n(p,r,hfb_density,gas=None):
    if gas:
        # Integral of profiles "f_q" as in Onsi et al. (2008)
        I_n = quad(f_profile, 0., box, args=(p[0],p[1],box))[0]
        
        # Neutron liquid density
        rho_n_liq = num_n/(4.*np.pi*I_n) + gas * (1. - box**3/(3.*I_n))
        
        # FD profile density
        fd_density = gas + (rho_n_liq-gas) * f_profile(r,p[0],p[1],box)/r**2
        
    else:
        # Integral of profiles "f_q" as in Onsi et al. (2008)
        I_n = quad(f_profile, 0., box, args=(p[1],p[2],box))[0]
        
        # Neutron gas density
        rho_n_gas = 3./box**3  * (num_n/(4.*np.pi) - I_n*p[0]) / (1. - 3*I_n/box**3)
        
        # FD profile density
        fd_density = rho_n_gas + (p[0]-rho_n_gas) * f_profile(r,p[1],p[2],box)/r**2
    
    return hfb_density - fd_density


def penalty_function_p(p,r,hfb_density):
    # Integral of profiles "f_q" as in Onsi et al. (2008)
    I_p = quad(f_profile, 0., box, args=(p[0],p[1],box))[0]
    
    # Proton liquid density
    rho_p_liq = float(num_p)/(4.*np.pi*I_p)
    
    # FD profile density
    fd_density = rho_p_liq * f_profile(r,p[0],p[1],box)/r**2
    
    return hfb_density - fd_density


## COMMAND LINE ARGUMENTS ##
parser = argparse.ArgumentParser()
parser.add_argument("--num_p", help="number of protons",type=int)
parser.add_argument("--box", help="WS cell radius",type=float)
parser.add_argument("--path", help="location of densities to fit",type=str)
args = parser.parse_args()

num_p = args.num_p
box = args.box
path = args.path

# Fix gas?
fix_gas = True

## HFB DENSITIES ##
# Load data
HFB_profile_n_full = np.loadtxt(path+'DensitiesN.dat')
HFB_profile_p_full = np.loadtxt(path+'DensitiesP.dat')

# Cut densities at radius "box"
HFB_profile_n_cut = HFB_profile_n_full[HFB_profile_n_full[:,0] < box,0:2]
HFB_profile_p_cut = HFB_profile_p_full[HFB_profile_p_full[:,0] < box,0:2]

if fix_gas:
    gas = HFB_profile_n_cut[-1,1]

else:
    gas = None

# Number of neutrons
num_n = simps(4.*np.pi*HFB_profile_n_cut[:,0]**2 * HFB_profile_n_cut[:,1], HFB_profile_n_cut[:,0])

# Initial parameter guesses
if fix_gas:
    x0_n = [5.0,1.0]

else:
    x0_n = [0.1,5.0,1.0]

x0_p = [5.0,1.0]

# Fit parameters
popt_n, pcov_n = leastsq(func=penalty_function_n, x0=x0_n, args=(HFB_profile_n_cut[:,0],HFB_profile_n_cut[:,1],gas), ftol=1e-15, maxfev=10000)
popt_p, pcov_p = leastsq(func=penalty_function_p, x0=x0_p, args=(HFB_profile_p_cut[:,0],HFB_profile_p_cut[:,1]), ftol=1e-15, maxfev=10000)

if fix_gas:
    # Integral of profiles "f_q" as in Onsi et al. (2008)
    I_n = quad(f_profile, 0., box, args=(popt_n[0],popt_n[1],box))[0]

    # Neutron liquid density
    rho_n_liq = num_n/(4.*np.pi*I_n) + gas * (1. - box**3/(3.*I_n))

    # All 5 fitted parameters
    p = np.concatenate(([rho_n_liq],popt_n,popt_p))
    
else:
    # All 5 fitted parameters
    p = np.concatenate((popt_n,popt_p))

# Create "profile.in"
write_profile(p,num_n,num_p,True,box,"profile.in")