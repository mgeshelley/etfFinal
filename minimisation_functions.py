#!/usr/bin/python3
import sys
import numpy as np

from scipy.optimize import curve_fit, minimize, leastsq, fmin_slsqp, Bounds
from scipy.integrate import quad

import subprocess
from multiprocessing import Pool


# Class for raising when exception raised for bad parameter values
class ParameterValueError(Exception):
    def __init__(self, pars):
        self.pars = pars


# Function to print to screen and log info in "min_log"
def info_log(info,log_file):
    # Print to screen
    print(info)
    
    # Write to file
    with open(log_file, 'a') as min_log:
        print(info, file=min_log)


# Function for integration of "f_q" profiles as in Onsi et al.
def f_profile(r,rq,aq,box):
    return r**2 / (1. + np.exp((r-rq)/aq))
    #print(np.exp(((rq-box)/(r-box))**2 - 1.))
    #if (((rq-box)/(r-box))**2 > 100.):
        #return 0.
    #else:
        #return r**2 / ( 1. + np.exp(((rq-box)/(r-box))**2 - 1.) * np.exp((r-rq)/aq) )


# Create array of parameters "p" with appropriate values added/ calculated
def parameter_array(p,num_n,num_p,gamOne,box):
    # Integrals of profiles "f_q" as in Onsi et al. (2008)
    I_n = quad(f_profile, 0., box, args=(p[1], p[2], box))[0]
    I_p = quad(f_profile, 0., box, args=(p[3], p[4], box))[0]
    
    # Neutron liquid density
    #rho_n_liq = num_n/(4.*np.pi*I_n) + p[0] * (1. - box**3/(3.*I_n))
    #p = np.insert(p,1,rho_n_liq)
    
    # Neutron gas density
    rho_n_gas = 3./box**3  * (num_n/(4.*np.pi) - I_n*p[0]) / (1. - 3*I_n/box**3)
    p = np.insert(p, 0, rho_n_gas)
    
    # Insert "0." in 4th element for proton gas to be fixed at 0
    p = np.insert(p, 4, 0.)
    
    # Proton liquid density
    rho_p_liq = float(num_p) / (4.*np.pi*I_p)
    p = np.insert(p, 5, rho_p_liq)
    
    # Set gamma = 1
    if gamOne:
        p = np.insert(p, 4, 1.)
        p = np.append(p, 1.)
    
    return p


# Function to write parameters to profile with "filename"
def write_profile(p,num_n,num_p,gamOne,box,filename):
    p = parameter_array(p,num_n,num_p,gamOne,box)
    
    # Create and begin writing to filename
    with open(filename, 'w') as profile:
        profile.write("&profiles\n")
        profile.write("    !!! Neutron density profile parameters\n")
        profile.write("    !           rho_gas                 rho_liq                 r                       a                       gamma\n")
        profile.write("    n_profile = ")
        for i in p[0:5]:
            profile.write("{:<24.16e}".format(i))
        profile.write(",\n    \n")
        
        profile.write("    !!! Proton density profile parameters\n")
        profile.write("    !           rho_gas                 rho_liq                 r                       a                       gamma\n")
        profile.write("    p_profile = ")
        for i in p[5:]:
            profile.write("{:<24.16e}".format(i))
        profile.write(",\n    \n")
        
        profile.write("    ! Number of neutrons\n")
        profile.write("    num_n = ")
        profile.write("{:<24.16e}".format(num_n))
        profile.write(",\n    \n")
        
        profile.write("    ! Number of protons\n")
        profile.write("    num_p = ")
        profile.write("{:<d}".format(num_p))
        profile.write(",\n    \n")
        
        profile.write("    ! WS cell radius\n")
        profile.write("    r_ws = ")
        profile.write("{:<24.16e}".format(box))
        profile.write("\n /\n")
     
    return p


# Make profile, run "etf" code, calculate penalty from outputs and biases
def run_etf(p,num_n,num_p,rho_b,gamOne,box,return_delta_mu=False):
    # Create full parameter array
    p = parameter_array(p,num_n,num_p,gamOne,box)
    
    # Sensible limits for parameters
    parameters_sensible = [(p[0] > rho_b/1.e4) & (p[0] < rho_b),    # Neutron gas limits
                           (p[1] > 0.) & (p[1] < 0.16),     # Neutron liquid limits
                           (p[2] > 3.) & (p[2] < box),      # Neutron radius limits
                           (p[3] > 0.) & (p[3] < p[2]),     # Neutron diffuseness limits
                           
                           (p[6] > 0.) & (p[6] < 0.16),     # Proton liquid limits
                           (p[7] > 3.) & (p[7] < box),      # Proton radius limits
                           (p[8] > 0.) & (p[8] < p[7])]     # Proton diffuseness limits
    
    # Raise error if not all parameters sensible
    if not all(parameters_sensible):
        #raise ParameterValueError(p)
        return 1.e9
    
    # Add r_ws to end of parameter array
    p = np.append(p, box)
    
    # Pad each element with spaces up to a length of 24, then join together into one string
    p_string = ''.join([str(par).ljust(24, ' ') for par in p])
    
    # Create command to run "etf" with profile parameters and r_WS as command line inputs
    command = "./etf "+p_string
    
    # Run Fortran program "etf", and pipe all output straight to "out_string"
    #out_string = subprocess.run([command], stdout=sys.stdout, shell=True) # For seeing output of "etf" for debugging
    out_string = subprocess.run([command], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8')
    
    # Split "out_string" into list of individual strings
    split_strings = out_string.split()
    
    # Retrieve outputs of "etf" (last 2 strings), and convert to floats
    etf_output = list(map(float, split_strings[-2:]))
    
    # Extract energy and overall chemical potential
    etf_energy = etf_output[0]
    delta_mu = etf_output[1]
    
    #print(etf_energy, p)
    
    if return_delta_mu:
        return delta_mu
    else:
        return etf_energy


# Determine density profiles for a given N, Z, and baryonic density
def determine_profiles(num_n,num_p,rho_b,gamOne,strut,pair_n_perturb,pair_p_perturb,beta,log_file):
    # WS cell radius
    box = ((num_n+num_p)/(4.*np.pi*rho_b/3))**(1./3)
    
    info_log('{:>13} {:>12.6f}'.format("Neutrons:",num_n),log_file)
    info_log('{:>13} {:>5}'.format("Protons:",num_p),log_file)
    info_log('{:>13} {:>14.8f}'.format("rho_b:",rho_b),log_file)
    info_log('{:>13} {:>12.6f}'.format("Cell radius:",box),log_file)
    
    
    ## MINIMISATION ##
    # Make sure "etf" runs in run_mode 1 for minimisation
    subprocess.run(["sed -i '/run_mode = /s/[1-4]/1/g' input.in"], shell=True)
    
    # Make sure "etf" runs without getting chemical potentials (i.e. getting s.p. states)
    subprocess.run(["sed -i '/calc_chem_pots = /s/true/false/g' input.in"], shell=True)
    
    # Make sure "etf" runs without Strutinsky
    subprocess.run(["sed -i '/strutinsky_on = /s/true/false/g' input.in"], shell=True)
    
    # Make sure "etf" runs without pairing if it will be added perturbatively at end of minimisation
    if pair_n_perturb:
        subprocess.run(["sed -i '/neutron_pairing = /s/[0-2]/0/g' input.in"], shell=True)
    
    if pair_p_perturb:
        subprocess.run(["sed -i '/proton_bcs = /s/true/false/g' input.in"], shell=True)
    
    # Create random x0, and then run minimisation, until no errors raised
    while True:
        try:
            ## INITIAL GUESS AT PARAMETERS ##
            # Initialise list of rules as all False
            rules = [False,False]
            
            # Loop until all rules are met
            while not all(rules):
                # Random values (between reasonable values) of each of the five parameters that are varied
                x0 = np.random.uniform([0.08,4.0,0.6, 4.0,0.4], [0.1,10.0,1.0, 10.0,0.6])
                #x0 = np.random.uniform([0.095,7.0,1.0, 7.0,1.0], [0.1,8.0,2.0, 8.0,2.0])
                
                # Create full array of ten parameters
                p = parameter_array(x0,num_n,num_p,gamOne,box)
                
                # Maximum neutron cluster density (liquid - gas)
                n_cluster = num_n / (4.*np.pi*quad(f_profile, 0.,box,args=(p[2],p[3],box))[0])
                
                # Constraints that parameters must follow
                rules = [   p[0] > 0.,                  # Neutron gas positive
                            p[1]-p[0] <= n_cluster  ]   # Maximum neutron cluster density (Onsi et al. (2008), Eq. 2.10)
            
            # Save start profile as "start_profile.in"
            x0_full = write_profile(x0,num_n,num_p,gamOne,box,"start_profile.in")
            
            info_log('{:>13}'.format("x_0_full:")+' '+' '.join(['{:>11.8f}'.format(i) for i in x0_full]),log_file)
            # Nelder-Mead, no adaptive
                # Energy precise to 1eV, parameters precise to 1e-8 (smallest typical parameter variation is rho_n_gas, ~1e-7)
            result = minimize(run_etf,
                            x0,
                            args=(num_n,num_p,rho_b,gamOne,box,False),
                            method='Nelder-Mead',
                            options={'disp': False, 'maxfev': 3000, 'xatol': 1e-8, 'fatol': 1e-6, 'adaptive': False}
                            #method='BFGS',
                            #options={'disp': True, 'maxiter': 3000},
                            #tol=2e-2
                            )
            
            # Continue if no error raised
            break
            
        # Print error message if overflow encountered, before restarting
        except FloatingPointError as message:
            info_log(
            '{}\n{}'.format("ERROR with x0:"," "*13)+' '+' '.join(['{:>11.8f}'.format(i) for i in x0])+'{} {}\n{}\n'.format(":",message,"Restarting minimisation with new x0"),
            log_file)
            
        # Print error message if bad parameter value(s) reached, before restarting
        except ParameterValueError as error:
            info_log(
            '{}\n{}'.format("Bad value for parameter(s) reached in minimisation:"," "*13)+' '+' '.join(['{:>11.8f}'.format(i) for i in error.pars])+'\n{}\n'.format("Restarting minimisation with new x0"),
            log_file)
    
    # System exit if profile not converged
    if not result['success']:
        info_log(result['message'],log_file)
        sys.exit()
    
    # Converged parameters
    x_final = result['x']
    
    # Save end profile as "profile.in"
    x_final_full = write_profile(x_final,num_n,num_p,gamOne,box,"profile.in")
    
    # Print results
    info_log('{:>13}'.format("x_final_full:")+' '+' '.join(['{:>11.8f}'.format(i) for i in x_final_full]),log_file)
    info_log('{:>13} {:>5}'.format("etf calls:",result['nfev']),log_file)
    
    
    ## FINAL RUN ##
    # Change "input.in" to run with Strutinsky
    if strut:
        subprocess.run(["sed -i '/strutinsky_on = /s/.false. .false./.false. .true./g' input.in"], shell=True)
    
    # Change "input.in" to run with pairing perturbatively
    if pair_n_perturb:
        subprocess.run(["sed -i '/neutron_pairing = /s/[0-2]/2/g' input.in"], shell=True)
    
    if pair_p_perturb:
        subprocess.run(["sed -i '/proton_bcs = /s/false/true/g' input.in"], shell=True)
    
    if beta:
        # Change "input.in" to calculate chemical potentials
        subprocess.run(["sed -i '/calc_chem_pots = /s/false/true/g' input.in"], shell=True)
        
        # Run "etf" one final time, returning overall chemical potential
        delta_mu = run_etf(x_final,num_n,num_p,rho_b,gamOne,box,return_delta_mu=True)
        
        info_log('{:>13} {:>15.6f}'.format("E/A:",result['fun']),log_file)
        info_log('{:>13} {:>15.6f}\n'.format("delta_mu:",delta_mu),log_file)
        
    else:
        # Run "etf" one final time, returning total energy per particle
        e = run_etf(x_final,num_n,num_p,rho_b,gamOne,box,return_delta_mu=False)
        
        info_log('{:>13} {:>12.6f}\n'.format("e:",e),log_file)
        
    if beta:
        return delta_mu
    else:
        return e