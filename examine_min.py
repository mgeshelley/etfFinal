import numpy as np
from matplotlib import pyplot as plt

import subprocess
from multiprocessing import Pool

import emcee
import corner

from minimisation_functions import f_profile, write_profile, penalty_function



# Definition of log-probability function for MCMC sampler
def log_prob(p,num_n,num_p,gamOne,box):
    energy = penalty_function(p,num_n,num_p,gamOne,box)
    
    # If "etf" returns NaN, log-probability should be -infinity
    if np.isnan(energy):
        print("With parameters:", p, "'etf' returned 'NaN'")
        return -np.inf
    
    # Otherwise return -1/2 * energy**2 (Gaussian likelihood)
    else:
        return -0.5 * energy**2



## OPTIONS ##
# Do sampling or plotting?
sample = True

# N and Z
num_n = 21850
num_p = 50

# WS cell radius
box = 60.

# Gamma = 1?
gamOne = True

# More walkers than parameters
nwalkers = 500

# Chain length
nchain = 20000

# Show progress bar
progress_bar = False


## WALKER START POSITIONS ##
# Minimum location (no Strutinsky, pairing, or electrons)
x0 = np.array([2.3946695867757706e-02, 7.9989076561187709e+00, 8.8631791966191342e-01,   7.5542465332541404e+00, 5.6803316934008130e-01])
#x0 = np.array([-7.5312228755220794, 0.55050433336346438])

# Number of parameters
ndim = x0.shape[0]

print('Dimensions = {:d}'.format(ndim))
print('Walkers = {:d}'.format(nwalkers))
print('Chain length = {:d}'.format(nchain))
print('Function calls = {:d}'.format(nwalkers*nchain))

# Initialise walkers in "small Gaussian ball" around minimum
p0 = x0 + 1e-5 * np.random.randn(nwalkers,ndim)


if sample:
    ## SAMPLE ##
    # Make sure "etf" runs in run_mode 1
    subprocess.run(["sed -i '/run_mode = /s/[1-4]/1/g' input.in"], shell=True)

    # Make sure "etf" runs without Strutinsky
    subprocess.run(["sed -i '/strutinsky_on = /s/true/false/g' input.in"], shell=True)

    # Recompile with "make fast" to ensure "etf" not compiled with -fopenmp
    subprocess.run(["make fast"], shell=True)

    # Enable (shared-memory) parallelisation
    with Pool() as pool:
        # Create sampler
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args=[num_n,num_p,gamOne,box], pool=pool)
        
        # Run sampler with progress bar
        sampler.run_mcmc(p0, nchain, progress=progress_bar)

    # Autocorrelation time
    tau = sampler.get_autocorr_time()

    # Number of samples to discard from beginning of chain
    burnin = int(2 * np.max(tau))

    # How much to thin chain
    thin = int(0.5 * np.min(tau))

    print('tau = ', tau)
    print('Burn-in = {:d}'.format(burnin))
    print('Thin = {:d}'.format(thin))

    # 3D array of samples
    samples = sampler.get_chain()

    # Save all samples for reuse
    np.save("samples",samples)

    # Remove "burn-in" samples, and thin samples so not too many correlated samples are plotted
    flat_samples = sampler.get_chain(discard=burnin, thin=thin, flat=True)

    # Save modified samples for plotting again
    np.savetxt("flat_samples.dat",flat_samples)

    # Extract log probabilites of chain
    log_prob_samples = sampler.get_log_prob(discard=burnin, thin=thin, flat=True)

    # Save log-probs for plotting
    np.savetxt("log_probs.dat",log_prob_samples)

else:
    ## PLOTS ##
    # Create figure of chain(s)
    fig, axes = plt.subplots(ndim, figsize=(10,ndim*3), sharex=True)

    # Load samples for plotting
    samples = np.load("samples.npy")

    # Labels
    if ndim == 5:
        # 5-par minimisation
        labels = ["rho_n_gas","r_n","a_n","r_p","a_p"]
    elif ndim == 2:
        # 2-par minimisation
        labels = ["r_p","a_p"]

    # Plot value(s) of chain(s) for each parameter
    for i in range(ndim):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        #ax.yaxis.set_label_coords(-0.1, 0.5)
        ax.set_rasterized(True)

    axes[-1].set_xlabel("step number")

    fig.savefig("chains.eps", format='eps')


    # Load samples for plotting
    flat_samples = np.loadtxt("flat_samples.dat")

    # Corner plot
    fig = corner.corner(flat_samples,
                        bins=20,
                        labels=labels,
                        label_kwargs={"fontsize": 22},
                        show_titles=True,
                        title_fmt='.4f',
                        title_kwargs={"fontsize": 16},
                        truths=x0,
                        quantiles=[0.16, 0.5, 0.84],
                        max_n_ticks=5,
                        levels=(1-np.exp(-0.5),)
                        )

    fig.savefig("corner.eps", format='eps')


    # Begin new figure
    fig, ax = plt.subplots(1)

    # Load energies
    log_prob_samples = np.loadtxt("log_probs.dat")

    # Create histogram of energies
    ax.hist(np.sqrt(log_prob_samples*-2.))

    fig.savefig("energy.eps", format='eps')


### 2D CONTOUR PLOT OF SURFACE ##
## Grid
#mesh_size = 20

## Determine how far from x0 bounds are
#sigma = [0.0003, 0.0001]

## Bounds for plot x0 +/- (sigma*n_sig)
##n_sig = 3000
##bounds = np.array([[x0[0]-sigma[0]*n_sig, x0[0]+sigma[0]*n_sig], [x0[1]-sigma[1]*n_sig, x0[1]+sigma[1]*n_sig]])
#bounds = np.array([[7.2, 7.8], [0.4, 0.8]])

## Grids for contour plot
#X1 = np.linspace(bounds[0][0], bounds[0][1], mesh_size)
#X2 = np.linspace(bounds[1][0], bounds[1][1], mesh_size)
#x1, x2 = np.meshgrid(X1, X2)

## 2D surface
#surf_2D = [[penalty_function(np.array([x1[i,j],x2[i,j]]),num_n,num_p,gamOne,box) for i in range(x1.shape[0])] for j in range(x1.shape[1])]
#surf_2D = np.array(surf_2D).T
#np.savetxt('surface_2D.dat',surf_2D)
#surf_2D = np.loadtxt('surface_2D.dat')
    
## Plot
#plt.xlabel(r'$r_p$')
#plt.ylabel(r'$a_p$')
#plt.contourf(x1, x2, surf_2D)
#plt.colorbar()
#plt.plot(x0[0], x0[1], 'ro')
#plt.savefig('surface_2D.eps', format='eps')