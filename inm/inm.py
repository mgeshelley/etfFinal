import numpy as np
import pandas as pd
from scipy.optimize import minimize

import matplotlib.pyplot as plt



# Dutra et al. 2012
def H(n,y):
    return 2**(n-1.) * (y**n + (1.-y)**n)


# Dutra et al. 2012
def E_ANM(rho, y, parList, BSk_extra=False):
    # Unpack parameter list
    t0, t1, t2, t3, x0, x1, x2, x3, b4, b4p, alpha, amfc, t4, t5, x4, x5, delta, gamma = parList
    
    # Convenience parameters
    a = t1*(x1+2.) + t2*(x2+2.)
    b = 0.5 * (t2*(2.*x2 + 1.) - t1*(2.*x1 + 1.))
    
    result = (3./5 * hbar2_2m * f_pi * rho**fr23 * H(fr53,y)
              + t0 / 8. * rho * (2.*(x0+2.) - (2.*x0 + 1.)*H(2.,y))
              + t3 / 48. * rho**(alpha+1.) * (2.*(x3+2.) - (2.*x3 + 1.)*H(2.,y))
              + 3./40 * f_pi * rho**fr53 * (a*H(fr53,y) + b*H(fr83,y)))
    
    if BSk_extra:
        result = (result
                  + 3./40 * f_pi * rho**(fr53+delta) * (t4*(x4+2.)*H(fr53,y) - t4*(x4+0.5)*H(fr83,y))
                  + 3./40 * f_pi * rho**(fr53+gamma) * (t5*(x5+2.)*H(fr53,y) + t5*(x5+0.5)*H(fr83,y)))
    
    return result


# Dutra et al. 2012
def S_SNM(rho, parList, BSk_extra=False):
    # Unpack parameter list
    t0, t1, t2, t3, x0, x1, x2, x3, b4, b4p, alpha, amfc, t4, t5, x4, x5, delta, gamma = parList
    
    # Convenience parameters
    a = t1*(x1+2.) + t2*(x2+2.)
    b = 0.5 * (t2*(2.*x2 + 1.) - t1*(2.*x1 + 1.))
    
    result = (hbar2_2m/3. * f_pi * rho**fr23
              - t0/ 8. * (2.*x0 + 1.) * rho
              - t3 / 48. * (2.*x3 + 1.) * rho**(alpha+1.)
              + f_pi/24. * (a + 4.*b) * rho**fr53)
    
    if BSk_extra:
        result = (result
                  - f_pi/8. * t4*x4 * rho**(fr53+delta)
                  + f_pi/24. * t5*(5.*x5 + 4.) * rho**(fr53+gamma))
    
    return result


# Dutra et al. 2012
def L_SNM(rho0, parList, BSk_extra=False):
    # Unpack parameter list
    t0, t1, t2, t3, x0, x1, x2, x3, b4, b4p, alpha, amfc, t4, t5, x4, x5, delta, gamma = parList
    
    # Convenience parameters
    a = t1*(x1+2.) + t2*(x2+2.)
    b = 0.5 * (t2*(2.*x2 + 1.) - t1*(2.*x1 + 1.))
    
    result = (fr23*hbar2_2m * f_pi * rho0**fr23
              - 3./8 * t0 * (2.*x0 + 1.) * rho0
              - t3 / 16. * (2.*x3 + 1.) * (alpha+1.)*rho0**(alpha+1.)
              + 5./24 * f_pi * (a + 4.*b) * rho0**fr53)
    
    if BSk_extra:
        result = (result
                  - f_pi/8. * (5. + 3.*delta) * t4*x4 * rho0**(fr53+delta)
                  + f_pi/24. * (5. + 3.*gamma) * t5*(5.*x5 + 4.) * rho0**(fr53+gamma))
    
    return result


def INM_properties(parList):
    # Check if BSk_gamma is NaN or not
    if np.isnan(parList[-1]):
        BSk_extra = False
    else:
        BSk_extra = True
    
    # Minimise E/rho w.r.t rho, to get rho0 and E0
        # Start guess of 0.16 for rho0, with bounds of (0.13, 0.18)
    e_min_result = minimize(E_ANM, x0=0.16, bounds=[(0.13, 0.18)], args=(0.5, parList, BSk_extra))
    rho0 = e_min_result['x'][0]
    E0 = e_min_result['fun'][0]
    
    # Symmetry energy in SNM
    J = S_SNM(rho0, parList, BSk_extra=BSk_extra)
    
    # Slope of symmetry energy
    L = L_SNM(rho0, parList, BSk_extra=BSk_extra)
    
    return rho0, E0, J, L



## CONSTANTS ##
# Constants
fr23, fr53, fr83 = (2./3, 5./3, 8./3)
f_pi = (1.5*np.pi**2)**fr23

# Nucleon mass
hbar2_2m = 20.7355


## LOAD DATA ##
#file_name = '../../ale_utils/skyrme_forces_noReferences.csv'
#df = pd.read_csv(file_name, na_values='-', index_col='Skyrme')
#df.to_pickle('par_df.pkl')

# Read from Pandas data frame file
df = pd.read_pickle('par_df.pkl')


## CLEAN DATA ##
# Filter out forces with extra parameters
extra_pars = ['t32','t33','t4','x32','x33','x4','alpha2','alpha3','te','to']
df_stand = df[df[extra_pars].isnull().all(1)]

# Remove columns for extra parameters
df_stand = df_stand.drop(extra_pars, axis=1)

# Remove columns with extra information
extra_info_columns = ['exchange', 'ls', 'cm', 'CHG', 'ipair', 'pairfp', 'pairfn', 'pairfun']
df_stand = df_stand.drop(extra_info_columns, axis=1)

# Filter out 2-parameter spin-orbit
df_stand = df_stand[df_stand['b4'] == df_stand['b4p']]

# Ssk causes problems
df_stand = df_stand.drop(['Ssk'])


## BSk FORCES WITH EXTRA PARAMETERS ##
# Add extra parameters in BSk forces
extra_pars_BSk = ['BSk_t4', 'BSk_t5', 'BSk_x4', 'BSk_x5', 'BSk_delta', 'BSk_gamma']
df_stand = df_stand.reindex(columns=[*df_stand.columns.tolist(), *extra_pars_BSk], fill_value=np.nan)

df_stand.loc['BSk21'] = [-3961.39, 396.131, 1.e-3, 22588.2,
                         0.885231, 0.0648452, -1390380, 1.03928,
                         54.811, 54.811,
                         1./12,
                         np.nan,
                         -100.000, -150.000, 2.00000, -11.0000,
                         0.5, 1./12]

df_stand.loc['BSk24'] = [-3970.28718809175, 395.7662092819614, 1.e-5, 22648.57756432157,
                         0.8943714063887188, 0.05635350781437905, -138961080.022762, 1.051191623994705,
                         54.20264386, 54.20264386,
                         1./12,
                         np.nan,
                         -100.000, -150.000, 2.00000, -11.0000,
                         0.5, 1./12]


## INDIVIDUAL FORCE PROPERTIES ##
# List of forces for which to print and save properties
force_strings = ['SLy4', 'BSk21', 'BSk24', 'LNS']

for force in force_strings:
    # Extract list of parameters for 'force'
    parList = df_stand.loc[force].values.tolist()

    # Print INM properties
    print(force, INM_properties(parList))

    # Grid of densities
    rho_min = 0.
    rho_max = 1.
    rho_step = 0.001
    rho_grid = np.linspace(rho_min, rho_max, int((rho_max-rho_min)/rho_step) + 1)

    # Check if BSk_gamma is NaN or not
    if np.isnan(parList[-1]):
        BSk_extra = False
    else:
        BSk_extra = True

    # Calculate E/A for PNM
    e_grid = E_ANM(rho_grid, 0., parList, BSk_extra=BSk_extra)

    # Save to file
    np.savetxt('PNM_'+force+'.dat', np.column_stack((rho_grid, e_grid)))


## INM FOR ALL FORCES ##
# New column with all INM properties in one tuple on each row
df_stand['inm_tuple'] = df_stand.iloc[:, :].apply(INM_properties, axis=1)

# Unpack tuple to four separate columns
INM_strings = ['rho0', 'E0', 'J', 'L']
df_stand[INM_strings] = pd.DataFrame(df_stand['inm_tuple'].tolist(), index=df_stand.index)

# Delete tuple column
df_stand = df_stand.drop(['inm_tuple'], axis=1)

# Print selection
print(df_stand[df_stand['rho0'] < 0.15][INM_strings])
print(df_stand[df_stand['E0'] < -16.][INM_strings])
print(df_stand[df_stand['J'] > 35.][INM_strings])
print(df_stand[df_stand['L'] > 100.][INM_strings])

# Histogram of INM properties
fig, axs = plt.subplots(2, 2)
df_stand.hist(column=INM_strings, ax=axs, bins=20)
fig.savefig('histogram.pdf')