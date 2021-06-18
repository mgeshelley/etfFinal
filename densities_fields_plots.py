#!/usr/bin/python3

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import plotting_functions as matt_plot

plt.style.use("style_1.mplstyle")


## OPTIONS ##
# Line style list
lineStyle41 = ['-',':','--','-.']

# Colours from colour map "tab10"; blue, green, red
colour_map_colours = [plt.cm.tab10(i) for i in [0,3,2]]

# Add black as first colour
colour_list = ['k'] + colour_map_colours


## LOAD DATA ##
densities_n = np.loadtxt('densities_n.dat')
densities_p = np.loadtxt('densities_p.dat')
fields_n = np.loadtxt('fields_n.dat')
fields_p = np.loadtxt('fields_p.dat')


## DENSITIES ##
xlimit = 15.0

# Neutrons
labelList = [r"$\rho_{n} \left[\mathrm{fm}^{-3}\right]$",r"$\tau_{n} \left[\mathrm{fm}^{-5}\right]$",r"$J_{n}(\times10) \left[\mathrm{fm}^{-4}\right]$"]
xdataList = [densities_n[:,0],densities_n[:,0],densities_n[:,0]]
ydataList = [densities_n[:,1],densities_n[:,2],densities_n[:,3]*10.]
matt_plot.cell_radius_plot(xdataList,ydataList,xlimit,"",labelList,lineStyle41,colour_list,'densities_n',ylims=(-0.005,0.13))

xlimit = 15.0

# Protons
labelList = [r"$\rho_{p} \left[\mathrm{fm}^{-3}\right]$",r"$\tau_{p} \left[\mathrm{fm}^{-5}\right]$",r"$J_{p}(\times10) \left[\mathrm{fm}^{-4}\right]$"]
xdataList = [densities_p[:,0],densities_p[:,0],densities_p[:,0]]
ydataList = [densities_p[:,1],densities_p[:,2],densities_p[:,3]*10.]
matt_plot.cell_radius_plot(xdataList,ydataList,xlimit,"",labelList,lineStyle41,colour_list,'densities_p',ylims=(-0.005,0.13))


## FIELDS ##
xlimit = 15.0

# Neutrons
labelList = [r"$U_{n} \left[\mathrm{MeV}\right]$",r"$\frac{\hbar^2}{2m^*_n} \left[\mathrm{MeV}\cdot \mathrm{fm}^2\right]$",r"$W_{n}(\times10) \left[\mathrm{MeV}\cdot \mathrm{fm}\right]$"]
xdataList = [fields_n[:,0],fields_n[:,0],fields_n[:,0]]
ydataList = [fields_n[:,1],fields_n[:,2],fields_n[:,3]*10.]
matt_plot.cell_radius_plot(xdataList,ydataList,xlimit,"",labelList,lineStyle41,colour_list,'fields_n',ylims=(-80.,35.),legPos='lower right')

xlimit = 15.0

# Protons
labelList = [r"$U_{p} \left[\mathrm{MeV}\right]$",r"$\frac{\hbar^2}{2m^*_p} \left[\mathrm{MeV}\cdot \mathrm{fm}^2\right]$",r"$W_{p}(\times10) \left[\mathrm{MeV}\cdot \mathrm{fm}\right]$"]
xdataList = [fields_p[:,0],fields_p[:,0],fields_p[:,0]]
ydataList = [fields_p[:,1],fields_p[:,2],fields_p[:,3]*10.]
matt_plot.cell_radius_plot(xdataList,ydataList,xlimit,"",labelList,lineStyle41,colour_list,'fields_p',ylims=(-80.,35.),legPos='lower right')