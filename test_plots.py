#!/usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import plotting_functions as matt_plot




## OPTIONS ##

# Plotting parameters
axisLabelSize = 25
tickLabelSize = 25
minorTickLength = 10
majorTickLength = 25
lineWidth = 2
tickWidth = 2
plotParList = [axisLabelSize,tickLabelSize,minorTickLength,majorTickLength,lineWidth,tickWidth]

# Line style list
lineList11 = ['-k']
lineList21 = ['-k','--b']
lineList31 = ['-k','--b','-.g']
lineList41 = ['-k','--b',':r','-.g']

lineList22 = ['-k','--k','-b','--b']
lineList42 = ['-k','--k','-b','--b','-r','--r','-g','--g']

lineList2n2p = ['-k','--k','-r','--r']
lineList2n = ['-k','--k']

# HF data to load for tests
    # 1 --> Z=40, N=40000, SLy4
    # 2 --> Pb208, SLy4
sample_profile = 2



## TESTS ##

# Load test data
etf_rho = np.loadtxt('test_densities.dat')
etf_tau = np.loadtxt('test_kinetic_densities.dat')
etf_spin = np.loadtxt('test_spin_current_densities.dat')
etf_tau_order = np.loadtxt('order_by_order_kinetic_densities.dat')

if sample_profile == 1:
    hfb_rho_n = np.loadtxt('ale_profiles/Z40_N40000_SLy4_Rbox_80/DensitiesN.dat')
    hfb_rho_p = np.loadtxt('ale_profiles/Z40_N40000_SLy4_Rbox_80/DensitiesP.dat')
    
elif sample_profile == 2:
    hfb_rho_n = np.loadtxt('ale_profiles/Pb208_SLy4/neutron_densities.dat')
    hfb_rho_p = np.loadtxt('ale_profiles/Pb208_SLy4/proton_densities.dat')
    hfb_n_fields = np.loadtxt('ale_profiles/Pb208_SLy4/neutron_fields.dat')
    hfb_p_fields = np.loadtxt('ale_profiles/Pb208_SLy4/proton_fields.dat')
    
    etf_u = np.loadtxt('test_central_potentials.dat')
    
    # Use correct form factor for HFB spin-orbit fields
    HFB_spin_n = hfb_n_fields[:,2]*hfb_n_fields[:,0]*0.5
    HFB_spin_p = hfb_p_fields[:,2]*hfb_p_fields[:,0]*0.5
    
    # ETF spin-orbit fields
    ETF_spin_n = etf_spin[:,3]
    ETF_spin_p = etf_spin[:,4]
    
    # HFB central potentials
    HFB_U_n = hfb_n_fields[:,1]
    HFB_U_p = hfb_p_fields[:,1]
    
    # ETF central potentials
    ETF_U_n = etf_u[:,1]
    ETF_U_p = etf_u[:,2]
    
else:
    print('No sample profile for comparisons')

HFB_J_n = hfb_rho_n[:,3]*-1
HFB_J_p = hfb_rho_p[:,3]*-1
ETF_J_n = etf_spin[:,1]
ETF_J_p = etf_spin[:,2]



xlimit = 12.0

## MATTER DENSITIES COMPARISON WITH HF ##

labelList = ['neutrons (HF)','neutrons (fitted)','protons (HF)','protons (fitted)']
xdataList = [hfb_rho_n[:,0],etf_rho[:,0],hfb_rho_p[:,0],etf_rho[:,0]]
ydataList = [hfb_rho_n[:,1],etf_rho[:,2],hfb_rho_p[:,1],etf_rho[:,3]]

matt_plot.cell_radius_plot_old(xdataList,ydataList,xlimit,r"$\rho_{q}$ [fm$^{-3}$]",plotParList,labelList,lineList2n2p,'test_densities')


## KINETIC DENSITIES COMPARISON WITH HF ##

labelList = ['neutrons (HF)','neutrons (ETF)','protons (HF)','protons (ETF)']
xdataList = [hfb_rho_n[:,0],etf_tau[:,0],hfb_rho_p[:,0],etf_tau[:,0]]
ydataList = [hfb_rho_n[:,2],etf_tau[:,2],hfb_rho_p[:,2],etf_tau[:,3]]

matt_plot.cell_radius_plot_old(xdataList,ydataList,xlimit,r"$\tau_{q}$ [fm$^{-5}$]",plotParList,labelList,lineList2n2p,'test_kinetic_densities')


## SPIN CURRENT DENSITIES COMPARISON WITH HF ##

labelList = ['neutrons (HF)','neutrons (ETF)','protons (HF)','protons (ETF)']
xdataList = [hfb_rho_n[:,0],etf_spin[:,0],hfb_rho_p[:,0],etf_spin[:,0]]
ydataList = [HFB_J_n,ETF_J_n,HFB_J_p,ETF_J_p]

matt_plot.cell_radius_plot_old(xdataList,ydataList,xlimit,r"$J_{q}$ [fm$^{-4}$]",plotParList,labelList,lineList2n2p,'test_spin_current_densities')

# Spin-orbit fields not available for WS cell test
if sample_profile == 2:
    xdataList = [hfb_n_fields[:,0],etf_spin[:,0],hfb_p_fields[:,0],etf_spin[:,0]]
    ydataList = [HFB_spin_n,ETF_spin_n,HFB_spin_p,ETF_spin_p]
    
    matt_plot.cell_radius_plot_old(xdataList,ydataList,xlimit,r"$W_{q}$ [fm$^{-1}$]",plotParList,labelList,lineList2n2p,'test_spin_orbit_fields')


## ORDER BY ORDER KINETIC DENSITIES COMPARISON WITH HF - NEUTRONS ##

labelList = ['HF','TF','ETF(2)','ETF(4)']
xdataList = [hfb_rho_n[:,0],etf_tau_order[:,0],etf_tau_order[:,0],etf_tau_order[:,0]]
ydataList = [hfb_rho_n[:,2],etf_tau_order[:,1],etf_tau_order[:,2],etf_tau_order[:,3]]

matt_plot.cell_radius_plot_old(xdataList,ydataList,xlimit,r"$\tau_{n}$[fm$^{-5}$]",
                 plotParList,labelList,lineList41,'test_order_by_order_kinetic_neutrons')


## ORDER BY ORDER KINETIC DENSITIES COMPARISON WITH HF - PROTONS ##

labelList = ['HF','TF','ETF(2)','ETF(4)']
xdataList = [hfb_rho_p[:,0],etf_tau_order[:,0],etf_tau_order[:,0],etf_tau_order[:,0]]
ydataList = [hfb_rho_p[:,2],etf_tau_order[:,4],etf_tau_order[:,5],etf_tau_order[:,6]]

matt_plot.cell_radius_plot_old(xdataList,ydataList,xlimit,r"$\tau_{p}$[fm$^{-5}$]",
                 plotParList,labelList,lineList41,'test_order_by_order_kinetic_protons')


## CENTRAL POTENTIALS COMPARISON WITH HF ##
# Central potentials not available for WS cell test
if sample_profile == 2:
    labelList = ['neutrons (HF)','neutrons (ETF)','protons (HF)','protons (ETF)']
    xdataList = [hfb_n_fields[:,0],etf_u[:,0],hfb_p_fields[:,0],etf_u[:,0]]
    ydataList = [HFB_U_n,ETF_U_n,HFB_U_p,ETF_U_p]
    
    matt_plot.cell_radius_plot_old(xdataList,ydataList,xlimit,r"$U_{q}$[MeV]",
                    plotParList,labelList,lineList2n2p,'test_central_potentials',legPos='lower right')



##Plots of the eigenvalues
#EigN = np.loadtxt('neutron_singleparticles.dat')
#EigP = np.loadtxt('proton_singleparticles.dat')
#fig,ax = plt.subplots(1,2)
#ax[0].hlines(EigN[:,0],0,1)
#ax[0].set_title('Neutron states')
#plt.setp(ax[0].get_xticklabels(), visible=False)
#ax[1].hlines(EigP[:,0],0,1)
#ax[1].set_title('Proton states')
#plt.setp(ax[1].get_xticklabels(), visible=False)
#plt.savefig('eigen.pdf')



## BENCHMARKS ##

    ## BARTEL AND BENCHEIKH (2002) ##
bart_ben = np.loadtxt('bartel_bencheikh_benchmark.dat')
bart_ben_extra = np.loadtxt('bartel_bencheikh_extra.dat')

hfbrad = np.loadtxt('benchmark_data/neutron_126_82.dens')

bb_paper_tau_n_2_rho = np.loadtxt('benchmark_data/bartel_bencheikh_tau_n(2)_rho.dat')
bb_paper_tau_n_2_f = np.loadtxt('benchmark_data/bartel_bencheikh_tau_n(2)_f.dat')
bb_paper_tau_n_2_rho_f = np.loadtxt('benchmark_data/bartel_bencheikh_tau_n(2)_rho_f.dat')
bb_paper_tau_n_2_so = np.loadtxt('benchmark_data/bartel_bencheikh_tau_n(2)_s-o.dat')

bb_paper_tau_n_4_rho = np.loadtxt('benchmark_data/bartel_bencheikh_tau_n(4)_rho.dat')
bb_paper_tau_n_4_f = np.loadtxt('benchmark_data/bartel_bencheikh_tau_n(4)_f.dat')
bb_paper_tau_n_4_rho_f = np.loadtxt('benchmark_data/bartel_bencheikh_tau_n(4)_rho_f.dat')
bb_paper_tau_n_4_so = np.loadtxt('benchmark_data/bartel_bencheikh_tau_n(4)_s-o.dat')

bb_paper_J_n_2 = np.loadtxt('benchmark_data/bartel_bencheikh_J_n(2).dat')
bb_paper_J_n_4 = np.loadtxt('benchmark_data/bartel_bencheikh_J_n(4).dat')

bb_paper_rho_n = np.loadtxt('benchmark_data/bartel_bencheikh_rho_n.dat')
bb_paper_rho_p = np.loadtxt('benchmark_data/bartel_bencheikh_rho_p.dat')
bb_paper_tau_n_TF = np.loadtxt('benchmark_data/bartel_bencheikh_tau_n_TF.dat')
bb_paper_tau_n_2 = np.loadtxt('benchmark_data/bartel_bencheikh_tau_n_2.dat')
bb_paper_tau_n_4 = np.loadtxt('benchmark_data/bartel_bencheikh_tau_n_4.dat')
bb_paper_U_n = np.loadtxt('benchmark_data/bartel_bencheikh_U_n.dat')
bb_paper_U_p = np.loadtxt('benchmark_data/bartel_bencheikh_U_p.dat')

# Maximum radius to plot
xlimit = 12.0

# Fig. 2. (a)
labelList = [r"$\tau^{(0)}_{TF}$ (MS)",r"$\tau^{(0)}_{TF}$ (Bartel)",r"$\tau^{(2)}_{ETF}$ (MS)",r"$\tau^{(2)}_{ETF}$ (Bartel)",r"$\tau^{(4)}_{ETF}$ (MS)",r"$\tau^{(4)}_{ETF}$ (Bartel)"]
xdataList = [bart_ben[:,0],bb_paper_tau_n_TF[:,0],bart_ben[:,0],bb_paper_tau_n_2[:,0],bart_ben_extra[:,0],bb_paper_tau_n_4[:,0]]
ydataList = [bart_ben[:,1],bb_paper_tau_n_TF[:,1],bart_ben[:,2],bb_paper_tau_n_2[:,1],bart_ben_extra[:,5]*10,bb_paper_tau_n_4[:,1]]

matt_plot.cell_radius_plot_old(xdataList,ydataList,xlimit,r"$\tau_{n}$ [fm$^{-5}$]",
                plotParList,labelList,lineList42,'bart_ben_tau_n_order_contributions',ylims=[-0.025,0.15],y0line=True)

# Fig. 2. (b)
labelList = [r"$\nabla\rho$ terms (MS)",r"$\nabla\rho$ terms (Bartel)",
            r"$\nabla f$ terms (MS)",r"$\nabla f$ terms (Bartel)",
            r"$\nabla\rho\cdot\nabla f$ terms (MS)",r"$\nabla\rho\cdot\nabla f$ terms (Bartel)",
            'spin-orbit terms (MS)','spin-orbit terms (Bartel)']
xdataList = [bart_ben[:,0],bb_paper_tau_n_2_rho[:,0],bart_ben[:,0],bb_paper_tau_n_2_f[:,0],bart_ben[:,0],bb_paper_tau_n_2_rho_f[:,0],bart_ben[:,0],bb_paper_tau_n_2_so[:,0]]
ydataList = [bart_ben[:,3],bb_paper_tau_n_2_rho[:,1],bart_ben[:,4],bb_paper_tau_n_2_f[:,1],bart_ben[:,5],bb_paper_tau_n_2_rho_f[:,1],bart_ben[:,6],bb_paper_tau_n_2_so[:,1]]

matt_plot.cell_radius_plot_old(xdataList,ydataList,xlimit,r"$\tau^{(2)}_{n}$ [fm$^{-5}$]",plotParList,labelList,lineList42,
                'bart_ben_tau_n_2nd_order_contributions',ylims=[-0.015,0.015],y0line=True,legPos='upper left')

# Fig. 2. (c)
labelList = [r"$\nabla\rho$ terms (MS)",r"$\nabla\rho$ terms (Bartel)",
            r"$\nabla f$ terms (MS)",r"$\nabla f$ terms (Bartel)",
            r"$\nabla\rho\cdot\nabla f$ terms (MS)",r"$\nabla\rho\cdot\nabla f$ terms (Bartel)",
            'spin-orbit terms (MS)','spin-orbit terms (Bartel)']
xdataList = [bart_ben_extra[:,0],bb_paper_tau_n_4_rho[:,0],bart_ben_extra[:,0],bb_paper_tau_n_4_f[:,0],bart_ben_extra[:,0],bb_paper_tau_n_4_rho_f[:,0],bart_ben_extra[:,0],bb_paper_tau_n_4_so[:,0]]
ydataList = [bart_ben_extra[:,6],bb_paper_tau_n_4_rho[:,1],bart_ben_extra[:,7],bb_paper_tau_n_4_f[:,1],bart_ben_extra[:,8],bb_paper_tau_n_4_rho_f[:,1],bart_ben_extra[:,9],bb_paper_tau_n_4_so[:,1]]

matt_plot.cell_radius_plot_old(xdataList,ydataList,xlimit,r"$\tau^{(4)}_{n}$ [fm$^{-5}$]",plotParList,labelList,lineList42,
                'bart_ben_tau_n_4th_order_contributions',ylims=[-0.002,0.001],y0line=True,legPos='upper left')

# Fig. 3.
labelList = [r"$J^{(2)}_{ETF}$ (MS)",r"$J^{(2)}_{ETF}$ (Bartel)",r"$J^{(4)}_{ETF}$ (MS)",r"$J^{(4)}_{ETF}$ (Bartel)"]
xdataList = [bart_ben[:,0],bb_paper_J_n_2[:,0],bart_ben[:,0],bb_paper_J_n_4[:,0]]
ydataList = [bart_ben[:,7],bb_paper_J_n_2[:,1],bart_ben[:,8],bb_paper_J_n_4[:,1]]

matt_plot.cell_radius_plot_old(xdataList,ydataList,xlimit,r"$J_{n}$ [fm$^{-4}$]",
                plotParList,labelList,lineList22,'bart_ben_j_n_order_contributions',ylims=[-0.01,0.02],y0line=True)

# Fig. 1.
labelList = ['neutrons (MS)','neutrons (Bartel)','protons (MS)','protons (Bartel)']
xdataList = [bart_ben_extra[:,0],bb_paper_rho_n[:,0],bart_ben_extra[:,0],bb_paper_rho_p[:,0]]
ydataList = [bart_ben_extra[:,1],bb_paper_rho_n[:,1],bart_ben_extra[:,2],bb_paper_rho_p[:,1]]

matt_plot.cell_radius_plot_old(xdataList,ydataList,xlimit,r"$\rho_{q}$ [fm$^{-3}$]",plotParList,labelList,lineList2n2p,
                'bart_ben_rho_comparison',ylims=[0.0,0.1],y0line=True,legPos='upper right')

# Fig. 6
labelList = ['neutrons (MS)','neutrons (Bartel)','protons (MS)','protons (Bartel)']
xdataList = [bart_ben_extra[:,0],bb_paper_U_n[:,0],bart_ben_extra[:,0],bb_paper_U_p[:,0]]
ydataList = [bart_ben_extra[:,3],bb_paper_U_n[:,1],bart_ben_extra[:,4],bb_paper_U_p[:,1]]

matt_plot.cell_radius_plot_old(xdataList,ydataList,xlimit,r"$U_{q}$ [MeV]",plotParList,labelList,lineList2n2p,
                'bart_ben_U_comparison',ylims=[-80.0,0.0],y0line=True,legPos='upper left')