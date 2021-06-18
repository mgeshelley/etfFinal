#!/usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter



# Plotting function for quantity across Wigner-Seitz cell
def cell_radius_plot(xdataList,ydataList,xlim,ylab,labelList,lineStyle,colourList,output,ylimbottom0=False,ylims=None,y0line=False,
                     legPos='upper right'):
    # Set up plot
        # Dimensions for plotting canvas
    xinch = 3.4
    yinch = 3.4
        # Margins: left, right, top, bottom
    margins = [0.11, 0.04, 0.11, 0.04]
        # Fraction of image's x- and y-dimension that canvas covers
    xfrac = 1. - margins[0] - margins [1]
    yfrac = 1. - margins[2] - margins [3]
    
    # Create figure
    #fig = plt.figure()
    fig = plt.figure(figsize=(xinch/xfrac,yinch/yfrac))
    
    # Add axes (plotting canvas)
    #ax = ax.add_subplot(1, 1, 1)
    ax = plt.axes([margins[0], margins[2], xfrac, yfrac])
    
    # Labels
    #rcParams['mathtext.default'] = 'regular' # Non-italic font for plots
    ax.set_xlabel(r"$r$ [fm]")
    ax.set_ylabel(ylab)
    
    # Line at y=0
    if y0line:
        ax.axhline(color='k')
    
    # Axes
    ax.minorticks_on()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g')) # No unnecessary trailing decimals
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    
    # Plot
    for i in range(len(xdataList)):
        ax.plot(xdataList[i],ydataList[i],linestyle=lineStyle[i],color=colourList[i],clip_on=True,label=labelList[i])
    
    ax.set_xlim(0.0,xlim)
    
    # y-limit
    if ylimbottom0:
        ax.set_ylim(bottom=0)
    
    if ylims:
        ax.set_ylim(ylims)
    
    legend = ax.legend(loc=legPos)
    
    ax.grid()
    
    #ax.yaxis.set_major_locator(MultipleLocator(0.004))
    #ax.yaxis.set_major_locator(plt.MaxNLocator(6))
    
    #fig.tight_layout()
    fig.savefig(output+'.pdf', format='pdf')



# Old plotting function for quantity across Wigner-Seitz cell, that works with "testplots.py"
def cell_radius_plot_old(xdataList,ydataList,xlim,ylab,plotParList,labelList,lineList,output,ylimbottom0=False,ylims=None,y0line=False,
                     legPos='upper right'):
    
    rcParams['axes.linewidth'] = plotParList[5]
    
    # Set up plot
    xinch=14
    yinch=14
    #fig = plt.figure()
    fig = plt.figure(figsize=(xinch/0.8,yinch/0.9))
    #ax = fig.add_subplot(1, 1, 1)
    ax = plt.axes([0.11, 0.1, 0.85, 0.85])
    
    # Labels
    #rcParams['mathtext.default'] = 'regular' # Non-italic font for plots
    plt.xlabel(r"$r$ [fm]",fontsize=plotParList[0])
    plt.ylabel(ylab,fontsize=plotParList[0])

    # Line at y=0
    if y0line:
        plt.axhline(color='k')
    
    # Axes
    plt.tick_params(which='both', labelsize=plotParList[1], width=plotParList[5])
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g')) # No unnecessary trailing decimals
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.minorticks_on()
    ax.tick_params(which='minor',length=plotParList[2])
    ax.tick_params(which='major',length=plotParList[3])
    
    # Plot
    for i in range(len(xdataList)):
        plt.plot(xdataList[i],ydataList[i],lineList[i],lw=plotParList[4],clip_on=True,label=labelList[i])
    
    plt.xlim(0.0,xlim)
    
    # y-limit
    if ylimbottom0:
        plt.ylim(bottom=0)
    
    if ylims:
        plt.ylim(ylims)
    
    plt.legend(loc=legPos, prop={'size': plotParList[0]}, frameon=False)
    
    #ax.yaxis.set_major_locator(MultipleLocator(0.004))
    #ax.yaxis.set_major_locator(plt.MaxNLocator(6))

    #plt.tight_layout()
    plt.savefig(output+'.eps', format='eps', dpi=1000)