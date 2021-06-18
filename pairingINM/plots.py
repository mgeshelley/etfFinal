#!/usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter



# Load data
data1 = np.loadtxt('gap_sharp.dat')
data2 = np.loadtxt('gap.dat')
#data3 = np.loadtxt('Campana_PNM.dat')

max_pair = np.argmax(data2[:,2])
print(data2[max_pair,1])
print(data2[max_pair,2])


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

# Labels
rcParams['mathtext.default'] = 'regular' # Non-italic font for plots
plt.xlabel(r"k$_{Fn}$ [fm$^{-1}$]",fontsize=25)
plt.ylabel(r"Pairing gap $\Delta_n$ [MeV]",fontsize=25)

# Axes
plt.xlim(0.0,np.max(data1[:,1]))
#plt.ylim(0.0,np.max(data1[:,2]*1.05))
plt.tick_params(labelsize=25)
ax.xaxis.set_major_formatter(FormatStrFormatter('%g')) # No unnecessary trailing decimals
ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
ax.minorticks_on()
ax.tick_params(which='minor',length=6)
ax.tick_params(which='major',length=15)

# Pairing gap
#plt.plot(data3[:,1],data3[:,2],'-b',lw=4,clip_on=False)
plt.plot(data2[:,1],data2[:,2],'--r',lw=4,clip_on=False)
#plt.plot(data3[:,1],data3[:,2],'-.g',lw=4,clip_on=False)
# Iterations
#plt.plot(data1[:,1],data1[:,3]/100.0,'-rx')

plt.tight_layout()
plt.savefig('gap.eps', format='eps', dpi=1000)
