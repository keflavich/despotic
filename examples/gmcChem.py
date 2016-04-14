"""
This is an example code for the DESPOTIC library. This code uses the
Nelson & Langer (1999) network to compute the chemical state of a
cloud as a function of its column density. The result is written out
to a file gmcChem.pdf.
"""

# Import the despotic library and the NL99 network; also import numpy
from despotic import cloud
from despotic.chemistry import NL99_GC
import matplotlib.pyplot as plt
import numpy as np

# Range in log column density to use
logNH = np.arange(21., 23.01, 0.025)

# Use the Milky Way GMC file as a base
gmc=cloud('cloudfiles/MilkyWayGMC.desp')

# Lower the CR ionization rate so that a fully CO composition becomes
# possible
gmc.rad.ionRate = 2e-17

# Raise the temperature to 20 K
gmc.Tg = 20.0

# Loop over column densities, getting chemical composition at each one
# and recording the value
abd=[]
for lNH in logNH:
    gmc.colDen = 10.**lNH
    gmc.setChemEq(network=NL99_GC, verbose=True)
    abd.append(gmc.chemnetwork.abundances)

# Make arrays for plotting
xCO = []
xC = []
xCp = []
for ab in abd:
    xCO.append(ab['CO'])
    xC.append(ab['C'])
    xCp.append(ab['C+'])
xCO=np.array(xCO)
xC=np.array(xC)
xCp=np.array(xCp)

# Plot
plt.plot(logNH, xCO+xC+xCp, 'k', linewidth=3, linestyle='dashed',
         label='Total C')
plt.plot(logNH, xCO, label=r'$x_{\rm CO}$', linewidth=3)
plt.plot(logNH, xC, label=r'$x_{\rm C}$', linewidth=3)
plt.plot(logNH, xCp, label=r'$x_{\rm C^+}$', linewidth=3)
plt.yscale('log')
plt.xlabel(r'$\log\,N_{\rm H}$')
plt.ylabel('Abundance')
plt.ylim([1e-7,1e-2])
plt.legend()
plt.savefig('gmcChem.pdf')
