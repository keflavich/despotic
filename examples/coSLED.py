########################################################################
# 
# coSLED.py
#
# This is an example code for the DESPOTIC library. This code computes
# the spectral line energy distributions for ^12CO and ^13CO in a Milky
# Way molecular cloud and for a ULIRG-like molecular cloud, then makes
# a graph of the result and saves it in the file coSLED.eps.
#
########################################################################

########################################################################
# Program code
########################################################################

# Import the despotic library
from despotic import cloud

# Import standard python libraries
from numpy import *
from matplotlib.pyplot import *
from datetime import datetime
from datetime import timedelta

# Read the Milky Way GMC cloud file
gmc = cloud(fileName='cloudfiles/MilkyWayGMC.desp')

# Read the ULIRG cloud file
ulirg = cloud(fileName='cloudfiles/ULIRG.desp')

# Compute the luminosity of the CO lines in both clouds
t1=datetime.now()
gmclines = gmc.lineLum('co')
ulirglines = ulirg.lineLum('co')
gmclines13 = gmc.lineLum('13co')
ulirglines13 = ulirg.lineLum('13co')
t2=datetime.now()
print('Execution time = '+str(t2-t1))

# Print out the CO X factor for both clouds. This is column density
# divided by velocity-integrated brightness temperature.
print("GMC X_CO = "+str(gmc.colDen/gmclines[0]['intTB']) + 
      " cm^-2 / (K km s^-1)")
print("ULIRG X_CO = "+str(ulirg.colDen/ulirglines[0]['intTB']) + 
      " cm^-2 / (K km s^-1)")

# Extract the integrated brightness temperatures
gmcTBint=array([line['intTB'] for line in gmclines])
ulirgTBint=array([line['intTB'] for line in ulirglines])
gmcTBint13=array([line['intTB'] for line in gmclines13])
ulirgTBint13=array([line['intTB'] for line in ulirglines13])

# Mask negative values to avoid warnings. Negative values are obtained
# for some of the very high J states due simply to roundoff error.
gmcTBint[gmcTBint <= 0] = 1e-50
ulirgTBint[ulirgTBint <= 0] = 1e-50
gmcTBint13[gmcTBint13 <= 0] = 1e-50
ulirgTBint13[ulirgTBint13 <= 0] = 1e-50

# Make bar plots
fig=figure(1)
ax=fig.add_subplot(111)
ax.bar(arange(len(ulirgTBint))+0.1, 1e20*ulirgTBint/ulirg.colDen, \
           width=0.2, \
           color='g', bottom=1e-6, label='ULIRG $^{12}$CO')
ax.bar(arange(len(ulirgTBint13))+0.1+0.2, 1e20*ulirgTBint13/ulirg.colDen, \
           width=0.2, \
           color='lightgreen', bottom=1e-6, label=r'ULIRG $^{13}$CO')
ax.bar(arange(len(gmcTBint))+0.1+0.4, 1e20*gmcTBint/gmc.colDen, \
           width=0.2, bottom=1e-6, \
           label='Milky Way $^{12}$CO')
ax.bar(arange(len(gmcTBint13))+0.1+0.6, 1e20*gmcTBint13/gmc.colDen, \
           width=0.2, bottom=1e-6, \
           color='lightblue', label=r'Milky Way $^{13}$CO')
legend()

# Label transitions
skip=1
ax.set_xticks(arange(0,40,skip)+0.5)
ticklabels=["$J={{{0}}}-{{{1}}}$".format(skip*j+1,skip*j) \
                for j in arange(40)]
ax.set_xticklabels(ticklabels)

# Label axes
ax.set_xlabel('Transition')
ax.set_ylabel(r'$\left[\int\, T_{\rm B}\, dv\right] / N_{\rm H,20}$  [K km s$^{-1}$]')

# Set plot ranges
ax.set_xlim([0,8])
ax.set_ylim([1e-3,10.**1.5])
ax.set_yscale('log')

# Save
savefig('coSLED.eps')
