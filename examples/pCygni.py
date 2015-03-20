########################################################################
# 
# pCygni.py
#
# This is an example code for the DESPOTIC library. This code
# calculates the line profiles through a collapsing protostellar core
# in the lines CS(3-2) and C^34S(3-2), under conditions where the
# former is marginally opticall thick and shows an inverse pCygi
# profile, while the latter is optically thin and does not.
#
########################################################################

########################################################################
# Program code
########################################################################

# Import general numpy and plotting stuff
from numpy import *
import matplotlib.pyplot as plt

# Import the line profile function and the emitter data class from
# DESPOTIC
from despotic.lineProfLTE import lineProfLTE
from despotic.emitterData import emitterData

# Read the data files for HCN and N2H+
HCN=emitterData('hcn')
N2Hp=emitterData('n2h+')

# Core radius; set to 0.02 pc
R = 0.02*3.09e18

# Abundances and absolute number densities of HCN and N2H+, in cm^-3;
# abundances are taken from the models of Lee et al. (2004, ApJ, 617,
# 360).
xHCN = 2e-8
xN2Hp = 2e-9
nH = 3e6
nHCN = nH*xHCN
nN2Hp = nH*xN2Hp

# Define a temperature profile so that we can get a nice p Cygni
# profile. The choice here is somewhat arbitrary. For simplicity we'll
# just take 8 K at large radii, added to a Gaussian component rising
# to 20 K at small radii. This function takes as an argument the
# radius r normalized to the core radius, and returns the temperature
# at that radius.
def TProf(r):
    return 8.0 + 12.0*exp(-r**2/(2.0*0.5**2))

# Choose a velocity versus radius. We'll use a velocity simple linear
# profile v ~ r, rising to a maximum amplitude of -0.4 km/s at the
# cloud edge.
def vProf(r):
    return -4.0e4*r

# Choose a non-thermal velocity dispersion. We'll also make this
# subsonic.
sigmaNT = 2.0e4

# Calculate the line profile for HCN(1-0). The functon returns two
# arrays, the first giving the brightness temperature and the second
# giving the velocity at which that brightness temperature is
# measured. The arguments of the function are, in order, the molecule
# to use, the upper state of the transition (where states are ordered
# by energy and the ground state is state 0), the lower state of the
# transition, the density, the temperature, the radial velocity, and
# the non-thermal velocity dispersion. The last four of this can be
# either constants or functions, and additional optional arguments
# also exist. See the User's Guide.
TB_HCN, vOut = lineProfLTE(HCN, 1, 0, R, nHCN, TProf, vProf, sigmaNT, \
                           vLim=[-2e5,2e5], nOut=400)

# Perform the same calculation for N2H+ (1-0)
TB_N2Hp, vOut1 = lineProfLTE(N2Hp, 1, 0, R, nN2Hp, TProf, vProf,
                             sigmaNT, vLim=[-2e5,2e5], nOut=400)

# Now make plots of both lines; divide by 10^5 to convert cm/s to km/s
plt.figure(1, figsize=(6,4))
plt.plot(vOut/1e5, TB_HCN, label='HCN(1-0)', linewidth=2)
plt.plot(vOut1/1e5, TB_N2Hp, label='N$_2$H$^+$(1-0)', linewidth=2)
plt.subplots_adjust(bottom=0.15)

# Adjust range and add labels
plt.xlim([-2,2])
plt.xlabel('v [km s$^{-1}$]')
plt.ylabel('$T_B$ [K]')

# Add legend
plt.legend()

# Save
plt.savefig('pCygni.eps')
