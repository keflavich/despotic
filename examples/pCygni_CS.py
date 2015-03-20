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

# Read the data files for CS and C^34S
cs=emitterData('cs')
c34s=emitterData('c34s')

# Core radius; set to 0.02 pc
R = 0.02*3.09e18

# Number densities of CS and C^34S, in cm^-3; the abundance ratio is
# set to the terrestrial value
ncs = 1e-1
nc34s = ncs / 22.0

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

# Calculate the line profile for CS(3-2). The functon returns two
# arrays, the first giving the brightness temperature and the second
# giving the velocity at which that brightness temperature is
# measured. The arguments of the function are, in order, the molecule
# to use, the upper state of the transition (where states are ordered
# by energy and the ground state is state 0), the lower state of the
# transition, the density, the temperature, the radial velocity, and
# the non-thermal velocity dispersion. The last four of this can be
# either constants or functions, and additional optional arguments
# also exist. See the User's Guide.
TBcs, vOut = lineProfLTE(cs, 2, 1, R, ncs, TProf, vProf, sigmaNT, \
                             vLim=[-2e5,2e5], nOut=400)

# Perform the same calculation for C^34S(3-2)
TBc34s, vOut1 = lineProfLTE(c34s, 2, 1, R, nc34s, TProf, vProf,
                            sigmaNT, vLim=[-2e5,2e5], nOut=400)

# Now make plots of both lines; divide by 10^5 to convert cm/s to km/s
plt.figure(1, figsize=(6,4))
plt.plot(vOut/1e5, TBcs, label='CS(3-2)', linewidth=2)
plt.plot(vOut1/1e5, TBc34s, label='C$\,^{34}$S(3-2)', linewidth=2)
plt.subplots_adjust(bottom=0.15)

# Adjust range and add labels
plt.xlim([-2,2])
plt.xlabel('v [km s$^{-1}$]')
plt.ylabel('$T_B$ [K]')

# Add legend
plt.legend()

# Save
plt.savefig('pCygni.eps')
