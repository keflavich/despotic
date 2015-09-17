########################################################################
# 
# radexComp.py
#
# This is an example code for the DESPOTIC library. This code
# compares DESPOTIC results to the results of the RADEX code of van
# der Tak et al. (2007). Note that, for this code to work properly,
# the RADEX code must be downloaded and compiled separately, and the
# executable radex must be present in the directory where this code is
# being run.
#
########################################################################

# Import libraries and constants
from despotic import cloud
from numpy import *
from datetime import datetime
from datetime import timedelta
import os
import scipy.constants as physcons
from copy import deepcopy
kB = physcons.k/physcons.erg
mH = physcons.m_p/physcons.gram

########################################################################
# User-settable options
########################################################################

# Set up a range of densities and column densities
lognHgrid = arange(2, 8.01, 0.2)
logcolDengrid = arange(14, 24.01, 0.2)

# Set total FWHM velocity to 2 km/s
FWHM = 2.0e5

########################################################################
# Program code
########################################################################

# Construct a test cloud
cloud = cloud()

# Assign temperature; we'll keep this fixed as we vary the volume and
# column density
cloud.Tg = 10.0

# Set dust opacity to 0, since RADEX doesn't include it
cloud.dust.sigma10 = 0.0

# Add CO as an emitter
cloud.addEmitter('co', 1e-4)

# Assign abundances of bulk constituents. To do a fair comparison, we
# make the same assumptions RADEX does about abundances, i.e. that the
# ratio of ortho-to-para H_2 equals the thermal ratio of the J = 0 to
# J = 1 level of H2. The Boltzmann factor of J = 1 state is 3(2J+1)
# exp(-J(J+1)Theta_rot / T), where Theta_rot = 85.3 K
boltzfac_oH2 = 9.0*exp(-2*85.3/cloud.Tg)
cloud.comp.xoH2 = 0.5 * boltzfac_oH2 / (1.0 + boltzfac_oH2)
cloud.comp.xpH2 = 0.5 / (1.0 + boltzfac_oH2)

# Compute sound speed
cloud.nH = 1.0  # Avoid annoying warning message
cloud.comp.computeDerived(cloud.nH)
cs = sqrt(kB*cloud.Tg/(cloud.comp.mu*mH))

# Compute non-thermal velocity dispersion to get desired FWHM
sigmaTot = FWHM/sqrt(8.0*log(2))
cloud.sigmaNT = sqrt(sigmaTot**2 - \
                         cs**2/cloud.emitters['co'].data.molWgt)

# Initialize execution time counters
timeDesp=timedelta(0)
timeRadex=timedelta(0)

# Open radex input file
fp = open('despotic_comp.inp', 'w')

# Loop over the volume and column density grids, writing to radex
# input file
for i, lognH in enumerate(lognHgrid):
    for j, logcolDen in enumerate(logcolDengrid):
        fp.write('LAMDA/co.dat\n')
        fp.write('despotic_comp.rdx\n')
        fp.write('0 0\n')
        fp.write('10.0\n')
        fp.write('1\n')
        fp.write('H2\n')
        fp.write(str(10.**lognH/2.0)+'\n')
        fp.write('2.73\n')
        fp.write(str(10.0**logcolDen*cloud.emitters['co'].abundance)+'\n')
        fp.write(str(FWHM/1e5)+'\n')
        if i==len(lognHgrid)-1 and j==len(logcolDengrid)-1:
            fp.write('0\n')
        else:
            fp.write('1\n')        
fp.close()

# Run radex, timing execution
t1 = datetime.now()
os.system('./radex < despotic_comp.inp > /dev/null')
t2 = datetime.now()
timeRadex = t2-t1
print "RADEX time = "+str(timeRadex)

# Create lists to save despotic results
cloudList = []
lineList = []

# Loop over volume and column density grids
for lognH in lognHgrid:

    # Assign density to cloud
    cloud.nH = 10.**lognH

    for logcolDen in logcolDengrid:

        # Assign column density to cloud
        cloud.colDen = 10.**logcolDen

        # Compute lines with despotic
        t1 = datetime.now()
        lines=cloud.lineLum('co', noClump=True, lumOnly=True, \
                                escapeProbGeom='slab', abstol=1e-4)
        t2 = datetime.now()
        timeDesp += t2 - t1

        # Save cloud
        cloudList.append(deepcopy(cloud))
        lineList.append(deepcopy(lines))

# Print time and ratio
print "DESPOTIC time = "+str(timeDesp)
print "DESPOTIC / RADEX = " + \
    str(timeDesp.total_seconds()/timeRadex.total_seconds())

# Prepare to read radex results
fp = open('despotic_comp.rdx', 'r')

# Arrays to hold radex results
levPop = zeros(40)
lineFlux = zeros(40)

# Loop over grid and compare results
ctr=0
maxLevPopDiff = zeros((len(lognHgrid), len(logcolDengrid)))
for i, lognH in enumerate(lognHgrid):
    for j, logcolDen in enumerate(logcolDengrid):

        # Skip 13 header lines
        for n in range(13):
            fp.readline()

        # Read 40 lines of data
        for n in range(40):
            line = fp.readline()
            split = line.split()
            levPop[n] = float(split[10])
            lineFlux[n] = float(split[12])

        # Compare to despotic result
        cloud = cloudList[ctr]
        lines = lineList[ctr]
        maxLevPopDiff[i,j] = \
            amax(abs(levPop - cloud.emitters['co'].levPop[:-1]))

        # Increment counter
        ctr = ctr+1

# Print difference
print 'Maximum level population difference = '+str(amax(maxLevPopDiff))
