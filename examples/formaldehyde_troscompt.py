# Import the despotic library
from despotic import cloud
from despotic import emitter

# Read the Milky Way GMC cloud file
gmc = cloud(fileName='cloudfiles/MilkyWayGMC.desp')

import numpy as np
import pylab as pl

gmc.sigmaNT = 1e5 # cm/s, instead of 2 as default
gmc.Tg = 20 # start at 20 instead of 15 K
gmc.Td = 20

# add ortho-h2co
gmc.addEmitter('o-h2co_troscompt', 1e-9)

# first plot: versus density
densities = np.logspace(1,6)
gmc.colDen = 5e21 # use a moderately high column, but not as high as the default

tau11 = np.empty(densities.shape)
tau22 = np.empty(densities.shape)
for ii in xrange(tau11.size):
    gmc.nH = densities[ii]
    line = gmc.lineLum('o-h2co_troscompt')
    tau11[ii] = line[0]['tau']
    tau22[ii] = line[2]['tau']

pl.rc('font',size=20)

pl.figure(1)
pl.clf()
pl.subplots_adjust(hspace=0)
pl.subplot(211)
pl.title("N(H) = %0.2g cm$^{-2}$  T=%0.1f K" % (gmc.colDen,gmc.Tg))
pl.semilogx(densities,tau11,label='1-1',linewidth=4,alpha=0.75)
pl.semilogx(densities,tau22,label='2-2',linewidth=4,alpha=0.75)
#pl.xlabel('Density $n(H)$ cm$^{-3}$')
pl.gca().set_xticks([])
pl.gca().set_yticks(pl.gca().get_yticks()[1:])
pl.ylabel('Optical Depth $\\tau$')
pl.legend(loc='best')
pl.subplot(212)
pl.semilogx(densities,tau11/tau22,linewidth=4,alpha=0.5)
pl.xlabel('Density $n(H)$ cm$^{-3}$')
pl.ylabel('Optical Depth Ratio')
pl.savefig('despotic_formaldehyde_22vs11_density.png',bbox_inches='tight')

# second plot: versus column
gmc.nH = 1e4 # use a high density, right in the interesting zone
columns = np.logspace(20,23)

tau11 = np.empty(columns.shape)
tau22 = np.empty(columns.shape)
for ii in xrange(tau11.size):
    gmc.colDen = columns[ii]
    line = gmc.lineLum('o-h2co_troscompt')
    tau11[ii] = line[0]['tau']
    tau22[ii] = line[2]['tau']

pl.figure(2)
pl.clf()
pl.subplots_adjust(hspace=0)
pl.subplot(211)
pl.title("n(H) = %0.2g cm$^{-3}$  T=%0.1f K" % (gmc.nH,gmc.Tg))
pl.semilogx(columns,tau11,label='1-1',linewidth=4,alpha=0.75)
pl.semilogx(columns,tau22,label='2-2',linewidth=4,alpha=0.75)
#pl.xlabel('Density $n(H)$ cm$^{-3}$')
pl.gca().set_xticks([])
pl.gca().set_yticks(pl.gca().get_yticks()[1:])
pl.ylabel('Optical Depth $\\tau$')
pl.legend(loc='best')
pl.subplot(212)
pl.semilogx(columns,tau11/tau22,linewidth=4,alpha=0.5)
pl.xlabel('Column Density $N(H)$ cm$^{-2}$')
pl.ylabel('Optical Depth Ratio')
pl.savefig('despotic_formaldehyde_22vs11_column.png',bbox_inches='tight')

# third plot: versus temperature
gmc.nH = 1e4 # use a high density, right in the interesting zone
gmc.colDen = 5e21 # use a moderately high column, but not as high as the default
temperatures = np.linspace(5,50)

tau11 = np.empty(temperatures.shape)
tau22 = np.empty(temperatures.shape)
for ii in xrange(tau11.size):
    # assume gas/dust have same temperature
    gmc.Tg = temperatures[ii]
    gmc.Td = temperatures[ii]
    line = gmc.lineLum('o-h2co_troscompt')
    tau11[ii] = line[0]['tau']
    tau22[ii] = line[2]['tau']

pl.figure(3)
pl.clf()
pl.subplots_adjust(hspace=0)
pl.subplot(211)
pl.title("n(H) = %0.2g cm$^{-3}$  N(H) = %0.2g cm$^{-2}$" % (gmc.nH,gmc.colDen))
pl.plot(temperatures,tau11,label='1-1',linewidth=4,alpha=0.75)
pl.plot(temperatures,tau22,label='2-2',linewidth=4,alpha=0.75)
#pl.xlabel('Density $n(H)$ cm$^{-3}$')
pl.gca().set_xticks([])
pl.gca().set_yticks(pl.gca().get_yticks()[1:])
pl.ylabel('Optical Depth $\\tau$')
pl.legend(loc='best')
pl.subplot(212)
pl.plot(temperatures,tau11/tau22,linewidth=4,alpha=0.5)
pl.xlabel('Temperature (K)$')
pl.ylabel('Optical Depth Ratio')
pl.savefig('despotic_formaldehyde_22vs11_temperature.png',bbox_inches='tight')

# fourth plot: versus OPR
gmc.nH = 1e4 # use a high density, right in the interesting zone
gmc.colDen = 5e21 # use a moderately high column, but not as high as the default
gmc.Tg = 20
gmc.Td = 20
oprs = np.logspace(-3,np.log10(3))

tau11 = np.empty(oprs.shape)
tau22 = np.empty(oprs.shape)
for ii in xrange(tau11.size):
    # assume gas/dust have same temperature
    gmc.comp.H2OPR = oprs[ii]
    line = gmc.lineLum('o-h2co_troscompt')
    tau11[ii] = line[0]['tau']
    tau22[ii] = line[2]['tau']

pl.figure(4)
pl.clf()
pl.subplots_adjust(hspace=0)
pl.subplot(211)
pl.title("n(H) = %0.2g cm$^{-3}$  N(H) = %0.2g cm$^{-2}$ T=%0.1f K" % (gmc.nH,gmc.colDen,gmc.Tg))
pl.semilogx(oprs,tau11,label='1-1',linewidth=4,alpha=0.75)
pl.semilogx(oprs,tau22,label='2-2',linewidth=4,alpha=0.75)
#pl.xlabel('Density $n(H)$ cm$^{-3}$')
pl.gca().set_xticks([])
pl.gca().set_yticks(pl.gca().get_yticks()[1:])
pl.ylabel('Optical Depth $\\tau$')
pl.legend(loc='best')
pl.subplot(212)
pl.semilogx(oprs,tau11/tau22,linewidth=4,alpha=0.5)
pl.xlabel('Ortho/Para Ratio')
pl.ylabel('Optical Depth Ratio')
pl.savefig('despotic_formaldehyde_22vs11_orthopararatio.png',bbox_inches='tight')
