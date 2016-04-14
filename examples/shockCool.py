########################################################################
# 
# shockCool.py
#
# This is an example code for the DESPOTIC library. This code
# considers a portion of a molecular cloud that has been shock-heated
# to 250 K, and then calculates its time-dependent cooling back to its
# equilibrium temperature. To illustrate the capability of pickling
# cloud states, it saves intermediate states of the calculation as
# python pickles.
#
########################################################################

########################################################################
# User-settable options
########################################################################

# Set total time for which to evolve cloud, in sec
kyr = 365.25*24.*3600.*1e3
tEvol = 40*kyr

# Set desired number of outputs
nOut = 41    # output once per kyr

# Specify whether verbose printing while running is desired
verbose = True

# Damping factor for calculation of level populations; normally
# despotic's default of 0.5 is a good choice, but in this case
# choosing something a bit smaller makes things run ~10% faster,
# though the code will still run successfully even for the default
# choice.
dampFactor = 0.25


########################################################################
# Program code
########################################################################

# Import the despotic library and various standard python libraries
from despotic import cloud
from numpy import *
import matplotlib.pyplot as plt
import pickle
from copy import deepcopy
from datetime import datetime
from datetime import timedelta


# Create array of times
times = arange(0, tEvol*(1.0+1.0e-6), tEvol/(nOut-1))

# Create list of states
stateList = []

# Read any existing data
restart=False
for i in arange(0, nOut, 1):
    try:
        fp = open('shockCool{:03d}.pkl'.format(i), 'rb')
        stateList.append(pickle.load(fp))
        fp.close()
        restart=True
    except IOError:
        break

# Did we find any existing data?
if restart:
    # Yes: copy last state to slab object to initialize it
    slab = deepcopy(stateList[-1])
else:
    # No
    # Read the cloud initialization file
    slab=cloud(fileName='cloudfiles/postShockSlab.desp', \
                   verbose=verbose)
    # Compute initial level populations and dust temperature
    slab.dEdt(overrideSkip=True, dampFactor=dampFactor, \
                  verbose=verbose)
    slab.setDustTempEq(verbose=verbose)
    # Save initial state, both locally and to disk
    stateList.append(deepcopy(slab))
    fp = open('shockCool000.pkl', 'wb')
    pickle.dump(slab, fp)
    fp.close()

# Calculate time evolution, starting at proper position
istart = len(stateList) - 1
extime = timedelta(0)
for i, t in enumerate(times[istart:-1]):

    # Evolve in time
    t1=datetime.now()
    slab.tempEvol(times[i+istart+1], tInit=t, \
                      escapeProbGeom='slab', isobaric=True, \
                      verbose=verbose, nOut=1, \
                      dampFactor=dampFactor)

    # Compute cooling rates from all coolants, so that we will have
    # these stored on disk
    slab.dEdt(overrideSkip=True, verbose=True)
    t2=datetime.now()
    extime += t2-t1

    # Save state, both locally and to disk
    stateList.append(deepcopy(slab))
    fp = open('shockCool{:03d}.pkl'.format(i+istart+1), 'wb')
    pickle.dump(slab, fp)
    fp.close()

# Print execution time
print("Execution time = "+str(extime))

# Get temperature history from saved states
Tg=array([s.Tg for s in stateList])
Td=array([s.Td for s in stateList])
nH=array([s.nH for s in stateList])

# Compute cooling rates from saved states
rates=[s.dEdt(fixedLevPop=True) for s in stateList]
PsiGD = array([r['PsiGD'] for r in rates])
PsiGD[PsiGD > 0.0] = 0.0
PsiGD = -PsiGD
LambdaSpec = {}
for k in rates[0]['LambdaLine'].keys():
    LambdaSpec[k] =  array([r['LambdaLine'][k] for r in rates])
LambdaLine = array([sum(array(list(r['LambdaLine'].values()))) for r in rates])

# Compute luminosities of individual CO lines
coLines = zeros((stateList[0].emitters['co'].data.nlev-1, nOut))
for i, s in enumerate(stateList):
    coLines[:,i] = s.lineLum('co', lumOnly=True, noRecompute=True)

# First plot: temperature, density versus time
fig, ax = plt.subplots(1, figsize=(6,4))
plt.subplots_adjust(bottom=0.15, right=0.85, left=0.15)
axes = [ax, ax.twinx()]
ln1 = axes[0].plot(times/kyr, Tg, 'r', linewidth=2, label=r'$T_g$')
ln2 = axes[0].plot(times/kyr, Td, 'g', linewidth=2, label=r'$T_d$')
axes[0].set_yscale('log')
axes[0].set_ylabel('T [K]')
axes[0].set_xlabel('t [kyr]')
axes[0].set_ylim([8,800])
ln3 = axes[1].plot(times/kyr, nH, 'b', linewidth=2, label=r'$n_{\rm H}$')
axes[1].set_ylabel(r'$n_{\rm H}$ [cm$^{-3}$]')
axes[1].set_yscale('log')
axes[1].set_ylim([1e3,1e5])
lns = ln1 + ln2 + ln3
labs = [l.get_label() for l in lns]
axes[0].legend(lns, labs, loc='upper right', ncol=3)
plt.savefig('shockCool_temp.eps')


# Second plot: contributions of various coolants
fig2=plt.figure(2, figsize=(6,4))
plt.subplots_adjust(bottom=0.15, left=0.15)
plt.plot(times/kyr, PsiGD+LambdaLine, color='0.5', linewidth=10, label='All')
plt.plot(times/kyr, LambdaSpec['co'], linewidth=2, label='CO')
plt.plot(times/kyr, LambdaSpec['13co'], linewidth=2, label=r'$^{13}$CO')
plt.plot(times/kyr, LambdaSpec['o'], linewidth=2, label='O')
plt.plot(times/kyr, PsiGD, linewidth=2, label='Dust')
plt.xlabel('t [kyr]')
plt.ylabel(r'$\Lambda$ [erg s$^{-1}$ H$^{-1}$]')
plt.yscale('log')
plt.ylim([1e-28,1e-24])
plt.legend()
plt.savefig('shockCool_coolants.eps')

# Third plot: CO SLED versus time
fig3=plt.figure(3, figsize=(6,4))
plt.subplots_adjust(bottom=0.15)
ax3=fig3.add_subplot(111)
plt.plot(coLines[:,0]/LambdaSpec['co'][0], linewidth=2, label='0 kyr')
plt.plot(coLines[:,20]/LambdaSpec['co'][20], linewidth=2, label='20 kyr')
plt.plot(coLines[:,40]/LambdaSpec['co'][40], linewidth=2, label='40 kyr')
plt.legend()
plt.ylim([2e-3,1])
plt.yscale('log')
plt.ylabel('CO luminosity fraction')
plt.xlabel('CO transition')

# Label transitions
skip=2
ax3.set_xticks(arange(0,40,skip))
ticklabels=["$J={{{0}}}-{{{1}}}$".format(skip*j+1,skip*j) \
                for j in arange(40)]
ax3.set_xticklabels(ticklabels)
plt.xlim([0,12])

# Save
plt.savefig('shockCool_coSLED.eps')
