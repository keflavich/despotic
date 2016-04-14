"""
This is an example code for the DESPOTIC library. This code
calculates the equilibrium temperature and the major coolants for
protostellar cores over a grid in density from 10^2 - 10^8
cm^-3. It then produces a plot of the results. 

This code saves its work as it runs, and before starting it searches
for saved work in the directory where it is running, in the files
coreTemp_Tg.pkl, coreTemp_Td.pkl, and coreTemp_rates.pkl. If the
user wants, he/she can prevent this behavior just by deleting or
renaming those files. If this code resumes from saved work, the
results will not be bitwise identical to the results produced if the
code is run from the start, due to differing starting guesses in the
iterative solvers. However, the results will agree to within the
tolerances of the solvers.
"""

# Import the despotic library and various standard python libraries
from despotic import cloud
import numpy as np
import matplotlib.pyplot as plt
import pickle
from datetime import datetime
from datetime import timedelta

########################################################################
# User-settable options
########################################################################

# Set up a range of densities
lognHgrid = np.arange(2,6.01,0.2)

# Specify whether verbose printing while running is desired
verbose = True


########################################################################
# Program code
########################################################################

# Read the protostellar core file
core = cloud(fileName='cloudfiles/protostellarCore.desp', verbose=True)

# Check if we have saved work
try:
    # Load pickle files
    inFile = open('coreTemp_Tg.pkl', 'rb')
    Tg = pickle.load(inFile)
    inFile.close()
    inFile = open('coreTemp_Td.pkl', 'rb')
    Td = pickle.load(inFile)
    inFile.close()
    inFile = open('coreTemp_rates.pkl', 'rb')
    rates = pickle.load(inFile)
    inFile.close()
    startIdx = len(Tg)    # Point at which to restart
except IOError:
    # If we're here, no saved work was found, so create empty lists
    # and set the starting index to 0
    Tg = []
    Td = []
    rates = []
    startIdx = 0

# Loop over grid, calculating gas and dust temperatures
extime=timedelta(0)
for i, logn in enumerate(lognHgrid[startIdx:]):

    # Set density
    core.nH = 10.0**logn

    # Recalculate derived quantities (this changes the CR energy
    # yield)
    core.comp.computeDerived(core.nH)

    # Print status message
    print("Calculating core with lognH = "+str(logn)+"...")

    # Calculate equilibrium temperature
    t1=datetime.now()
    core.setTempEq(c1Grav=1.0, verbose=verbose)
    t2=datetime.now()
    extime += t2-t1

    # Print status message
    print("Converged to Tg = "+str(core.Tg)+" K, Td = "+str(core.Td)+" K")
    print("")

    # Record result; note that we use the fixedLevPop option for the
    # rates method, so that the level populations are not
    # unnecessarily recomputed.
    Tg.append(core.Tg)
    Td.append(core.Td)
    rates.append(core.dEdt(c1Grav=1.0, fixedLevPop=True))

    # Write intermediate work to disk
    outFile = open('coreTemp_Tg.pkl', 'wb')
    pickle.dump(Tg, outFile)
    outFile.close()
    outFile = open('coreTemp_Td.pkl', 'wb')
    pickle.dump(Td, outFile)
    outFile.close()
    outFile = open('coreTemp_rates.pkl', 'wb')
    pickle.dump(rates, outFile)
    outFile.close()


# Print execution time
print("Execution time = "+str(extime))

# Construct arrays of various cooling terms from the rates dict, for
# convenience of plotting
GammaCR = np.array([r['GammaCR'] for r in rates])
PsiGD = np.array([r['PsiGD'] for r in rates])
GammaDustIR = np.array([r['GammaDustIR'] for r in rates])
GammaDustISRF = np.array([r['GammaDustISRF'] for r in rates])
GammaDustLine = np.array([r['GammaDustLine'] for r in rates])
LambdaDust = np.array([r['LambdaDust'] for r in rates])
LambdaCO = np.array([r['LambdaLine']['co'] for r in rates])
LambdaC = np.array([r['LambdaLine']['c'] for r in rates])
Lambda13CO = np.array([r['LambdaLine']['13co'] for r in rates])
LambdaC18O = np.array([r['LambdaLine']['c18o'] for r in rates])
LambdaHCOp = np.array([r['LambdaLine']['hco+'] for r in rates])
LambdaCS = np.array([r['LambdaLine']['cs'] for r in rates])
LambdaO = np.array([r['LambdaLine']['o'] for r in rates])
LambdaoH2CO = np.array([r['LambdaLine']['o-h2co'] for r in rates])
LambdapH2CO = np.array([r['LambdaLine']['p-h2co'] for r in rates])
LambdaoNH3 = np.array([r['LambdaLine']['o-nh3'] for r in rates])
LambdapNH3 = np.array([r['LambdaLine']['p-nh3'] for r in rates])
LambdaoH2O = np.array([r['LambdaLine']['oH2O'] for r in rates])
LambdapH2O = np.array([r['LambdaLine']['pH2O'] for r in rates])
LambdaLine = np.array([np.sum(np.array(list(r['LambdaLine'].values())))
                       for r in rates])


# First plot: gas and dust temperature versus density
plt.figure(1, figsize=(4,3))
plt.subplot(111)
plt.plot(lognHgrid, Tg, linewidth=2, label='Gas')
plt.plot(lognHgrid, Td, linewidth=2, label='Dust')
plt.xlim([2,6])
plt.ylim([0,30])
plt.legend()
plt.xlabel(r'$\log\,n_{\rm H}$ [cm$^{-3}$]')
plt.ylabel('T [K]')
plt.subplots_adjust(bottom=0.2, left=0.15)
plt.savefig('coreTemp.eps')


# Second plot: heating and cooling processes. We will have 4 frames,
# showing gas heating terms, gas cooling terms, dust heating terms,
# and dust cooling terms
plt.figure(2, figsize=(8,6))
plt.subplots_adjust(wspace=0, hspace=0)
xlim=[2,6]
xlim1=[2,5.999]
ylim=[-31,-21.5]
ylim1=[-31,-21.5]

# Gas heating
plt.subplot(221)
plt.plot(lognHgrid, np.log10(GammaCR), linewidth=2, label='Ionization')
leg=plt.legend(title='Gas heating')
leg.get_title().set_fontsize('14')
plt.xlim(xlim1)
plt.ylim(ylim)
plt.ylabel(r'$\log\,\Gamma$ [erg s$^{-1}$ H$^{-1}$]')
plt.setp(plt.gca().get_xticklabels(), visible=False)

# Dust heating
plt.subplot(222)
plt.plot(lognHgrid, np.log10(GammaDustIR), linewidth=2, label='IR')
plt.plot(lognHgrid, np.log10(GammaDustISRF), linewidth=2, label='ISRF')
plt.plot(lognHgrid, np.log10(GammaDustLine), linewidth=2, label='Line')
plt.plot(lognHgrid, np.log10(-PsiGD), linewidth=2, label='Gas-Dust')
leg=plt.legend(ncol=2, title='Dust heating')
leg.get_title().set_fontsize('14')
plt.xlim(xlim)
plt.ylim(ylim)
plt.setp(plt.gca().get_xticklabels(), visible=False)
plt.setp(plt.gca().get_yticklabels(), visible=False)

# Gas cooling
plt.subplot(223)
plt.plot(lognHgrid, np.log10(LambdaLine), linewidth=2, label='Lines (sum)')
plt.plot(lognHgrid, np.log10(-PsiGD), linewidth=2, label='Gas-Dust')
leg=plt.legend(title='Gas cooling')
leg.get_title().set_fontsize('14')
plt.xlim(xlim1)
plt.ylim(ylim1)
plt.ylabel(r'$\log\,\Lambda$ [erg s$^{-1}$ H$^{-1}$]')
plt.xlabel(r'$\log\,n_{\rm H}$ [cm$^{-3}$]')

# Dust cooling
plt.subplot(224)
plt.plot(lognHgrid, np.log10(LambdaDust), linewidth=2, label='Thermal')
leg=plt.legend(title='Dust cooling')
leg.get_title().set_fontsize('14')
plt.xlim(xlim)
plt.ylim(ylim1)
plt.xlabel(r'$\log\,n_{\rm H}$ [cm$^{-3}$]')
plt.setp(plt.gca().get_yticklabels(), visible=False)

# Save figure
plt.savefig('coreHeatCool.eps')


# Third plot: contributions of individual lines to cooling
plt.figure(3, figsize=(8,6))
plt.plot(lognHgrid, np.log10(LambdaLine), 'k', linewidth=8, label='Lines (sum)')
plt.plot(lognHgrid, np.log10(LambdaCO), linewidth=2, label='CO')
plt.plot(lognHgrid, np.log10(Lambda13CO), linewidth=2, label=r'$^{13}$CO')
plt.plot(lognHgrid, np.log10(LambdaC18O), linewidth=2, label=r'C$^{18}$O')
plt.plot(lognHgrid, np.log10(LambdaC), linewidth=2, label='C')
plt.plot(lognHgrid, np.log10(LambdaO), linewidth=2, label='O')
plt.plot(lognHgrid, np.log10(LambdaCS), linewidth=2, label='CS')
plt.plot(lognHgrid, np.log10(LambdaHCOp), '--', linewidth=2, label=r'HCO$^+$')
plt.plot(lognHgrid, np.log10(LambdaoH2CO), '--', linewidth=2, label=r'oH$_2$CO')
plt.plot(lognHgrid, np.log10(LambdapH2CO), '--', linewidth=2, label=r'pH$_2$CO')
plt.plot(lognHgrid, np.log10(LambdaoNH3), '--', linewidth=2, label=r'oNH$_3$')
plt.plot(lognHgrid, np.log10(LambdapNH3), '--', linewidth=2, label=r'pNH$_3$')
plt.plot(lognHgrid, np.log10(LambdaoH2O), '--', linewidth=2, label=r'oH$_2$O')
plt.plot(lognHgrid, np.log10(LambdapH2O), '--', linewidth=2, label=r'pH$_2$O')
plt.xlim([2,6])
plt.ylim([-33, -27])
plt.ylabel(r'$\log\,\Lambda$ [erg s$^{-1}$ H$^{-1}$]')
plt.xlabel(r'$\log\,n_{\rm H}$ [cm$^{-3}$]')
plt.legend(ncol=3, loc='lower right')
plt.savefig('coreLines.eps')
