"""
This is an example code for the DESPOTIC library. This code
calculates the chemical state and equilibrium temperature, and the
major heating and cooling processes, for molecular clouds as a
function of gas density.
"""

# Import the despotic library and various standard python libraries
from despotic import cloud
from despotic.chemistry import NL99
from numpy import *
import matplotlib.pyplot as plt

########################################################################
# User-settable options
########################################################################

# Set up a range of volume densities
lognHgrid = arange(0.9,4.01,0.05)

########################################################################
# Program code
########################################################################

# Constants
import scipy.constants as physcons
kB = physcons.k/physcons.erg
mH = physcons.m_p/physcons.gram
G = physcons.G*1e3

# Only recompute if we haven't already done the computation
try:
    gmc
except NameError:

    # Use the Milky Way GMC file as a base, but add C+ and O as emitting species
    gmc=cloud('../cloudfiles/MilkyWayGMC.desp')
    gmc.addEmitter('c+', 1e-10)
    gmc.addEmitter('o', 1e-4)

    # Set CR ionization rate
    gmc.rad.ionRate = 3e-17

    # Set IR temp to 10 K
    gmc.rad.TradDust = 10.0
    gmc.Td = 10.0

    # Set column density to 10^22 cm^-2
    gmc.colDen = 1.0e22

    # Set initial temperature guess to 20 K
    gmc.Tg = 20.0

    # Turn on extrapolation
    gmc.emitters['co'].extrap = True
    gmc.emitters['c+'].extrap = True
    gmc.emitters['13co'].extrap = True
    gmc.emitters['o'].extrap = True

    # Loop over grid, calculating composition and gas and dust temperatures
    abd=[]
    rates=[]
    Tg=[]
    Td=[]
    tcool=[]
    for i, logn in enumerate(lognHgrid):

        # Set density
        gmc.nH = 10.0**logn

        # Calculate chemical abundance at T = 10 K
        gmc.setChemEq(network=NL99, verbose=True)
        abd.append(gmc.chemnetwork.abundances)

        # Recalculate derived quantities (this changes the CR energy
        # yield)
        gmc.comp.computeDerived(gmc.nH)

        # Set virial ratio to unity
        gmc.setVirial()

        # Print status message
        print "Calculating with lognH = "+str(logn)+"..."

        # Calculate equilibrium temperature
        gmc.setTempEq(verbose=True)

        # Print status message
        print "Converged to Tg = "+str(gmc.Tg)+" K, Td = "+str(gmc.Td)+" K"
        print ""

        # Record result; note that we use the fixedLevPop option for the
        # rates method, so that the level populations are not
        # unnecessarily recomputed.
        Tg.append(gmc.Tg)
        Td.append(gmc.Td)
        rates.append(gmc.dEdt(fixedLevPop=True))
        tcool.append(gmc.comp.computeEint(gmc.Tg)*kB*gmc.Tg / 
                     (rates[-1]['GammaCR']+rates[-1]['GammaPE']))

    # Construct arrays of various terms for convenience of plotting
    GammaCR = array([r['GammaCR'] for r in rates])
    GammaPE = array([r['GammaPE'] for r in rates])
    PsiGD = array([r['PsiGD'] for r in rates])
    GammaDustIR = array([r['GammaDustIR'] for r in rates])
    GammaDustISRF = array([r['GammaDustISRF'] for r in rates])
    GammaDustLine = array([r['GammaDustLine'] for r in rates])
    LambdaDust = array([r['LambdaDust'] for r in rates])
    LambdaCO = array([r['LambdaLine']['co'] for r in rates])
    LambdaCp = array([r['LambdaLine']['c+'] for r in rates])
    Lambda13CO = array([r['LambdaLine']['13co'] for r in rates])
    LambdaLine = array([sum(array(r['LambdaLine'].values())) for r in rates])
    xCO = array([ab['CO'] for ab in abd])
    xCp = array([ab['C+'] for ab in abd])
    xC = array([ab['C'] for ab in abd])

# Plots
plt.figure(1, figsize=(5,10))

# Gas temperature
ax1 = plt.subplot(411)
ax1.plot(lognHgrid, Tg, 'g', linewidth=2, label='T')
ax1.set_xlim([1,4])
ax1.set_ylim([0,20])
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel('T [K]')

# Abundances
ax2=plt.subplot(412)
ax2.plot(lognHgrid, xCO/2e-4, 'r--', linewidth=2, label=r'$f_{\rm CO}$')
ax2.plot(lognHgrid, xCp/2e-4, 'b--', linewidth=2, label=r'$f_{\rm C^+}$')
ax2.set_ylabel(r'$f_{\rm C^+}, f_{\rm CO}$')
ax2.set_xlim([1,4])
ax2.set_ylim([0,1.79])
plt.setp(ax2.get_xticklabels(), visible=False)
leg2=ax2.legend(loc=2)

# Heating and cooling terms
ax3=plt.subplot(413)
p1,=ax3.plot(lognHgrid, log10(GammaCR), 'b', linewidth=2, label='CR')
p2,=ax3.plot(lognHgrid, log10(GammaPE), 'b--', linewidth=2, label='PE')
p3,=ax3.plot(lognHgrid, log10(LambdaCO), 'r', linewidth=2, label='CO')
p4,=ax3.plot(lognHgrid, log10(LambdaCp), 'r--', linewidth=2, label=r'C$^+$')
p5,=ax3.plot(lognHgrid, log10(PsiGD), 'b:', linewidth=2, label='Dust')
ax3.set_xlim([1,4])
ax3.set_ylim([-29,-25.1])
leg3a=ax3.legend([p1,p2,p5], ['CR', 'PE', 'Dust'], loc=2)
leg3b=ax3.legend([p3,p4], ['CO', r'C$^+$'], loc=1)
plt.gca().add_artist(leg3a)
ax3.set_ylabel(r'$\Gamma, \Lambda$ [erg s$^{-1}$]')
plt.setp(ax3.get_xticklabels(), visible=False)

# Cooling time
tff=sqrt(3*pi/(32*G*10.**lognHgrid*gmc.comp.muH*mH))
ax4=plt.subplot(414)
p1,=ax4.plot(lognHgrid, log10(tcool/tff), linewidth=2)
ax4.set_ylim([-2,-0.01])
ax4.set_xlim([1,4])
ax4.set_xlabel(r'$\log\,n_{\rm H}$ [cm$^{-3}$]')
ax4.set_ylabel(r'$\log\,t_{\rm eq}/t_{\rm ff}$')

plt.subplots_adjust(hspace=0, top=0.97, bottom=0.05, left=0.18, right=0.95)

plt.show()
plt.savefig('despotic_sfreview.pdf')
