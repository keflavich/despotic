########################################################################
# 
# matrixCondition.py
#
# This is an example code for the DESPOTIC library. This code
# demonstrates how DESPOTIC handles ill-conditioned matrices, and
# plots graphical representations of the matrices at various points in
# the procedues DESPOTIC uses to render them calculable.
#
########################################################################

# Import libraries
from despotic import cloud
import numpy as np
import matplotlib.pyplot as plt

########################################################################
# Program code
########################################################################

# Construct a test cloud
cloud = cloud()

# Assign density, gas temperature, abundances; they're all we need for
# this test. Just use pure para-H2 for simplicity.
cloud.nH = 1e3
cloud.comp.xpH2 = 0.5
cloud.Tg = 10.0

# Add three exmaple emitters; abundance values don't matter for this
# example, so just set them to 1
cloud.addEmitter('co', 1.0)
cloud.addEmitter('c+', 1.0, extrap=True, emitterURL='c+@uv.dat')
cloud.addEmitter('o-nh3', 1.0, extrap=True)

# Compute level populations for optically thin cloud with no clumping;
# get back the dict containing diagnostic information
codict = cloud.emitters['co'].setLevPop(cloud, thin=True, noClump=True)
cdict = cloud.emitters['c+'].setLevPop(cloud, thin=True, noClump=True)
nh3dict = cloud.emitters['o-nh3'].setLevPop(cloud, thin=True, noClump=True)

# Print condition numbers, before and after reduction
print "CO condition number = "+str(np.linalg.cond(codict['m']))
print "C condition number = "+str(np.linalg.cond(cdict['m']))
print "NH3 condition number = "+str(np.linalg.cond(nh3dict['m']))

# Plot matrices graphically
mco=codict['m'][:-1,:]
mco += np.identity(mco.shape[0])
mc=cdict['m'][:-1,:]
mc += np.identity(mc.shape[0])
mnh3=nh3dict['m'][:-1,:]
mnh3 += np.identity(mnh3.shape[0])
fig=plt.figure(1, figsize=(7,2.8))
plt.subplot(131)
plt.imshow(np.log10(mco/(np.amax(mco))+1e-20), \
               vmin=-15, vmax=0.0, interpolation='nearest')
plt.xlabel('j')
plt.ylabel('i')
plt.title(r'CO')
plt.subplot(132)
plt.imshow(np.log10(mc/(np.amax(mc))+1e-20), \
               vmin=-15, vmax=0.0, interpolation='nearest')
plt.title(r'C$^+$')
plt.xlabel('j')
plt.subplot(133)
plt.imshow(np.log10(mnh3/(np.amax(mnh3))+1e-20), \
               vmin=-15, vmax=0.0, interpolation='nearest')
plt.title(r'oNH$_3$')
plt.xlabel('j')
plt.subplots_adjust(bottom=0.05, top=0.95, left=0.075, right=0.95)

# Save figure
plt.savefig('matrix1.eps')

# Now repeat for reduced matrices

# Print condition numbers, before and after reduction
print "C condition number (reduced) = "+str(np.linalg.cond(cdict['mRed']))
print "NH3 condition number (reduced) = "+str(np.linalg.cond(nh3dict['mRed']))
print "NH3 condition number (re-reduced) = "+str(np.linalg.cond(nh3dict['mRed2']))

# Graphics
mc=cdict['mRed'][:-1,:]
mc += np.identity(mc.shape[0])
mnh3=nh3dict['mRed2'][:-1,:]
mnh3 += np.identity(mnh3.shape[0])
fig=plt.figure(2, figsize=(5.5,2.8))
plt.subplot(121)
plt.imshow(np.log10(mc/(np.amax(mc))+1e-20), \
               vmin=-15, vmax=0.0, interpolation='nearest')
plt.xlabel('j')
plt.ylabel('i')
plt.title(r'C$^+$')
plt.gca().set_xticks([0,1])
plt.gca().set_yticks([0,1])
plt.subplot(122)
plt.imshow(np.log10(mnh3/(np.amax(mnh3))+1e-20), \
               vmin=-15, vmax=0.0, interpolation='nearest')
plt.xlabel('j')
plt.ylabel('i')
plt.title(r'oNH$_3$')
plt.subplots_adjust(bottom=0.2, top=0.85, left=0.1, right=0.95)

# Save figure
plt.savefig('matrix2.eps')
