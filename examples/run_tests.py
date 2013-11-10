import time
import os
import matplotlib
matplotlib.use('TkAgg')
import pylab

# run all examples except RADEX
# skip this one execfile('examples/radexComp.py')

from despotic import fetchLamda
for species in ('cs','c','co','13co','c18o','hco+','o','cs','oh2co-h2','h2co-h2','o-nh3','p-nh3','oh2o@daniel','ph2o@daniel'):
    fetchLamda(species+".dat")

# keep track of timing
start_time = time.time()
for fn in ['coSLED','coreTemp','formaldehyde_tests','matrixCondition','pCygni','shockCool']:

    this_time = time.time()
    print "Running {0}... started at {1}".format(fn,start_time)
    scriptfile = '{0}.py'.format(fn)
    if os.path.exists('tests/'+scriptfile):
        execfile('tests/'+scriptfile)
    elif os.path.exists('examples/'+scriptfile):
        execfile('examples/'+scriptfile)
    else:
        raise IOError("Missing script file for {0}".format(str(fn)))
    print "Completed {0} in {1} seconds.  Elapsed time is {2} seconds.".format(fn,
            time.time()-this_time, time.time()-start_time)
    print
