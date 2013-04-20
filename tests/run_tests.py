import time

# skip this one execfile('examples/radexComp.py')

# run all examples except RADEX

t0 = time.time()
for fn in 'coSLED','coreTemp','formaldehyde_tests','matrixCondition','pCygni','shockCool']:

    print "Running {0}...".format(fn),
    execfile('examples/{0}.py'.format(fn))
    print "Completed {0} in {1} seconds".format(fn,time.time()-t0)
    print
