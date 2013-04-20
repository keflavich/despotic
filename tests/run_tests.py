import time

# run all examples except RADEX
# skip this one execfile('examples/radexComp.py')

# keep track of timing
start_time = time.time()
for fn in ['coSLED','coreTemp','formaldehyde_tests','matrixCondition','pCygni','shockCool']:

    this_time = time.time()
    print "Running {0}... started at {1}".format(fn,start_time)
    execfile('examples/{0}.py'.format(fn))
    print "Completed {0} in {1} seconds.  Elapsed time is {2} seconds.".format(fn,
            time.time()-this_time, time.time()-start_time)
    print
