"""
Script to generate a hot wind output table
"""

import subprocess
import shlex
import numpy as np
import time
from multiprocessing import cpu_count
import os.path as osp
from shutil import copyfileobj

# Set to overwrite existing grid
overwrite = False

# Set directory into which to write output
tmpdir = "hot_gas_data_src"
outdir = "hot_gas_data"

# Resolution of grid to compute
nu = 1024
nq = 256
ngex = 512
qmin = -4.0
qmax = 4.0
loggex_max = 3.0
ngex_lo = 64

# Range of parameters
m = np.array([0, 1], dtype=int)
y = np.array([0, 1, 2], dtype=int)
uh =  np.append(
    np.append(np.arange(1.0, 10.01, 0.1), np.arange(10.5, 20.01, 0.5)),
    np.arange(21.0, 50.01, 1.0))

# Construct array of values to loop over
ma, ya, uha = np.meshgrid(m, y, uh)
ma = ma.flatten()
ya = ya.flatten()
uha = uha.flatten()

# Calculate all cases in parallal
ncpus = cpu_count()
ptr = 0
procs = []
logfiles = []
while ptr < ya.size:

    # Check if there is a free slot in the process list; if there is,
    # start a new job
    while len(procs) < ncpus:

        # Build command string for this case
        cmd = (
            "./hot_wind_tab {:f} {:d} {:d} {:d} {:d} {:f} {:f} {:d} "
            "{:f} {:d} --dir {:s} --verbose").format(
                uha[ptr], ya[ptr], ma[ptr],
                nu, nq, qmin, qmax, ngex, loggex_max, ngex_lo,
                tmpdir)
        if overwrite:
            cmd += " --overwrite"
        args = shlex.split(cmd)

        # Set name of log file
        logfile = "uh{:04.1f}_y{:d}_m{:d}.log".\
                  format(uha[ptr], ya[ptr], ma[ptr])
        logfiles.append(logfile)

        # Start process
        print(cmd)
        procs.append(subprocess.Popen(args, stdout=subprocess.PIPE))

        # Increment pointer; if we are out of jobs to do, stop
        ptr += 1
        if ptr >= ya.size: break

    # Check if any existing processes have finished, and, if so, read
    # their output and then pop them off the list
    polls = [p.poll() for p in procs]
    for pr, po, lo in zip(procs, polls, logfiles):
        if po is not None:
            outs, errs = pr.communicate()
            fp = open(lo, "wb")
            if outs is not None: fp.write(outs)
            if errs is not None: fp.write(errs)
            fp.close()
    procs = [pr for pr, po in zip(procs, polls) if po is None]
    logfiles = [lo for lo, po in zip(logfiles, polls) if po is None]

# Wait for remaining processes to finish
while len(procs) > 0:
    polls = [p.poll() for p in procs]
    for pr, po, lo in zip(procs, polls, logfiles):
        if po is not None:
            outs, errs = pr.communicate()
            fp = open(lo, "wb")
            if outs is not None: fp.write(outs)
            if errs is not None: fp.write(errs)
            fp.close()
    procs = [pr for pr, po in zip(procs, polls) if po is None]
    logfiles = [lo for lo, po in zip(logfiles, polls) if po is None]
    
# Consolidate the output data into a smaller number of files for
# easier IO; to do this, we basically just stack all the different uh
# values into a single file
for m_ in m:
    for y_ in y:

        # Metadata files
        metaname = osp.join(outdir, "table_params_y{:1d}_m{:1d}.txt".
                            format(y_, m_))
        fp = open(metaname, 'w')
        fp.write(str(len(uh))+"\n")
        fp.write(str(nu)+"\n")
        fp.write(str(nq)+"\n")
        fp.write(str(ngex)+"\n")
        fp.write(str(ngex_lo)+"\n")
        fp.write(str(qmin)+"\n")
        fp.write(str(qmax)+"\n")
        fp.write(str(loggex_max)+"\n")
        for u_ in uh:
            fp.write("   {:5.2f}".format(u_))
        fp.write("\n")
        fp.close()

        # gex files
        dstname = osp.join(outdir,
                           "gextab_gex_y{:1d}_m{:1d}.bin".format(y_, m_))
        if overwrite or not osp.isfile(dstname):
            print("Consolidating to {:s}...".format(dstname))
            fpout = open(dstname, 'wb')
            for u_ in uh:
                inname = osp.join(tmpdir,
                                  "gextab_gex_uh{:f}_y{:1d}_m{:1d}.bin".
                                  format(u_, y_, m_))
                fpin = open(inname, 'rb')
                copyfileobj(fpin, fpout)
                fpin.close()
                print("   processed {:s}".format(inname))

        # gex_q files
        dstname = osp.join(outdir,
                           "gextab_q_y{:1d}_m{:1d}.bin".format(y_, m_))
        if overwrite or not osp.isfile(dstname):
            print("Consolidating to {:s}...".format(dstname))
            fpout = open(dstname, 'wb')
            for u_ in uh:
                inname = osp.join(tmpdir,
                                  "gextab_q_uh{:f}_y{:1d}_m{:1d}.bin".
                                  format(u_, y_, m_))
                fpin = open(inname, 'rb')
                copyfileobj(fpin, fpout)
                fpin.close()
                print("   processed {:s}".format(inname))

        # gex_u files
        dstname = osp.join(outdir,
                           "gextab_u_y{:1d}_m{:1d}.bin".format(y_, m_))
        if overwrite or not osp.isfile(dstname):
            print("Consolidating to {:s}...".format(dstname))
            fpout = open(dstname, 'wb')
            for u_ in uh:
                inname = osp.join(tmpdir,
                                  "gextab_u_uh{:f}_y{:1d}_m{:1d}.bin".
                                  format(u_, y_, m_))
                fpin = open(inname, 'rb')
                copyfileobj(fpin, fpout)
                fpin.close()
                print("   processed {:s}".format(inname))

        # q_u files
        dstname = osp.join(outdir,
                           "qtab_u_y{:1d}_m{:1d}.bin".format(y_, m_))
        if overwrite or not osp.isfile(dstname):
            print("Consolidating to {:s}...".format(dstname))
            fpout = open(dstname, 'wb')
            for u_ in uh:
                inname = osp.join(tmpdir,
                                  "qtab_u_uh{:f}_y{:1d}_m{:1d}.bin".
                                  format(u_, y_, m_))
                fpin = open(inname, 'rb')
                copyfileobj(fpin, fpout)
                fpin.close()
                print("   processed {:s}".format(inname))

        # q_q files
        dstname = osp.join(outdir,
                           "qtab_q_y{:1d}_m{:1d}.bin".format(y_, m_))
        if overwrite or not osp.isfile(dstname):
            print("Consolidating to {:s}...".format(dstname))
            fpout = open(dstname, 'wb')
            for u_ in uh:
                inname = osp.join(tmpdir,
                                  "qtab_q_uh{:f}_y{:1d}_m{:1d}.bin".
                                  format(u_, y_, m_))
                fpin = open(inname, 'rb')
                copyfileobj(fpin, fpout)
                fpin.close()
                print("   processed {:s}".format(inname))

        if (y_ <= m_):
            # q_gexlo_u files
            dstname = osp.join(outdir,
                               "qtab_gexlo_u_y{:1d}_m{:1d}.bin".format(y_, m_))
            if overwrite or not osp.isfile(dstname):
                print("Consolidating to {:s}...".format(dstname))
                fpout = open(dstname, 'wb')
                for u_ in uh:
                    inname = osp.join(tmpdir,
                                      "qtab_gexlo_u_uh{:f}_y{:1d}_m{:1d}.bin".
                                      format(u_, y_, m_))
                    fpin = open(inname, 'rb')
                    copyfileobj(fpin, fpout)
                    fpin.close()
                    print("   processed {:s}".format(inname))

            # q_gexlo_q files
            dstname = osp.join(outdir,
                               "qtab_gexlo_q_y{:1d}_m{:1d}.bin".format(y_, m_))
            if overwrite or not osp.isfile(dstname):
                print("Consolidating to {:s}...".format(dstname))
                fpout = open(dstname, 'wb')
                for u_ in uh:
                    inname = osp.join(tmpdir,
                                      "qtab_gexlo_q_uh{:f}_y{:1d}_m{:1d}.bin".
                                      format(u_, y_, m_))
                    fpin = open(inname, 'rb')
                    copyfileobj(fpin, fpout)
                    fpin.close()
                    print("   processed {:s}".format(inname))
                
