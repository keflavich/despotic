#!python
#cython: boundscheck=False

"""
This module provides some improved speed for some of the functions
that hang off the collPartner class.
"""

import numpy as np
cimport numpy as np
NPFLOAT = np.float_
NPINT = np.int
ctypedef np.float_t NPFLOAT_T
ctypedef np.int_t NPINT_T


########################################################################
# These functions do interpolation on the collision rates; four
# versions exist:
# all_scalar = 1 temperature, all collisions
# all_vector = 1d array of temperatures, all collisions
# some_scalar = 1 temperature, specified list of states
# some_vector = 1d array of temperatures, specified list of states
########################################################################

#@cython.boundscheck(False)
cpdef np.ndarray[NPFLOAT_T, ndim=1] colRates_all_scalar(
    np.ndarray[NPFLOAT_T, ndim=2] logColRate, 
    np.ndarray[NPFLOAT_T, ndim=1] logTempTable,
    double temp):

    # Declare variables
    cdef int ntrans = logColRate.shape[0]
    cdef int ntemp = logColRate.shape[1]
    cdef np.ndarray[NPFLOAT_T, ndim=1] logrates = np.zeros(ntrans)
    cdef int i, j
    cdef double logtemp = np.log(temp)

    # Check for temperatures above or below the ends of the table
    if logtemp <= logTempTable[0]:
        logrates = logColRate[:,0]
        return np.exp(logrates)
    elif logtemp >= logTempTable[ntemp-1]:
        logrates = logColRate[:,ntemp-1]
        return np.exp(logrates)

    # Loop over transitions
    for i in range(ntrans):
        for j in range(1,ntemp):
            if logtemp < logTempTable[j]:
                break
        logrates[i] = logColRate[i,j-1] + \
                      (logColRate[i,j]-logColRate[i,j-1]) / \
                      (logTempTable[j]-logTempTable[j-1]) * \
                      (logtemp - logTempTable[j-1])

    # Return
    return np.exp(logrates)


cpdef np.ndarray[NPFLOAT_T, ndim=2] colRates_all_vector(
    np.ndarray[NPFLOAT_T, ndim=2] logColRate, 
    np.ndarray[NPFLOAT_T, ndim=1] logTempTable,
    np.ndarray[NPFLOAT_T, ndim=1] temp):

    # Declare variables
    cdef int ntrans = logColRate.shape[0]
    cdef int ntemp = logColRate.shape[1]
    cdef np.ndarray[NPFLOAT_T, ndim=2] logrates \
        = np.zeros((ntrans, temp.size))
    cdef int i, j, k
    cdef np.ndarray[NPFLOAT_T, ndim=1] logtemp = np.log(temp)

    # Loop over input temperatures
    for k in range(temp.size):

        # Check for temperatures above or below the ends of the table
        if logtemp[k] <= logTempTable[0]:
            logrates[:,k] = logColRate[:,0]
        elif logtemp[k] >= logTempTable[ntemp-1]:
            logrates[:,k] = logColRate[:,ntemp-1]
        else:
            # Loop over transitions
            for i in range(ntrans):
                for j in range(1,ntemp):
                    if logtemp[k] < logTempTable[j]:
                        break
                logrates[i,k] \
                    = logColRate[i,j-1] + \
                    (logColRate[i,j]-logColRate[i,j-1]) / \
                    (logTempTable[j]-logTempTable[j-1]) * \
                    (logtemp[k] - logTempTable[j-1])

    # Return
    return np.exp(logrates)


cpdef np.ndarray[NPFLOAT_T, ndim=1] colRates_some_scalar(
    np.ndarray[NPINT_T, ndim=1] colUpper,
    np.ndarray[NPINT_T, ndim=1] colLower,
    np.ndarray[NPFLOAT_T, ndim=2] logColRate, 
    np.ndarray[NPFLOAT_T, ndim=1] logTempTable,
    double temp,
    np.ndarray[NPINT_T, ndim=2] transition):

    # Declare variables
    cdef int i, j
    cdef int ntrans = transition.shape[1]
    cdef int ntemp = logTempTable.size
    cdef np.ndarray[NPFLOAT_T, ndim=1] rates = np.zeros(ntrans)
    cdef double logtemp = np.log(temp)
    cdef np.ndarray[NPINT_T, ndim=1] idxtmp
    cdef int idx

    # Loop over transitions
    for i in range(ntrans):

        # Find the entry in the collision rate table for this pair
        idxtmp = np.where(colUpper == transition[0,i])[0]
        if idxtmp.size == 0:
            continue

        # Check for temperatures above or below the ends of the table
        idx = idxtmp[0]
        if logtemp <= logTempTable[0]:
            rates[i] = np.exp(logColRate[idx,0])
        elif logtemp >= logTempTable[ntemp-1]:
            rates[i] = np.exp(logColRate[idx,ntemp-1])
        else:
            # Interpolate in temperature
            for j in range(1,ntemp):
                if logtemp < logTempTable[j]:
                    break
            rates[i] \
                = np.exp(logColRate[idx,j-1] + 
                         (logColRate[idx,j]-logColRate[idx,j-1]) / 
                         (logTempTable[j]-logTempTable[j-1]) * 
                         (logtemp - logTempTable[j-1]))

    # Return
    return rates


cpdef np.ndarray[NPFLOAT_T, ndim=1] colRates_some_vector(
    np.ndarray[NPINT_T, ndim=1] colUpper,
    np.ndarray[NPINT_T, ndim=1] colLower,
    np.ndarray[NPFLOAT_T, ndim=2] logColRate, 
    np.ndarray[NPFLOAT_T, ndim=1] logTempTable,
    np.ndarray[NPFLOAT_T, ndim=1] temp,
    np.ndarray[NPINT_T, ndim=2] transition):

    # Declare variables
    cdef int i, j, k
    cdef int ntrans = transition.shape[1]
    cdef int ntemp = logTempTable.size
    cdef np.ndarray[NPFLOAT_T, ndim=2] rates = np.zeros((ntrans, temp.size))
    cdef np.ndarray[NPFLOAT_T, ndim=1] logtemp = np.log(temp)
    cdef np.ndarray[NPINT_T, ndim=1] idxtmp
    cdef int idx

    # Loop over transitions
    for i in range(ntrans):

        # Find the entry in the collision rate table for this pair
        idxtmp = np.where(colUpper == transition[0,i])[0]
        if idxtmp.size == 0:
            continue

        # Loop over temperatures
        for k in range(temp.size):

            # Check for temperatures above or below the ends of the table
            idx = idxtmp[0]
            if logtemp[k] <= logTempTable[0]:
                rates[i,k] = np.exp(logColRate[idx,0])
            elif logtemp[k] >= logTempTable[ntemp-1]:
                rates[i,k] = np.exp(logColRate[idx,ntemp-1])
            else:
                # Interpolate in temperature
                for j in range(1,ntemp):
                    if logtemp[k] < logTempTable[j]:
                        break
                rates[i,k] \
                    = np.exp(logColRate[idx,j-1] + 
                             (logColRate[idx,j]-logColRate[idx,j-1]) / 
                             (logTempTable[j]-logTempTable[j-1]) * 
                             (logtemp[k] - logTempTable[j-1]))

    # Return
    return rates
