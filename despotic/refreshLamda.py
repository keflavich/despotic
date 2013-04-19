"""
This module contains the routine refreshLamda, which checks the
local installation of the Leiden Atomic and Molecular Database and
to if there are any files older than a specified age (default 6
months) and attempts to update old files from LAMDA on the web.
"""

########################################################################
# Copyright (C) 2013 Mark Krumholz
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
########################################################################

import os
import glob
import datetime as dt
from despoticError import despoticError
from fetchLamda import fetchLamda
from urlparse import urljoin

def refreshLamda(path=None, cutoffDate=None, cutoffAge=None, \
                     LamdaURL=None):
    """
    Refreshes LAMDA files by fetching new ones from the web

    Parameters
    ----------
    path : string
        path to the local LAMDA database; defaults to getting this
        information from the environment variable DESPOTIC_HOME
    cutoffDate : class datetime.date or class datetime.datetime
        a date or datetime specifying the age cutoff for updating
        files; files older than cutoffDate are updated, newer ones are
        not
    cutoffAge : class datetime.timedelta
        a duration between the present instant and the point in the
        past separating files that will be updated from files that
        will not be
    LamdaURL : string
        URL where LAMDA is located; defaults to the default value in
        fetchLamda

    Returns
    -------
    Nothing

    Remarks
    -------
    If the user sets both a cutoff age and a cutoff date, the date is
    used. If neither is set, the default cutoff age is 6 months.
    """

    # First check if we were given a cutoff date or age; if neither,
    # set a default age of 6 months
    if cutoffDate == None:
        if cutoffAge == None:
            # Set age to half a year ago
            cutoffDate = dt.datetime.today() - dt.timedelta(182.5)
        else:
            # Set cutoff age based on date
            cutoffDate = dt.datetime.today() - cutoffAge

    # Now locate the LAMDA path
    if path==None:
        if 'DESPOTIC_HOME' in os.environ:
            # Use environment variable if available
            path = os.path.join(os.environ['DESPOTIC_HOME'], \
                                    'LAMDA')
        else:
            # No path given, so just use relative location
            path = 'LAMDA'

    # Now go through all the .dat files in the directory
    for fname in glob.glob(os.path.join(path, '*.dat')):

        # Get time of last modification of this file
        mtime = dt.datetime.fromtimestamp(os.stat(fname).st_mtime)

        # Compare to cutoff date, and skip if the file is younger than
        # that; otherwise try to fetch file from the web
        if mtime < cutoffDate:
            fetchLamda(urljoin(LamdaURL, os.path.basename(fname)), \
                           path='', fileName=fname)

