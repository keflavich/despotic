"""
This module provides the routine fetchLambda, a utility for
downloading LAMDA files from the web.
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
from urllib2 import *
from urlparse import urljoin
from despoticError import despoticError

# Default location of LAMDA database
lamdaURL = 'http://home.strw.leidenuniv.nl/~moldata/'

def fetchLamda(inputURL, path=None, fileName=None):
    """
    Routine to download LAMDA files from the web.

    Parameters
    ----------
    inputURL : string
        URL of LAMDA file containing data on this species; if this
        does not begin with "http://", indicating it is a URL, then
        this is assumed to be a filename within LAMDA, and a default
        URL is appended
    path : string
        relative or absolute path at which to store the file; if not
        set, the current directory is used; if the specified path does
        not exist, it is created
    fileName : string
        name to give to file; if not set, defaults to the same as the
        name in LAMDA

    Returns
    -------
    fname : string
        local file name to which downloaded file was written; if URL
        cannot be opened, None is returned instead

    Raises
    ------
    despoticError if output file cannot be written
    """

    # Check if the URL we've been passed starts with http://. If not,
    # prepend a URL, taken from the environment or just using
    # the default.
    if inputURL[:7] == 'http://':
        emitterURL = inputURL
    else:
        if 'DESPOTIC_LAMDAURL' in os.environ:
            # Set URL from environment variable
            baseURL = urljoin(os.environ['DESPOTIC_LAMDAURL'], \
                                  'datafiles/')
        else:
            # Use default
            baseURL = urljoin(lamdaURL, 'datafiles/')
        emitterURL = urljoin(baseURL, inputURL)

    # Try to open URL
    try:
        urlPtr = urlopen(emitterURL)
    except HTTPError:
        return None

    # Print message that we're downloading
    print "Fetching LAMDA datafile from "+emitterURL+"..."

    # First make sure there's a valid path to where we want to
    # write the file, and create one if necessary
    if path == None:
        if 'DESPOTIC_HOME' in os.environ:
            # Use environment variable if available
            path = \
                os.path.join(os.environ['DESPOTIC_HOME'], \
                                 'LAMDA')
        else:
            # No path given, so just use relative location
            path = 'LAMDA'
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except:
            despoticError, "could not create path "+path

    # Construct file name to use to store file, unless one was given
    # as an option
    if fileName == None:
        fileName = emitterURL.rsplit('/',1)[1]

    # Now open an emitter file, and write from the URL to it
    try:
        fpWrite = open(os.path.join(path, fileName), 'w')
        fpWrite.write(urlPtr.read())
    except:
        raise despoticError, "could not write to file " + \
            os.path.join(path, fileName)

    # Close file and URL
    fpWrite.close()
    urlPtr.close()

    # Return success
    return os.path.join(path, fileName)
