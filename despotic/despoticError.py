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

########################################################################
# class despoticError
#
# A generic exception class for DESPOTIC-specific errors
########################################################################

class despoticError(Exception):
    """Exception raised by DESPOTIC-specific errors.

    Attributes
    ----------
    __init__ -- initialize
    """

    def __init__(self, message):
        Exception.__init__(self, message)
        self.message = message
