"""
This module provides helper classes to handle the computation of
reaction rates in chemical networks. They provide generic ways of
computing reaction rates from 2-body reactions, photoreactions, and
cosmic ray-induced rections, and allow users to create chemical
networks by registering sets of reactions and reactants.
"""

########################################################################
# Copyright (C) 2015 Mark Krumholz
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

import numpy as np
import scipy.sparse as sp
import inspect
from .abundanceDict import abundanceDict

########################################################################
# Useful little function, taken directly from
# http://stackoverflow.com/questions/4373631/sum-array-by-number-in-numpy. 
# It sums up arrays by a group index.
########################################################################
def sum_by_group(values, groups):
    order = np.argsort(groups)
    groups = groups[order]
    values = values[order]
    values.cumsum(out=values)
    index = np.ones(len(groups), 'bool')
    index[:-1] = groups[1:] != groups[:-1]
    values = values[index]
    groups = groups[index]
    values[1:] = values[1:] - values[:-1]
    return values, groups

########################################################################
# A general reaction class; it knows about stoichiometric factors, and
# how to compute reaction rates and rates of change of various species
# given a set of rate coefficients and abundances
########################################################################

class reaction_matrix(object):
    """
    This class provides a generic driver for computing rates of change
    of chemical species from chemical reactions. This class does the
    work of mapping from the reaction rates to rates of change in
    species abundances.

    Parameters
       specList: listlike of string
          List of all chemical species in the full reaction network,
          including those that are derived from conservation laws
          instead of being computed directly
       reactions : list of dict
          A list listing all the reactions to be registered; each
          entry in the list must be a dict containing the keys 'spec'
          and 'stoich', which list the species involved in the
          reaction and the stoichiometric factors for each species,
          respectively. Sign convention is that reactants on the left
          hand side have negative stoichiometric factors, those on the
          right hand side have positive factors.
       sparse : bool
          If True, the reaction rate matrix is represented as a
          sparse matrix; otherwise it is a dense matrix. This has no
          effect on results, but depending on the chemical network it
          may lead to improved execution speed and/or reduced memory
          usage.

    Examples
       To describe the reaction C + O -> CO, the correct dict entry is::

          { 'spec' : ['C', 'O', 'CO'], 'stoich' : [-1, -1, 1] }

       To describe the reaction H + H -> H2, the dict is::

          { 'spec' : ['H', 'H2'], 'stoich' : [-2, 1] }
    """

    def __init__(self, specList, reactions, sparse=False):
        """
        This creates a reactions object.

        Parameters
           specList: listlike of string
              List of all chemical species in the full reaction network,
              including those that are derived from conservation laws
              instead of being computed directly
           reactions : list of dict
              A list listing all the reactions to be registered; each
              entry in the list must be a dict containing the keys 'spec'
              and 'stoich', which list the species involved in the
              reaction and the stoichiometric factors for each species,
              respectively. Sign convention is that reactants on the left
              hand side have negative stoichiometric factors, those on the
              right hand side have positive factors.
           sparse : bool
              If True, the reaction rate matrix is represented as a
              sparse matrix; otherwise it is a dense matrix. This has no
              effect on results, but depending on the chemical network it
              may lead to improved execution speed and/or reduced memory
              usage.

        Examples
           To describe the reaction C + O -> CO, the correct dict entry is::

              { 'spec' : ['C', 'O', 'CO'], 'stoich' : [-1, -1, 1] }

           To describe the reaction H + H -> H2, the dict is::

              { 'spec' : ['H', 'H2'], 'stoich' : [-2, 1] }
        """

        # Count number of species and number of reactions
        nspec = len(specList)
        nreact = len(reactions)

        # Create a dummy abundanceDict we can use to get the index
        # mapping from strings to species numbers
        self.specDict = abundanceDict(specList, np.zeros(nspec))

        # Create the reaction rate matrix; this will have the property
        # that dx/dt = M * R, where dx/dt is a vector containing the
        # rate of change of the abundance of each species, M is the
        # reaction rate matrix, and R is the rate of every individual
        # reaction; at the same time, store the index and
        # stoichiometric factor for every reactant on the LHS of the
        # reaction; NB: we need to be careful to handle properly the
        # case where the same species appears on both sides of the
        # reaction, to ensure that we properly cancel the
        # stoichiometric factors
        lhs = []
        lhsstoich = []
        nlhsmax = 0
        if sparse:

            # Case where we treat the matrix as sparse

            # Pack the data into vectors to give matrix craetion
            # routine
            mat_i = np.zeros(0, dtype=int)
            mat_j = np.zeros(0, dtype=int)
            stoich = np.zeros(0, dtype=int)
            for i, r in enumerate(reactions):
                spec = self.specDict.index(r['spec'])
                stoich_tmp = np.array(r['stoich'], dtype=int)
                # These next two lines handle the case where the same
                # species appears on the LHS and RHS of the reactions
                spec_unq, inv = np.unique(spec, return_inverse=True)
                stoich_unq = sum_by_group(stoich_tmp, inv)[0]
                mat_i = np.append(mat_i, spec_unq)
                mat_j = np.append(mat_j, np.zeros(len(spec_unq), dtype=int)+i)
                stoich = np.append(stoich, stoich_unq)
                lhs.append(spec[stoich_tmp < 0])
                lhsstoich.append(-stoich_tmp[stoich_tmp < 0])
                if nlhsmax < np.sum(stoich_tmp < 0):
                    nlhsmax = np.sum(stoich_tmp < 0)

            # Build matrix
            self.mat = sp.coo_matrix((stoich, (mat_i, mat_j)), 
                                     shape=(nspec, nreact), dtype='int')
            self.mat = sp.csr_matrix(self.mat)

        else:

            # Dense matrix case
            self.mat = np.matrix(np.zeros((nspec, nreact), dtype='int'))

            # Load the reaction matrix with the reactions
            for i, r in enumerate(reactions):
                spec = self.specDict.index(r['spec'])
                np.add.at(self.mat[:, i], 
                          (spec, np.zeros(len(spec), dtype=int)), 
                          r['stoich'])
                spec_stoich = np.array(r['stoich'])
                lhs.append(spec[spec_stoich < 0])
                lhsstoich.append(-spec_stoich[spec_stoich < 0])
                if nlhsmax < np.sum(spec_stoich < 0):
                    nlhsmax = np.sum(spec_stoich < 0)

        # For each reaction, record the indices and stoichiometric
        # factors on the left hand side; record these in a square
        # matrix; also record the sum of the LHS stoichiometric factors
        self.lhsidx = np.zeros((nreact, nlhsmax), dtype=int)
        self.lhsstoich = np.zeros((nreact, nlhsmax), dtype=int)
        for i in np.arange(nreact):
            self.lhsidx[i,:len(lhs[i])] = lhs[i]
            self.lhsstoich[i,:len(lhs[i])] = lhsstoich[i]
        self.lhssum = np.sum(self.lhsstoich, axis=1)

    def dxdt(self, x, n, ratecoef):
        """
        This returns the rate of change of species abundances given a
        set of rate coefficients.

        Parameters
           x : array(N_species)
              array of current species abundances
           n : float
              number density of H nuclei
           ratecoef : array(N_reactions)
              rate coefficients for each reaction; reaction rate per unit
              volume = ratecoef * product of densities of reactants;
              dxdt = reaction rate / unit volume / n

        Returns
           dxdt : array(N)
              rate of change of all species abundances
        """

        # Get rates from rate coefficients
        rates = ratecoef * \
                np.prod(x[self.lhsidx]**self.lhsstoich, axis=1) * \
                n**(self.lhssum-1)

        # Compute dxdt
        dxdt = self.mat.dot(rates)

        # Return
        return dxdt
        


########################################################################
# Cosmic ray reactions
########################################################################
class cr_reactions(reaction_matrix):
    """
    The cr_reactions class is used to handle computation of
    cosmic ray-induced reaction rates. In addition to the constructor,
    the class has only a single method: dxdt, which returns the
    reaction rates.

    Parameters
       specList: listlike of string
          List of all chemical species in the full reaction network,
          including those that are derived from conservation laws
          instead of being computed directly
       reactions : list of dict
          A list listing all the reactions to be registered; each
          entry in the list must be a dict containing the keys 'spec'
          'stoich', and 'rate', which list the species involved in the
          reaction, the stoichiometric factors for each species, and
          the reaction rate per primary CR ionization,
          respectively. Sign convention is that reactants on the left 
          hand side have negative stoichiometric factors, those on the
          right hand side have positive factors.
       sparse : bool
          If True, the reaction rate matrix is represented as a
          sparse matrix; otherwise it is a dense matrix. This has no
          effect on results, but depending on the chemical network it
          may lead to improved execution speed and/or reduced memory
          usage.

    Examples
       To list the reaction cr + H -> H+ + e-, the dict entry should be::

          { 'spec' : ['H', 'H+', 'e-'], 'stoich' : [-1, 1, 1],
            'rate' : 1.0 }
    """

    def __init__(self, specList, reactions, sparse=False):
        """
        This creates a cosmic ray reaction rate computation class, for
        computing rates for reactions of the form CR + ... -> ....

        Parameters
           specList: listlike of string
              List of all chemical species in the full reaction network,
              including those that are derived from conservation laws
              instead of being computed directly
           reactions : list of dict
              A list listing all the reactions to be registered; each
              entry in the list must be a dict containing the keys 'spec'
              'stoich', and 'rate', which list the species involved in the
              reaction, the stoichiometric factors for each species, and
              the reaction rate per primary CR ionization,
              respectively. Sign convention is that reactants on the left 
              hand side have negative stoichiometric factors, those on the
              right hand side have positive factors.
           sparse : bool
              If True, the reaction rate matrix is represented as a
              sparse matrix; otherwise it is a dense matrix. This has no
              effect on results, but depending on the chemical network it
              may lead to improved execution speed and/or reduced memory
              usage.

        Examples
           To list the reaction cr + H -> H+ + e-, the dict entry
           should be::

              { 'spec' : ['H', 'H+', 'e-'], 'stoich' : [-1, 1, 1],
                'rate' : 1.0 }
        """

        # Call the parent constructor
        super(cr_reactions, self).__init__(specList, reactions, 
                                           sparse=sparse)

        # Record the rates
        self.rate_fac = np.array([r['rate'] for r in reactions])


    def dxdt(self, x, n, ionrate):
        """
        This function returns the time derivative of the abundances x
        for a given cosmic ray ionization rate.

        Parameters
           x : array(N)
              array of current species abundances
           n : float
              number density of H nuclei; only used if some reactions
              have multiple species on the LHS, otherwise can be set to
              any positive value
           ionrate : float
              cosmic ray primary ionization rate

        Returns
           dxdt: array(N)
              rate of change of all species abundances
        """

        # Get rate coefficients and rates
        ratecoef = self.rate_fac * ionrate

        # Compute dxdt
        dxdt = super(cr_reactions, self).dxdt(x, n, ratecoef)

        # Return
        return dxdt


########################################################################
# Photoreactions
########################################################################
class photoreactions(reaction_matrix):
    """
    The photoreactions class is used to handle computation of
    photon-induced reaction rates. Generally, it returns reaction
    rates for any reaction of the form
    gamma + ... -> ...
    with a rate that scales as the ISRF strength, parameterized in
    units of the Habing (1968) field, multiplied by dust and gas
    shielding factors. In addition to the constructor, the class has
    only a single method: dxdt, which returns the reaction rates.

    Parameters
       specList: listlike of string
          List of all chemical species in the full reaction network,
          including those that are derived from conservation laws
          instead of being computed directly
       reactions : list of dict
          A list listing all the reactions to be registered; each
          entry in the list must be a dict containing the keys:

          * 'spec' : list 
             list of strings giving the species involved in the reaction
          * 'stoich' : list
             list of int stoichiometric factor for each species
          * 'rate' : float
             reaction rate per target in a chi = 1 radiation field
          * 'av_fac' : float
             optical depth per unit A_V
          * 'shield_fac' : (optional) callable 
             callable to compute the
             shielding factor; see the dxdt method for an
             explanation of how to specify its arguments

          Reaction rates per target are given by 
             chi * rate * shield_fac * exp(-av_fac * A_V)
       sparse : bool
          If True, the reaction rate matrix is represented as a
          sparse matrix; otherwise it is a dense matrix. This has no
          effect on results, but depending on the chemical network it
          may lead to improved execution speed and/or reduced memory
          usage.
    """

    def __init__(self, specList, reactions, sparse=False):
        """
        This creates a photoreaction rate computation class, for
        computing rates for reactions of the form gamma + ... ->
        ...

        Parameters
           specList: listlike of string
              List of all chemical species in the full reaction network,
              including those that are derived from conservation laws
              instead of being computed directly
           reactions : list of dict
              A list listing all the reactions to be registered; each
              entry in the list must be a dict containing the keys:

              * 'spec' : list 
                 list of strings giving the species involved in the reaction
              * 'stoich' : list
                 list of int stoichiometric factor for each species
              * 'rate' : float
                 reaction rate per target in a chi = 1 radiation field
              * 'av_fac' : float
                 optical depth per unit A_V
              * 'shield_fac' : (optional) callable 
                 callable to compute the
                 shielding factor; see the dxdt method for an
                 explanation of how to specify its arguments

              Reaction rates per target are given by 
                 chi * rate * shield_fac * exp(-av_fac * A_V)
           sparse : bool
              If True, the reaction rate matrix is represented as a
              sparse matrix; otherwise it is a dense matrix. This has no
              effect on results, but depending on the chemical network it
              may lead to improved execution speed and/or reduced memory
              usage.
        """

        # Call the parent constructor
        super(photoreactions, self).__init__(specList, reactions,
                                             sparse=sparse)

        # Record the rates
        self.rate_fac = np.array([r['rate'] for r in reactions])

        # Save the A_V factors
        self.av_fac = np.array([r['av_fac'] for r in reactions])

        # Save the shielding factor functions
        self.shield_func = []
        self.shield_func_idx = []
        for i, r in enumerate(reactions):
            if 'shield_fac' in r.keys():
                self.shield_func.append(r['shield_fac'])
                self.shield_func_idx.append(i)

        # For each shielding function, record its call signature; we
        # need to know if it expects positional arguments and keyword
        # arguments so we can give it those when it is called
        self.shield_func_args = []
        self.shield_func_kw_args = []
        for f in self.shield_func:
            argspec = inspect.getargspec(f)
            self.shield_func_args.append(
                (argspec.args is not None) or 
                (argspec.varargs is not None))
            self.shield_func_kw_args.append(
                argspec.keywords is not None)


    def dxdt(self, x, n, chi, AV, shield_args=None,
             shield_kw_args=None):
        """
        This function returns the time derivative of the abundances x
        for a given ISRF, extinction, and gas shielding factor.

        Parameters
           x : array(N)
              array of current species abundances
           n : float
              number density of H nuclei; only used if some reactions
              have multiple species on the LHS, otherwise can be set to
              any positive value
           chi : float
              ISRF strength normalized to the solar neighborhood value
           AV : float
              visual extinction to apply to the ISRF
           shield_args : list
              list of argument lists to be passed to the shielding
              functions; *(shield_args[0]) will be passed to the first
              shielding function, *(shield_args[1]) to the second
              shielding function, etc. Arguments must be supplied for
              shielding functions that expect them, or an error is
              raised.
           shield_kw_args : list
              same as shield_args, but instead of the list items being
              lists themselves, they are dicts, and are passed as
              keywords to the shielding functions as **(shield_kw_args[0]),
              **(shield_kw_args[1]), etc. If this is not None, the list
              must have the same number of elements as the number of
              shielding functions.

        Returns
           dxdt: array(N)
              rate of change of all species abundances
        """

        # Evaluate the shielding function with the appropriate call
        # signature
        if len(self.shield_func) > 0:
            fshield = np.zeros(len(self.shield_func))
            for i, f, has_args, has_kw in \
                zip(np.arange(len(self.shield_func)), self.shield_func, 
                    self.shield_func_args,
                    self.shield_func_kw_args):
                has_kw = has_kw and (shield_kw_args is not None)
                if has_kw:
                    has_kw = has_kw and (shield_kw_args[i] is not None)
                if has_args and has_kw:
                    fshield[i] = f(*(shield_args[i]), **(shield_kw_args[i]))
                elif has_args:
                    fshield[i] = f(*(shield_args[i]))
                elif has_kw:
                    fshield[i] = f(**(shield_kw_args[i]))
                else:
                    fshield[i] = f()
 
        # Get rate coefficients
        ratecoef = self.rate_fac * chi * np.exp(-AV*self.av_fac)
        if len(self.shield_func) > 0:
            ratecoef[self.shield_func_idx] *= fshield

        # Compute dxdt
        dxdt = super(photoreactions, self).dxdt(x, n, ratecoef)

        # Return
        return dxdt


########################################################################
# Reactions with grains
########################################################################
class gr_reactions(reaction_matrix):
    """
    The gr_reactions class is used to handle computation of reaction
    rates where some of the reactions are grain-catalyzed, and thus
    have rates that scale as the dust abundance.
    """

    def __init__(self, specList, reactions, sparse=False):
        """
        This creates a photoreaction rate computation class, for
        computing rates for reactions of the form gamma + ... ->
        ...

        Parameters
        ----------
        specList: listlike of string
           List of all chemical species in the full reaction network,
           including those that are derived from conservation laws
           instead of being computed directly
        reactions : list of dict
           A list listing all the reactions to be registered; each
           entry in the list must be a dict containing the keys:
              'spec' : list of the species involved in the reaction
              'stoich' : stoichiometric factor for each species
              'grain' : (optional) if this key exists and is True, the
                 reaction is grain catalyzed and the rate scales as
                 the dust abundance; if this key is absence, the
                 reaction is assumed not to be grain-mediated
        sparse : bool
           If True, the reaction rate matrix is represented as a
           sparse matrix; otherwise it is a dense matrix. This has no
           effect on results, but depending on the chemical network it
           may lead to improved execution speed and/or reduced memory
           usage.

        Returns
        -------
        Nothing
        """

        # Call the parent constructor
        super(gr_reactions, self).__init__(specList, reactions,
                                             sparse=sparse)

        # Record which reactions, if any, are grain catalyzed
        grain = []
        for r in reactions:
            if 'grain' in r.keys():
                if r['grain']:
                    grain.append(True)
                else:
                    grain.append(False)
            else:
                grain.append(False)
        self.grain = np.where(grain)[0]

    def dxdt(self, x, n, ratecoef, Zd):
        """
        This returns the rate of change of species abundances given a
        set of rate coefficients.

        Parameters
        ----------
        x : array(Nspec)
           array of current species abundances
        n : float
           number density of H nuclei
        ratecoef : array(Nreactions)
           rate coefficients for each reaction, at Solar neighborhood
           grain abundances
        Zd : float
           dust abundances scaled to Solar neighborhood;
           grain-mediated reaction rate coefficients will be scaled by
           this value

        Returns
        -------
        dxdt: array(N)
           rate of change of all species abundances
        """

        # Scale grain-mediated reaction rates
        rcoef = np.copy(ratecoef)
        rcoef[self.grain] *= Zd

        # Compute dxdt
        dxdt = super(gr_reactions, self).dxdt(x, n, rcoef)

        # Return
        return dxdt
