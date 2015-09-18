import numpy as np
from abc import ABCMeta, abstractmethod

class _AARS_Space(object):
    __metaclass__ = ABCMeta

    def __init__(self, aars_names, aa_names):
        """ This function takes a list or indexable ([] operator) class
            of determinable length of aars names and aa names. """
        self._aarss = aars_names
        self._aas = aa_names
        
    @abstractmethod
    def get_aars_to_aa_mapping(self, aars):
        """ This function should return a list of probabilities that
            the passed aars will charge each aa. """
        pass

    @abstractmethod
    def mutations_between(self, from_, to):
        """ This function should return an integer > 0 corresponding to 
            the numer of mutations it will take to go from the first aars
            to the second. """
        pass

    def mutational_neigbhors(self, from_):
        """ An optional function that should return a list of all
            the aarss that are one mutation away from the passed aars.
            This function is used to optimize getting the mutational probabilities
            and if not overridden will cause mutations_between to be called for
            aars which is bad if the number of aarss is large. """
        pass

    @abstractmethod
    def mutation_probability(self, from_, to, mu):
        """ This function should return the probability that the first
            aars will mutate into the second aars in one time step. """
        pass

    def get_aars_aa_map(self, aarss):
        """ This function takes a list of aarss and returns a
            a matrix of probabilities mapping aars to aa. """
        aars_aa_mapping = np.array([self.get_aars_to_aa_mapping(aars)
                                    for aars in aarss], dtype=np.float)

        row_sums = np.array([[row_sum] if row_sum > 0 else [1.]
                             for row_sum in aars_aa_mapping.sum(axis=1)])

        # normalize the matrix's rows to sum to 1 or 0
        return aars_aa_mapping / row_sums

    @property
    def names(self):
        return self._aarss

    @property
    def aa_names(self):
        return self._aas

    @property
    def count(self):
        return len(self._aarss)

    @property
    def aas(self):
        return len(self._aas)
