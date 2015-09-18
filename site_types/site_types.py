from abc import ABCMeta, abstractmethod
import numpy as np
from collections import OrderedDict
from copy import copy

class _Site_Types(object):
    __metaclass__ = ABCMeta

    def __init__(self, phi, amino_acids, sites, weights):
        if not 0 < phi < 1:
            raise ValueError("The domain of phi is (0,1).")

        if len(weights) != len(sites):
            raise ValueError("There must be a weight for every site.")

        if any(map(lambda weight: weight < 1, weights)):
            raise ValueError("Every site weight must be greater than 1.")

        if len(set(map(lambda site: site[0], sites))) != len(sites):
            raise ValueError("All sites must have unique names.")

        if len(set(map(lambda aa: aa[0], amino_acids))) != len(amino_acids):
            raise ValueError("All amino acids must have unique names.")

        self._phi = phi
        self._weights = weights
        self._sites = OrderedDict(sites)
        self._aas = OrderedDict(amino_acids)
        
        self._build_matrices()

    @property
    def fitness_matrix(self):
        return self._fitness_matrix

    @property
    def distance_matrix(self):
        return self._distance_matrix

    @property
    def weights(self):
        return self._weights

    @property
    def phi(self):
        return self._phi

    @property
    def sites(self):
        return copy(self._sites)

    @property
    def aas(self):
        return copy(self._aas)

    def _build_matrices(self):
        self._distance_matrix = np.array([[self.distance(site, aa) for aa in self._aas.values()]
                                          for site in self._sites.values()], dtype=np.float)

        self._fitness_matrix = self.phi ** self._distance_matrix

        # Make the matrices read only so that they can be passed
        # without needing to make a copy
        self._distance_matrix.setflags(write=False)
        self._fitness_matrix.setflags(write=False)

    @abstractmethod
    def distance(self, site, amino_acid):
        pass
