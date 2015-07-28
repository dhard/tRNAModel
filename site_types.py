from abc import ABCMeta, abstractmethod
import numpy as np

class Site_Types(object):
    __metaclass__ = ABCMeta

    def __init__(self, phi, sites, site_weights):
        assert 0 < phi < 1, "Intensity of selection (phi) must be between 0 and 1."
        assert all([weight >= 1 for weight in site_weights]), "Site weights must 1 or larger"

        assert len(sites) == len(site_weights), "There must be the same number of site weights as sites."

        self._sites = tuple(sites)
        self._site_weights = tuple(site_weights)
        self._phi = phi

        self._distance_matrix = None
        self._fitness_matrix = None

    def _calculate_distance_matrix(self):
        self._distance_matrix = np.array([[self.calculate_distance(alpha, beta) for beta in self.sites]
                                          for alpha in self.sites]) # Ring space from cmcpy

        self._fitness_matrix = self.phi ** self._distance_matrix # Eq1 SellaArdell06

    @abstractmethod
    def calculate_distance(self, alpha, beta):
        pass

    @property
    def distance_matrix(self):
        if self._distance_matrix is None:
            self._calculate_distance_matrix()

        return self._distance_matrix.copy()

    @property
    def fitness_matrix(self):
        if self._fitness_matrix is None:
            self._calculate_distance_matrix()

        return self._fitness_matrix.copy()

    @property
    def site_weights(self):
        return self._site_weights

    @property
    def sites(self):
        return self._sites

    @property
    def phi(self):
        return self._phi

class Ring_Site_Types(Site_Types):
    def __init__(self, phi, sites, site_weights):
        super(Ring_Site_Types, self).__init__(phi, sites, site_weights)

    def calculate_distance(self, alpha, beta):
        return min(abs(alpha - beta), (1 - abs(alpha - beta))) # Ring space from cmcpy

