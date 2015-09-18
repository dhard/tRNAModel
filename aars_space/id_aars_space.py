import numpy as np

class Id_AARS_Space(_AARS_Space):
    """ This class maps aars i to aa i. 
        AARS are unmutatable. """

    def __init__(self, aars_names, aa_names):
        super(self, _AARS_Space).__init__(aars_names, aa_names)

    @abstractmethod
    def get_aars_to_aa_mapping(self, aars):
        mapping = np.zeros(self.aas, dtype=np.float)
        mapping[aars] = 1.0

        return mapping

    @abstractmethod
    def mutations_between(self, from_, to):
        return 0

    def mutational_neigbhors(self, from_):
        return []

    @abstractmethod
    def mutation_probability(self, from_, to, mu):
        return 0.0
