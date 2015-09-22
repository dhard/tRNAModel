import numpy as np

from aars_space import AARS_Space

class Id_AARS_Space(AARS_Space):
    """ This class maps aars i to aa i. 
        AARS are unmutatable. """

    def __init__(self, aars_names, aa_names):
        super(Id_AARS_Space, self).__init__(aars_names, aa_names)

    def get_aars_to_aa_mapping(self, aars):
        mapping = np.zeros(self.aas, dtype=np.float)
        mapping[aars] = 1.0

        return mapping

    def mutations_between(self, from_, to):
        return 0

    def mutational_neighbors(self, from_):
        return []

    def mutation_probability(self, from_, to, mu):
        if from_ == to:
            return 1.0
        return 0.0
