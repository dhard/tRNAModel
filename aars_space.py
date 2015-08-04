import numpy as np

class AARS_Space(object):
    def __init__(self, aars_names, aars_aa_mapping, aars_mutation_matrix):
        self._aarss = len(aars_names)

        if self._aarss != aars_aa_mapping.shape[0]:
            raise ValueError("There must be as many AARS names as there are AARSs in the "
                             "AARS-AARS mapping matrix.")

        if self._aarss != aars_mutation_matrix.shape[0] or self._aarss != aars_mutation_matrix.shape[1]:
            raise ValueError("There must be as many AARS names as there are AARSs in the "
                             "mutation matrix.")
            
        if not np.allclose(1.0, aars_mutation_matrix.sum(axis=0)):
            raise ValueError("All rows must sum to 1 in the mutation matrix.")

        self._aars_aa_mapping = aars_aa_mapping
        self._mutation_matrix = aars_mutation_matrix
        self._names = tuple(aars_names) # Force ordering and immutability

        self._aars_aa_mapping.setflags(write=False)
        self._mutation_matrix.setflags(write=False)


    @property
    def aars_aa_mapping(self):
        return self._aars_aa_mapping

    @property
    def mutation_matrix(self):
        return self._mutation_matrix

    @property
    def names(self):
        return self._names

    @property
    def count(self):
        return self._aarss
