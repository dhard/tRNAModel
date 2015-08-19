import numpy as np

class AARS_Space(object):
    def __init__(self, aars_names, aa_names, aars_aa_map, mutations_between, aars_mutation_matrix):
        self._aarss = len(aars_names)
        self._aas = len(aa_names)

        #if self._aarss != aars_aa_mapping.shape[0]:
        #    raise ValueError("There must be as many AARS names as there are AARSs in the "
        #                     "AARS-AARS mapping matrix.")

        if self._aarss != aars_mutation_matrix.shape[0] or self._aarss != aars_mutation_matrix.shape[1]:
            raise ValueError("There must be as many AARS names as there are AARSs in the "
                             "mutation matrix.")
            
        #if not np.allclose(1.0, aars_mutation_matrix.sum(axis=1)):
        #    raise ValueError("All rows must sum to 1 in the mutation matrix.")

        self._aars_aa_map = aars_aa_map
        self._mutation_matrix = aars_mutation_matrix
        self._mutations_between = mutations_between
        self._names = tuple(map(str, aars_names)) # Force ordering and immutability
        self._aa_names = tuple(map(str, aa_names))

        self._mutation_matrix.setflags(write=False)

    def get_aa_from_aars(self, aars):
        return self._aars_aa_map(aars)

    def get_aars_aa_map(self, aarss):
        aars_aa_mapping = np.zeros((len(aarss), self.aas), dtype=np.float)

        for i, aars in enumerate(aarss):
            aa = self.get_aa_from_aars(aars)
            aars_aa_mapping[i, aa] = 1.0

        row_sums = np.array([[row_sum] if row_sum > 0 else [1.]
                             for row_sum in aars_aa_mapping.sum(axis=1)])
                            
        return aars_aa_mapping / row_sums


    def get_aars_aa_probability(self, aars, aa):
        prob = self._aars_aa_map(aars, aa)

        assert 0 <= prob <= 1, "AARS to amino acid mapping probability must be in the range [0,1]."

        return prob

    def mutations_between(self, aars1, aars2):
        mutations = self._mutations_between(aars1, aars2)

        assert mutations >= 0, "Mutations between two AARSs must be greater than 0."
        
        return mutations

    @property
    def mutation_matrix(self):
        return self._mutation_matrix

    @property
    def names(self):
        return self._names

    @property
    def aa_names(self):
        return self._aa_names

    @property
    def count(self):
        return self._aarss

    @property
    def aas(self):
        return self._aas
