import numpy as np
from itertools import permutations
from code import Code

class Code_Mutator(object):
    def __init__(self, initial_code, aars_mutation_matrix, trna_mutation_matrix):
        if not np.allclose(1.0, aars_mutation_matrix.sum(axis=1)):
            raise ValueError("All rows in the aars mutation matrix must sum to 1.")

        if not np.allclose(1.0, trna_mutation_matrix.sum(axis=1)):
            raise ValueError("All rows in the trna mutation matrix must sum to 1.")

        self._initial_code = initial_code
        self._aars_mutation_matrix = aars_mutation_matrix
        self._trna_mutation_matrix = trna_mutation_matrix

        self._mutation_probabilities = None

    def get_mutation_probabilities(self):
        if self._mutation_probabilities is None:
            aarss_perm = permutations(range(self._aars_mutation_matrix.shape[0]))
            trnas_perm = permutations(range(self._trna_mutation_matrix.shape[0]))

            possible_codes = tuple([Code(trnas, aarss, self._initial_code._trna_space,
                                         self._initial_code._aars_space)
                                    for aarss in aarss_perm for trnas in trnas_perm])

            self._mutation_probabilities = {to_code:self.mutation_probability(to_code)
                                            for to_code in possible_codes}

            print len(self._mutation_probabilities.values())
            assert np.isclose(1.0, sum(self._mutation_probabilities.values())), \
                "All mutation probabilities must sum up to 1."

        return self._mutation_probabilities

    def mutation_probability(self, to):
        mu_prob = 0.0

        for from_aars, to_aars in zip(self._initial_code.aarss, to.aarss):
            mu_prob += np.log(self._aars_mutation_matrix[from_aars][to_aars])

        for from_trna, to_trna in zip(self._initial_code.trnas, to.trnas):
            mu_prob += np.log(self._trna_mutation_matrix[from_trna][to_trna])

        mu_prob = np.exp(mu_prob)

        assert 0 <= mu_prob <= 1, "Mutation probability was not in the range [0,1]."

        return mu_prob
