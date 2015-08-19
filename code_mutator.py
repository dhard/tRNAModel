import numpy as np
from itertools import product
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

        self._possible_codes = None
        self._mutation_probabilities = None

        self._one_gene_possible_codes = None
        self._one_gene_mutation_probabilities = None

    def get_possible_codes(self):
        if self._possible_codes is None:
            num_aarss = self._aars_mutation_matrix.shape[0]
            num_trnas = self._trna_mutation_matrix.shape[0]
            encodable_trnas = len(self._initial_code.trnas)
            encodable_aarss = len(self._initial_code.aarss)
            aarss_prod = list(product(range(num_aarss), repeat=encodable_aarss))
            trnas_prod = list(product(range(num_trnas), repeat=encodable_trnas))

            self._possible_codes = [Code(trnas, aarss, self._initial_code._trna_space, 
                                         self._initial_code._aars_space)
                                    for trnas in trnas_prod for aarss in aarss_prod]

            potential_codes = (num_trnas ** encodable_trnas) * (num_aarss ** encodable_aarss)
                              
            assert len(self._possible_codes) == potential_codes, \
                "An incorrect number of potential codes was generated {} {}.".format(potential_codes, len(self._possible_codes))

        return self._possible_codes

    def get_mutation_probabilities(self):
        if self._mutation_probabilities is None:
            self._mutation_probabilities = {to_code:self.mutation_probability(to_code)
                                            for to_code in self.get_possible_codes()}

            assert np.isclose(1.0, sum(self._mutation_probabilities.values())), \
                "All mutation probabilities must sum up to 1."

        return self._mutation_probabilities

    def get_one_gene_mutation_probabilities(self):
        if self._one_gene_mutation_probabilities is None:
            self._one_gene_mutation_probabilities = {to_code:self.mutation_probability(to_code)
                                                     for to_code in self.get_one_gene_possible_codes()}

            self._one_gene_mutation_probabilities[self._initial_code] = 1 - sum(self._one_gene_mutation_probabilities.values())

        return self._one_gene_mutation_probabilities

    def get_one_gene_possible_codes(self):
        if self._one_gene_possible_codes is None:
            trna_space = self._initial_code._trna_space
            aars_space = self._initial_code._aars_space

            self._one_gene_possible_codes = []

            for i, trna in enumerate(self._initial_code.trnas):
                for new_trna in xrange(trna_space.count):
                    if trna_space.mutations_between(trna, new_trna) == 1:
                        trnas = list(self._initial_code.trnas)
                        trnas[i] = new_trna
                        self._one_gene_possible_codes.append(Code(trnas, self._initial_code.aarss, trna_space,
                                                                  aars_space))

            for i, aars in enumerate(self._initial_code.aarss):
                for new_aars in xrange(aars_space.count):
                    if aars_space.mutations_between(aars, new_aars) == 1:
                        aarss = list(self._initial_code.aarss)
                        aarss[i] = new_aars
                        self._one_gene_possible_codes.append(Code(self._initial_code.trnas, aarss, trna_space,
                                                                  aars_space))

        return self._one_gene_possible_codes

    def mutation_probability(self, to):
        mu_prob = 0.0

        for from_aars, to_aars in zip(self._initial_code.aarss, to.aarss):
            if not np.isclose(0, self._aars_mutation_matrix[from_aars][to_aars]):
                mu_prob += np.log(self._aars_mutation_matrix[from_aars][to_aars])
            else:
                mu_prob = 0.0
                break

        if mu_prob == 0.0:
            return 0

        for from_trna, to_trna in zip(self._initial_code.trnas, to.trnas):
            if not np.isclose(0, self._trna_mutation_matrix[from_trna][to_trna]):
                mu_prob += np.log(self._trna_mutation_matrix[from_trna][to_trna])
            else:
                mu_prob = 0.0
                break

        if mu_prob == 0.0:
            return 0

        mu_prob = np.exp(mu_prob)

        assert 0 <= mu_prob <= 1, "Mutation probability was not in the range [0,1]."

        return mu_prob
