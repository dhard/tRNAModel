import numpy as np
from itertools import product
from code import Code

class Code_Mutator(object):
    def __init__(self, initial_code, mu):
        if not 0 <= mu <= 1:
            raise RuntimeError("Mu must be in the range [0,1].")

        self._mu = mu
        self._initial_code = initial_code

        # Set values to None so they can be built as needed
        self._possible_codes = None
        self._mutation_probabilities = None

        self._one_gene_possible_codes = None
        self._one_gene_mutation_probabilities = None

        self._trna_neighbors = initial_code.trna_space.mutational_neighbors(initial_code.trnas[0]) != None
        self._aars_neighbors = initial_code.aars_space.mutational_neighbors(initial_code.aarss[0]) != None

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
            trna_space = self._initial_code.trna_space
            aars_space = self._initial_code.aars_space

            self._one_gene_possible_codes = []

            if self._trna_neighbors:
                for i, trna in enumerate(self._initial_code.trnas):
                    neighbors = trna_space.mutational_neighbors(trna)
                    for new_trna in neighbors:
                        trnas = list(self._initial_code.trnas)
                        trnas[i] = new_trna
                        self._one_gene_possible_codes.append(Code(trnas, self._initial_code.aarss, trna_space,
                                                                  aars_space))
            else:
                for i, trna in enumerate(self._initial_code.trnas):
                    for new_trna in xrange(trna_space.count):
                        if trna_space.mutations_between(trna, new_trna) == 1:
                            trnas = list(self._initial_code.trnas)
                            trnas[i] = new_trna
                            self._one_gene_possible_codes.append(Code(trnas, self._initial_code.aarss, trna_space,
                                                                      aars_space))

            if self._aars_neighbors:
                for i, aars in enumerate(self._initial_code.aarss):
                    neighbors = aars_space.mutational_neighbors(aars)
                    for new_aars in neighbors:
                        aarss = list(self._initial_code.aarss)
                        aarss[i] = new_aars
                        self._one_gene_possible_codes.append(Code(self._initial_code.trnas, aarss, trna_space,
                                                                  aars_space))
            else:
                for i, aars in enumerate(self._initial_code.aarss):
                    for new_aars in xrange(aars_space.count):
                        if aars_space.mutations_between(aars, new_aars) == 1:
                            aarss = list(self._initial_code.aarss)
                            aarss[i] = new_aars
                            self._one_gene_possible_codes.append(Code(self._initial_code.trnas, aarss, trna_space,
                                                                      aars_space))

        return self._one_gene_possible_codes

    def mutation_probability(self, to):
        trna_space = self._initial_code.trna_space
        aars_space = self._initial_code.aars_space

        mu_prob = 1.0

        for from_aars, to_aars in zip(self._initial_code.aarss, to.aarss):
            prob = aars_space.mutation_probability(from_aars, to_aars, self._mu)
            mu_prob *= prob

        if mu_prob == 0.0:
            return 0.0

        for from_trna, to_trna in zip(self._initial_code.trnas, to.trnas):
            prob = trna_space.mutation_probability(from_trna, to_trna, self._mu)
            mu_prob *= prob

        if mu_prob == 0.0:
            return 0.0

        assert 0 <= mu_prob <= 1, "Mutation probability was not in the range [0,1]."

        return mu_prob
