import numpy as np

class TRNA_Space(object):
    def __init__(self, trna_names, codon_names, trna_aars_map, trna_codon_map, mutations_between, trna_mutation_matrix):
        self._trnas = len(trna_names)
        self._codons = len(codon_names)

        #if self._trnas != trna_aars_mapping.shape[0]:
        #    raise ValueError("There must be as many tRNA names as there are tRNAs in the "
        #                     "tRNA-AARS mapping matrix.")

        #if self._trnas != codon_trna_mapping.shape[1]:
        #    raise ValueError("There must be as many tRNA names as there are tRNAs in the "
        #                     "codon-tRNA mapping matrix.")

        if self._trnas != trna_mutation_matrix.shape[0] or self._trnas != trna_mutation_matrix.shape[1]:
            raise ValueError("There must be as many tRNA names as there are tRNAs in the "
                             "adjacency matrix.")

        #if not np.allclose(1.0, codon_trna_mapping.sum(axis=1)):
        #    raise ValueError("All rows must sum to 1 in the codon-tRNA mapping.")
            
        if not np.allclose(1.0, trna_mutation_matrix.sum(axis=1)):
            raise ValueError("All rows must sum to 1 in the mutation matrix.")

        self._trna_codon_map = trna_codon_map
        self._trna_aars_map = trna_aars_map
        self._mutations_between = mutations_between
        self._mutation_matrix = trna_mutation_matrix
        self._names = tuple(map(str, trna_names)) # Force ordering and immutability
        self._codon_names = tuple(map(str, codon_names))

        self._mutation_matrix.setflags(write=False)


    def get_codon_from_trna(self, trna):
        return self._trna_codon_map(trna)

    def get_trna_aars_probability(self, trna, aars):
        prob = self._trna_aars_map(trna, aars)

        assert 0 <= prob <= 1, "tRNA to AARS mapping probability must be in the range [0,1]."

        return prob

    def get_codon_trna_map(self, trnas):
        codon_trna_mapping = np.zeros((self.codons, len(trnas)), dtype=np.float)

        for i, trna in enumerate(trnas):
            codon = self.get_codon_from_trna(trna)
            codon_trna_mapping[codon, i] = 1.0

        # Avoid divide by zero
        row_sums = np.array([[row_sum] if row_sum > 0 else [1.]
                             for row_sum in codon_trna_mapping.sum(axis=1)])

        return codon_trna_mapping / row_sums

    def get_trna_aars_map(self, trnas, aarss):
        trna_aars_mapping = np.array([[self.get_trna_aars_probability(trna, aars)
                                       for aars in aarss]
                                      for trna in trnas], dtype=np.float)

        row_sums = np.array([[row_sum] if row_sum > 0 else [1.]
                             for row_sum in trna_aars_mapping.sum(axis=1)])

        return trna_aars_mapping / row_sums

    def mutations_between(self, trna1, trna2):
        mutations = self._mutations_between(trna1, trna2)

        assert mutations >= 0, "Mutations between two tRNAs must be greater than 0."
        
        return mutations

    @property
    def mutation_matrix(self):
        return self._mutation_matrix

    @property
    def names(self):
        return self._names

    @property
    def codon_names(self):
        return self._codon_names

    @property
    def count(self):
        return self._trnas

    @property
    def codons(self):
        return self._codons
