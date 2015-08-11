import numpy as np

class TRNA_Space(object):
    def __init__(self, trna_names, trna_aars_mapping, codon_trna_mapping, trna_mutation_matrix):
        self._trnas = len(trna_names)

        if self._trnas != trna_aars_mapping.shape[0]:
            raise ValueError("There must be as many tRNA names as there are tRNAs in the "
                             "tRNA-AARS mapping matrix.")

        if self._trnas != codon_trna_mapping.shape[1]:
            raise ValueError("There must be as many tRNA names as there are tRNAs in the "
                             "codon-tRNA mapping matrix.")

        if self._trnas != trna_mutation_matrix.shape[0] or self._trnas != trna_mutation_matrix.shape[1]:
            raise ValueError("There must be as many tRNA names as there are tRNAs in the "
                             "mutation matrix.")

        if not np.allclose(1.0, codon_trna_mapping.sum(axis=0)):
            raise ValueError("All rows must sum to 1 in the codon-tRNA mapping.")
            
        if not np.allclose(1.0, trna_mutation_matrix.sum(axis=0)):
            raise ValueError("All rows must sum to 1 in the mutation matrix.")

        self._codon_trna_mapping = codon_trna_mapping
        self._trna_aars_mapping = trna_aars_mapping
        self._mutation_matrix = trna_mutation_matrix
        self._names = tuple(map(str, trna_names)) # Force ordering and immutability

        self._codon_trna_mapping.setflags(write=False)
        self._trna_aars_mapping.setflags(write=False)
        self._mutation_matrix.setflags(write=False)

    @property
    def codon_trna_mapping(self):
        return self._codon_trna_mapping

    @property
    def trna_aars_mapping(self):
        return self._trna_aars_mapping

    @property
    def mutation_matrix(self):
        return self._mutation_matrix

    @property
    def names(self):
        return self._names

    @property
    def count(self):
        return self._trnas
