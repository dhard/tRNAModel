import numpy as np

# TODO incorporate misreading

class Code(object):
    def __init__(self, trnas, aarss, trna_space, aars_space):
        self._trnas = tuple(trnas)
        self._aarss = tuple(aarss)

        self._encoded_aarss = frozenset(aarss)
        self._encoded_trnas = frozenset(trnas)

        possible_aarss = frozenset(range(aars_space.count))
        possible_trnas = frozenset(range(trna_space.count))
        
        self._unencoded_aarss = frozenset(possible_aarss - self._encoded_aarss)
        self._unencoded_trnas = frozenset(possible_trnas - self._encoded_trnas)

        if not self._encoded_aarss.issubset(possible_aarss):
            raise ValueError("All encoded AARSs must exist in the AARS space.")

        if not self._encoded_trnas.issubset(possible_trnas):
            raise ValueError("All encoded tRNAs must exist in the tRNA space.")

        self._code_matrix = None
        self._effective_code_matrix = None
        self._aars_space = aars_space
        self._trna_space = trna_space

    @property
    def code_matrix(self):
        if self._code_matrix is None:
            self._build_code_matrix()

        return self._code_matrix

    @property
    def effective_code_matrix(self):
        if self._effective_code_matrix is None:
            self._build_code_matrix()

        return self._effective_code_matrix

    @property
    def encoded_aarss(self):
        return self._encoded_aarss

    @property
    def encoded_trnas(self):
        return self._encoded_trnas

    @property
    def unencoded_aarss(self):
        return self._unencoded_aarss

    @property
    def unencoded_trnas(self):
        return self._unencoded_trnas

    @property
    def aarss(self):
        return self._aarss

    @property
    def trnas(self):
        return self._trnas

    def _build_code_matrix(self):
        codon_trna_mapping = np.zeros((self._trna_space.codon_trna_mapping.shape[0],
                                       len(self._trnas)))
        for i in xrange(len(self._trnas)):
            codon_trna_mapping[:,i] = self._trna_space.codon_trna_mapping[:,self._trnas[i]]

        trna_aars_mapping = np.zeros((len(self._trnas), len(self._aarss)))

        for i in xrange(len(self._aarss)):
            trna_aars_mapping[:,i] = self._trna_space.trna_aars_mapping[:,self._aarss[i]]

        aars_aa_mapping = np.zeros((len(self._aarss), self._aars_space.aars_aa_mapping.shape[1]))

        for i in xrange(len(self._aarss)):
            aars_aa_mapping[i,:] = self._aars_space.aars_aa_mapping[self._aarss[i],:]

        self._code_matrix = np.dot(np.dot(codon_trna_mapping,
                                          trna_aars_mapping),
                                   aars_aa_mapping)

        # TODO Implement misreading
        self._effective_code_matrix = self._code_matrix

        assert self.code_matrix.shape == (self._trna_space.codon_trna_mapping.shape[0], \
                                          self._aars_space.aars_aa_mapping.shape[1]) # #codons X #amino acids

        assert np.allclose(1.0, self._coding_matrix.sum(axis=1)), \
            "The coding matrix's rows do not sum up to 1."

        assert np.allclose(1.0, self._effective_coding_matrix.sum(axis=1)), \
            "The effective coding matrix's rows do not sum up to 1."

    def __hash__(self):
        # TODO Investigate if this has any obvious collisions
        return hash(hash(self._aarss) + hash(self._trnas))

    def __eq__(self, other):
        return hash(self) == hash(other)
