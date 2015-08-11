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
        # Grab all the rows/codons currently encoded from the full matrix
        codon_trna_mapping = self._trna_space.codon_trna_mapping[self.trnas,:]

        # The mapping should have the same number of rows as encoded trnas and the same
        # number of columns as potentially encodable trnas
        assert codon_trna_mapping.shape == (len(self.trnas), self._trna_space.codon_trna_mapping.shape[1]), \
            "The codon trna mapping has the wrong size."

        self._code_matrix = np.dot(np.dot(codon_trna_mapping,
                                          self._trna_space.trna_aars_mapping),
                                   self._aars_space.aars_aa_mapping)

        # TODO Implement misreading
        self._effective_code_matrix = self._code_matrix

        # The number of encoded trnas serves as a proxy for the max number
        # of encodable codons, i.e. the number of encodable codons if each 
        # encoded trna encoded a unique codon
        # The number of columns is just the # of amino acids in the system.
        assert self.code_matrix.shape == (len(self.trnas), \
                                          self._aars_space.aars_aa_mapping.shape[1]), \
            "The code matrix is of the wrong dimensions."


        codon_trna_mapping = np.zeros((self._trna_space.codon_trna_mapping.shape[0],
                                       len(self._trnas)))

        assert np.allclose(1.0, self._code_matrix.sum(axis=1)), \
            "The coding matrix's rows do not sum up to 1."

        assert np.allclose(1.0, self._effective_coding_matrix.sum(axis=1)), \
            "The effective coding matrix's rows do not sum up to 1."
        
    def __hash__(self):
        # TODO Investigate if this has any obvious collisions
        return hash(str(self))

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __str__(self):
        return ','.join([self._trna_space.names[trna_id] for trna_id in self._trnas] + \
                        [self._aars_space.names[aars_id] for aars_id in self._aarss])
