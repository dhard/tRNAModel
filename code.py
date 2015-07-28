from abc import ABCMeta, abstractmethod
from copy import copy
import numpy as np
from evolvable import Evolvable


class Code(object):
    __metaclass__ = ABCMeta

    def __init__(self, genetic_code, misreading_matrix):
        self._genetic_code = genetic_code
        self._misreading_matrix = misreading_matrix

        self._effective_code_matrix = None
        self._code_matrix = None
        
    @abstractmethod
    def build_code_matrix(self):
        assert np.allclose(self._code_matrix.sum(axis=1), 1), "Code matrix rows did not sum to 1."
        assert np.allclose(self._effective_code_matrix.sum(axis=1), 1), \
            "Effective code matrix rows did not sum to 1."

    @property
    def effective_matrix(self):
        if self._effective_code_matrix is None:
            self.build_code_matrix()
        return self._effective_code_matrix.copy()

    @property
    def code_matrix(self):
        if self._code_matrix is None:
            self.build_code_matrix()
        return self._code_matrix.copy()

    @property
    def misreading_matrix(self):
        return copy(self._misreading_matrix)

    @property
    def genetic_code(self):
        return self._genetic_code

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __hash__(self):
        return hash(self._genetic_code)


class Bit_Code(Code):
    def __init__(self, genetic_code, tRNA_count, tRNA_length, AARS_count, AARS_length, misreading_matrix):
        assert tRNA_count > 1, "There must be more than one tRNA."
        assert tRNA_length > 1, "tRNA length must be greater than one."
        assert AARS_count > 1, "There must be more than one AARS."
        assert AARS_length > 1, "AARS length must be greater than one."

        self._tRNAs = []
        self._AARSs = []
        if genetic_code is not None:
            assert len(genetic_code) == (tRNA_count * tRNA_length + AARS_count * AARS_length), \
                "Inital code must be equal in length to the sum of tRNAs and AARSs."

            read_index = 0
            for _ in xrange(tRNA_count):
                self._tRNAs.append(Evolvable(genetic_code[read_index:read_index + tRNA_length]))
                read_index += tRNA_length

            for _ in xrange(AARS_count):
                self._AARSs.append(Evolvable(genetic_code[read_index:read_index + AARS_length]))
                read_index += AARS_length
        else: # Assume inital code is all zeros
            for _ in xrange(tRNA_count):
                self._tRNAs.append(Evolvable('0' * tRNA_length))

            for _ in xrange(AARS_count):
                self._AARSs.append(Evolvable('0' * AARS_length))

        self._tRNAs = tuple(self._tRNAs)
        self._AARSs = tuple(self._AARSs)

        super(Bit_Code, self).__init__(genetic_code, misreading_matrix)

    def build_code_matrix(self):
        self._code_matrix = np.array([[self.similarities(tRNA, AARS) for AARS in self._AARSs]
                                      for tRNA in self._tRNAs], dtype=np.float)
        row_sums = self._code_matrix.sum(axis=1)

        self._code_matrix = self._code_matrix / row_sums
        if self._misreading_matrix == None:
            self._effective_code_matrix = self._code_matrix
        else:
            self._effective_code_matrix = self._code_matrix.dot(self._misreading_matrix)

        super(Bit_Code, self).build_code_matrix()


    def similarities(self, tRNA, AARS):
        return sum([i == j for i, j in zip(tRNA.sequence, AARS.sequence)])

    @property
    def tRNAs(self):
        return self._tRNAs

    @property
    def AARSs(self):
        return self._AARSs
