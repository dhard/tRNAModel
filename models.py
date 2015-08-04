import numpy as np
from abc import ABCMeta, abstractmethod

from site_types import Ring_Site_Types
from aars_space import AARS_Space
from trna_space import TRNA_Space
from code import Code

def get_ring_mutation_matrix(size, mu):
    mut_matrix = np.diag([1 - (2 * mu)] * size)

    for i in xrange(size):
        mut_matrix[i][i - 1] = mu
        mut_matrix[i][i + 1 if (i + 1) < size else 0 ] = mu # Handles wraping around to 0

    return mut_matrix

def get_uniform_mutation_matrix(size, mu):
    mut_matrix = np.diag([1 - ((size-1) * mu)] * size)

    for i in xrange(size):
        for j in xrange(size):
            if i != j:
                mut_matrix[i][j] = mu

    return mut_matrix

class _Model(object):
    __metaclass__ = ABCMeta

    def __init__(self, args, rng):
        pass

    @abstractmethod
    def get_site_types(self):
        pass

    @abstractmethod
    def get_message_mutation_matrix(self):
        pass

    @abstractmethod
    def get_initial_code(self):
        pass


class Bit_Model(_Model):
    def __init__(self, args, rng):
        super(Bit_Model, self).__init__(self, args)

        self.trnas = 2**args.trna_length
        self.aarss = 2**args.aars_length

        trna_ids = ["trna_" + str(i) for i in xrange(self.trnas)]
        aars_ids = ["aars_" + str(i) for i in xrange(self.aarss)]

        aas_ids = ["aas_" + str(i) for i in xrange(args.amino_acids)]
        aa_vals = sorted([rng.random() for _ in xrange(args.amino_acids)])

        site_ids = ["site_" + str(i) for i in xrange(args.site_types)]
        site_vals = sorted([rng.random() for _ in xrange(args.site_types)])

        weights = args.site_weights or ['1'] * len(site_ids)
        weights = map(int, list(weights))

        self._site_types = Ring_Site_Types(args.phi, zip(aas_ids, aa_vals),
                                           zip(site_ids, site_vals), weights)

        self._message_mutation_matrix = get_ring_mutation_matrix(len(site_ids), args.message_mu)

        assert len(aa_vals) == len(aars_ids), "The number of AARSs and amino acids must be equal."

        aars_aa_map = np.diag([1] * len(aa_vals))
        aars_mutation_matrix = get_uniform_mutation_matrix(len(aars_ids), args.message_mu)
        
        aars_space = AARS_Space(aars_ids, aars_aa_map, aars_mutation_matrix)

        self.max_len = args.trna_length
        trna_aars_mapping = np.array([[len(trna) - self._hamming_distance(trna, aars)
                                       for aars in map(self.bin_pad, xrange(len(aars_ids)))]
                                      for trna in map(self.bin_pad, xrange(len(trna_ids)))], dtype=np.float)

        row_sums = trna_aars_mapping.sum(axis=1, keepdims=True)
        trna_aars_mapping = trna_aars_mapping / row_sums

        codon_trna_mapping = np.diag([1] * len(trna_ids))

        trna_mutation_matrix = get_uniform_mutation_matrix(len(trna_ids), args.message_mu)
        trna_space = TRNA_Space(trna_ids, trna_aars_mapping, codon_trna_mapping, trna_mutation_matrix)

        self._initial_code = Code([0] * args.trnas, [0] * args.aarss, trna_space, aars_space)

    def get_site_types(self):
        return self._site_types

    def get_message_mutation_matrix(self):
        return self._message_mutation_matrix

    def get_initial_code(self):
        return self._initial_code

    def _hamming_distance(self, arg1, arg2):
        assert len(arg1) == len(arg2), "Objects must be same length."
        return sum(x != y for x, y in zip(arg1, arg2))

    def bin_pad(self, arg): # TODO Fix this hacky thing
        arg = bin(arg)[2:]

        delta = self.max_len - len(arg)
        arg = ''.join(["0"] * delta) + arg

        return arg

        
