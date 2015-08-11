import unittest
import numpy as np
from code import Code
from trna_space import TRNA_Space
from aars_space import AARS_Space


class Code_Test(unittest.TestCase):
    def test_cmcpy_comparison(self):
        """ Checks if the coding matrix generated is the same as
            the one cmcpy generates. """

        # cmcpy equivalent code
        # aas = cmcpy.amino_acid_spaces.RegionAminoAcidSpace(coords = [0.0085,0.0086,0.0886,0.0988,0.1005,0.1363,0.2879,0.3254,0.3424,0.3425,0.3816,0.3817,0.4497,0.5213,0.5963,0.6048,0.7513,0.8637,0.9608,0.9659])
        # codons = cmcpy.codon_spaces.WordCodonSpace(num_bases = 4, num_positions = 2, mu = .01)
        # initial_code = cmcpy.genetic_codes.InitiallyAmbiguousGeneticCode(codons = codons,amino_acids = aas)
        # initial_code.get_code_matrix()
        # Which is just
        # numpy.ones((nc,na)) / na
        # or in this case
        # numpy.ones((16,20)) / 20

        trna_space = TRNA_Space(range(20), np.diag([1] * 20), np.array([[1. / 20] * 20] * 16), np.diag([1] * 20))
        aars_space = AARS_Space(range(20), np.diag([1] * 20), np.diag([1] * 20))
        code = Code([1] * 16, range(20), trna_space, aars_space)

        self.assertTrue(np.allclose(code.code_matrix, np.ones((16,20)) / 20))

if __name__ == '__main__':
    unittest.main()
