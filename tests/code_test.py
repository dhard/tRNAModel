import unittest
import numpy as np
from code import Code
from trna_space import TRNA_Space
from aars_space import AARS_Space
from messages import Messages
from site_types import Ring_Site_Types

def hamming_dist(from_, to, length):
    dist = 0
    for i in xrange(length):
        if from_ & (2**i) != to & (2**i):
            dist += 1

    return dist

def get_ring_mutation_matrix(size, mu):
    mut_matrix = np.diag([1 - (2 * mu)] * size)

    for i in xrange(size):
        mut_matrix[i][i - 1] = mu
        mut_matrix[i][i + 1 if (i + 1) < size else 0 ] = mu # Handles wraping around to 0

    return mut_matrix

class Code_Test(unittest.TestCase):
    def test_SellaArdell02(self):
        """ This test tries to replicate the results from Sella Ardell 02 Fig 3/4 """

        def id_map(val):
            return val

        def code_one(trna, aars):
            if trna == aars:
                return 1.0

            return 0.0
        
        code2_dict = {0:0, 1:3, 2:1, 3:4, 4:2}
        def code_two(trna, aars):
            if code2_dict[trna] == aars:
                return 1.0

            return 0.0

        def dist(gene1, gene2):
            return hamming_dist(gene1, gene2, 3)

        trnas = 5
        aarss = 5

        msg_mu = .01
        phi = .32768

        weights = [20] * aarss

        pchem_vals = [0.0, 0.2, 0.4, 0.6, 0.8]

        starting_trnas = range(5)
        starting_aarss = range(5)

        trna_ids = ["trna_" + str(i) for i in xrange(trnas)]
        codon_ids = ["codon_" + str(i) for i in xrange(trnas)]
        aars_ids = ["aars_" + str(i) for i in xrange(aarss)]

        aa_ids = ["aas_" + str(i) for i in xrange(aarss)]
        aa_vals = pchem_vals

        site_ids = ["site_" + str(i) for i in xrange(aarss)]
        site_vals = pchem_vals

        trna_space_one = TRNA_Space(range(5), range(5), code_one, id_map, dist, np.diag([1] * 5))
        trna_space_two = TRNA_Space(range(5), range(5), code_two, id_map, dist, np.diag([1] * 5))
        aars_space = AARS_Space(range(5), range(5), id_map, dist, np.diag([1] * 5))

        code_one = Code(starting_trnas, starting_aarss, trna_space_one, aars_space)
        code_two = Code(starting_trnas, starting_aarss, trna_space_two, aars_space)

        site_types = Ring_Site_Types(phi, zip(aa_ids, aa_vals),
                                     zip(site_ids, site_vals), weights)

        msgs_one = Messages(code_one.effective_code_matrix, site_types.fitness_matrix,
                            get_ring_mutation_matrix(trnas, msg_mu))

        msgs_two = Messages(code_two.effective_code_matrix, site_types.fitness_matrix,
                            get_ring_mutation_matrix(trnas, msg_mu))
        
        fit_one = 1
        for val in msgs_one.fitness_contributions:
            fit_one *= val ** 20

        fit_two = 1
        for val in msgs_two.fitness_contributions:
            fit_two *= val ** 20

        code_one_expected_fit = 0.144094106475
        code_two_expected_fit = 0.137619106921

        self.assertAlmostEqual(code_one_expected_fit, fit_one)
        self.assertAlmostEqual(code_two_expected_fit, fit_two)
        

    """def test(self):
        def id_map(val):
            return val
        
        def trna_aars(trna, aars):
            dist = 0
            for i in xrange(3):
                if trna & (2**i) == aars & (2**i):
                    dist += 1

            return dist / 3.
        
        def dist(gene1, gene2):
            return hamming_dist(gene1, gene2, 3)

        trna_space = TRNA_Space(range(5), range(5), trna_aars, id_map, dist, np.diag([1] * 5))
        aars_space = AARS_Space(range(5), range(5), id_map, dist, np.diag([1] * 5))

        code = Code(range(5), range(5), trna_space, aars_space)
        print code.code_matrix

        from code_mutator import Code_Mutator

        mutator = Code_Mutator(code, np.diag([1] * 5), np.diag([1] * 5))

        for code in mutator.get_one_gene_possible_codes():
            print code

        print len(mutator.get_one_gene_possible_codes())"""

if __name__ == '__main__':
    unittest.main()
