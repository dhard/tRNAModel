import unittest
import numpy as np
from random import Random 
from evolver import Evolver
from code import Code
from trna_space import TRNA_Space
from aars_space import AARS_Space
from site_types import Ring_Site_Types

def get_uniform_mutation_matrix(size, mu):
    mut_matrix = np.diag([1 - ((size-1) * mu)] * size)

    for i in xrange(size):
        for j in xrange(size):
            if i != j:
                mut_matrix[i][j] = mu

    return mut_matrix

correct_nonequilibrium_fitness = 0.042617989729540708

class Evolver_Test(unittest.TestCase):
    def setUp(self):
        self.aars_space = AARS_Space(range(2), np.diag([1] * 2),
                                     get_uniform_mutation_matrix(2, .01))

        # mu is zero so that there is only 3 potential codes available
        # and so the code is forced to handle mutants with a fixation
        # probability of zero.
        self.trna_space = TRNA_Space(range(2), np.diag([1] * 2), np.diag([1] * 2),
                                     get_uniform_mutation_matrix(2, .01))

        pchem_vals = [.30, .70]
        self.site_types = Ring_Site_Types(.05, zip(range(2), pchem_vals),
                                          zip(range(2), pchem_vals), [1] * 2)

        message_mutation_matrix = get_uniform_mutation_matrix(2, .01)
        #print message_mutation_matrix.shape
        #print np.matrix([[.5, .5], [.5, .5]]).shape

        self.code = Code([1], [1], self.trna_space, self.aars_space)

        rng = Random()
        rng.seed(42)
        self.evolver = Evolver(self.code, self.site_types, message_mutation_matrix,
                               100, rng)

        
    def test_overall_fitness(self):
        """ Check if Eq5 SellaArdell 06 is implemented correctly. """
        self.assertAlmostEqual(.5 ** 2, self.evolver._overall_fitness([.5, .5]))

    def test_transition_probabilities(self):
        self.assertTrue(False, "Still need to implement this test.")

    def test_step_time(self):
        self.evolver.step_time()
        self.assertEqual(self.evolver._current_code, self.code)


if __name__ == '__main__':
    unittest.main()
