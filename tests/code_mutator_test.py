import unittest
import numpy as np
from code_mutator import Code_Mutator
from code import Code

def get_uniform_mutation_matrix(size, mu):
    mut_matrix = np.diag([1 - ((size-1) * mu)] * size)

    for i in xrange(size):
        for j in xrange(size):
            if i != j:
                mut_matrix[i][j] = mu

    return mut_matrix

class Test_Space(object): # Serves as both tRNA and AARS Spaces
    def __init__(self, mutation_matrix):
        self._mutation_matrix = mutation_matrix
        self._names = map(str, range(mutation_matrix.shape[0]))

    @property
    def mutation_matrix(self):
        return self._mutation_matrix

    @property
    def count(self):
        return self._mutation_matrix.shape[0]

    @property
    def names(self):
        return self._names

class Code_Mutator_Test(unittest.TestCase):
    def setUp(self):
        self.trna_space = Test_Space(get_uniform_mutation_matrix(3, .0001))
        self.aars_space = Test_Space(get_uniform_mutation_matrix(3, .0001))

        self.code = Code([1,2], [1,2], self.trna_space, self.aars_space)

        self.mutator = Code_Mutator(self.code, self.aars_space.mutation_matrix,
                                    self.trna_space.mutation_matrix)


    def test_mutant_codes(self):
        """ Tests if all possible code variants are created. """
        possible_codes = self.mutator.get_possible_codes()

        # Each code contains 4 items with 3 variants for each item
        # variants can be reused
        # so 3^4 = 81
        self.assertEqual(81, len(possible_codes))

    def test_mutant_probability(self):
        """ Does basic checks on the mutation probability from one code to another. """
        # No mutations
        to_code = Code([1,2], [1,2], self.trna_space, self.aars_space)

        self.assertAlmostEqual(.9998 ** 4, self.mutator.mutation_probability(to_code))

        # 1 mutation
        to_code = Code([1,1], [1,2], self.trna_space, self.aars_space)

        self.assertAlmostEqual(.9998 ** 3 * .0001, self.mutator.mutation_probability(to_code))

        # 4 mutations
        to_code = Code([0,0], [0,0], self.trna_space, self.aars_space)

        self.assertAlmostEqual(.0001 ** 4, self.mutator.mutation_probability(to_code))

    def test_custom_mutant_probability(self):
        """ Checks if the mutation probability uses the correct indices. """
        aars_mutation = self.aars_space.mutation_matrix.copy()

        aars_mutation[0,:] = [.75, .15, .10]

        from_code = Code([1,2], [0,1], self.trna_space, self.aars_space)

        custom_mutator = Code_Mutator(from_code, aars_mutation,
                                    self.trna_space.mutation_matrix)

        to_code = Code([1,2], [1,1], self.trna_space, self.aars_space)

        self.assertAlmostEqual(.9998 ** 3 * .15, custom_mutator.mutation_probability(to_code))


if __name__ == '__main__':
    unittest.main()
