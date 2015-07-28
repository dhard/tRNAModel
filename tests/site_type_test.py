import unittest
import numpy as np
from site_types import Ring_Site_Types

class Ring_Site_Type_Test(unittest.TestCase):
    def setUp(self):
        self.phi = .5
        self.site_type = Ring_Site_Types(self.phi, [0.156, 0.375, 0.599, 0.732, 0.951], [1] * 5)

    def test_distance_matrix(self):
        correct_matrix = np.array([[ 0.   ,  0.219,  0.443,  0.424,  0.205],
                                   [ 0.219,  0.   ,  0.224,  0.357,  0.424],
                                   [ 0.443,  0.224,  0.   ,  0.133,  0.352],
                                   [ 0.424,  0.357,  0.133,  0.   ,  0.219],
                                   [ 0.205,  0.424,  0.352,  0.219,  0.   ]])
        self.assertTrue(np.allclose(correct_matrix, self.site_type.distance_matrix))

    def test_fitness_matrix(self):
        correct_matrix = np.array([[ 0.   ,  0.219,  0.443,  0.424,  0.205],
                                   [ 0.219,  0.   ,  0.224,  0.357,  0.424],
                                   [ 0.443,  0.224,  0.   ,  0.133,  0.352],
                                   [ 0.424,  0.357,  0.133,  0.   ,  0.219],
                                   [ 0.205,  0.424,  0.352,  0.219,  0.   ]])
        correct_matrix = self.phi ** correct_matrix
        self.assertTrue(np.allclose(correct_matrix, self.site_type.fitness_matrix))


if __name__ == '__main__':
    unittest.main()
