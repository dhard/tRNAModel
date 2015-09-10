import unittest
import numpy as np
from site_types import Ring_Site_Types

"""
From CMCPy SellaArdell02Fig3 demo

codons = cmcpy.codon_spaces.RingCodonSpace(num_codons = 5, mu = 0.01)
aas = cmcpy.amino_acid_spaces.RingAminoAcidSpace(coords = [0,0.2,0.4,0.6,0.8])
site_types = cmcpy.site_type_spaces.MirroringSiteTypeSpace(amino_acids = aas,
                                                           phi = 0.32768,
                                                           weights = [20,20,20,20,20])

correct_fitness_matrix = site_types.get_fitness_matrix()
correct_distance_matrix = aas.get_distance_matrix()
"""

class Ring_Site_Type_Test(unittest.TestCase):
    def setUp(self):
        # From CMCPy SellaArdell02Fig3 demo
        phi = .32768
        vals = [0,0.2,0.4,0.6,0.8]

        aas = [("aa_" + str(i), val) for i, val in enumerate(vals)]
        sites = [("site_" + str(i), val) for i, val in enumerate(vals)]

        weights = [1] * len(vals)

        self.site_type = Ring_Site_Types(phi, aas, sites, weights)

    def test_distance_matrix(self):
       correct_distance_matrix = np.array([[ 0. ,  0.2,  0.4,  0.4,  0.2],
                                           [ 0.2,  0. ,  0.2,  0.4,  0.4],
                                           [ 0.4,  0.2,  0. ,  0.2,  0.4],
                                           [ 0.4,  0.4,  0.2,  0. ,  0.2],
                                           [ 0.2,  0.4,  0.4,  0.2,  0. ]])
        
       self.assertTrue(np.allclose(correct_distance_matrix, self.site_type.distance_matrix))

    def test_fitness_matrix(self):
        correct_fitness_matrix = np.array([[ 1.  ,  0.8 ,  0.64,  0.64,  0.8 ],
                                           [ 0.8 ,  1.  ,  0.8 ,  0.64,  0.64],
                                           [ 0.64,  0.8 ,  1.  ,  0.8 ,  0.64],
                                           [ 0.64,  0.64,  0.8 ,  1.  ,  0.8 ],
                                           [ 0.8 ,  0.64,  0.64,  0.8 ,  1.  ]])
        
        self.assertTrue(np.allclose(correct_fitness_matrix, self.site_type.fitness_matrix))


if __name__ == '__main__':
    unittest.main()
