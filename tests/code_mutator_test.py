import unittest
import numpy as np
from code_mutator import Bit_Code_Mutator
from code import Bit_Code

class Bit_Code_Mutator_Test(unittest.TestCase):
    def test_build_code(self):
        code = Bit_Code('11111111', 2, 2, 2, 2, None) 
        mutator = Bit_Code_Mutator(code, .5)

        self.assertTrue(code == mutator.build_code(code.genetic_code))

    def test_mutation_probability(self):
        code1 = Bit_Code('11111111', 2, 2, 2, 2, None)
        code2 = Bit_Code('00000000', 2, 2, 2, 2, None)
        code3 = Bit_Code('11110000', 2, 2, 2, 2, None) 
        code4 = Bit_Code('11110001', 2, 2, 2, 2, None) 
        mutator = Bit_Code_Mutator(code1, .25)

        self.assertAlmostEqual(2**8 * .25, mutator.get_mutation_probability(code1.genetic_code,
                                                                           code2.genetic_code))

        self.assertAlmostEqual(2**4 * .25 * 2**4 * .75, mutator.get_mutation_probability(code1.genetic_code,
                                                                                         code3.genetic_code))
        
        self.assertAlmostEqual(2**5 * .25 * 2**3 * .75, mutator.get_mutation_probability(code1.genetic_code,
                                                                                         code3.genetic_code))
                


if __name__ == '__main__':
    unittest.main()
