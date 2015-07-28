import unittest
import numpy as np
from code import Bit_Code

class Code_Test(unittest.TestCase):
    def test_nonambiguous_code_matrix(self):
        code1 = Bit_Code('11010011', 2, 2, 2, 2, None) 
        code1_matrix = np.array([[0, 1], [.5, .5]])

        self.assertTrue(np.allclose(code1.code_matrix, code1_matrix))
        self.assertTrue(np.allclose(code1.effective_matrix, code1_matrix))

    def test_ambiguous_code_matrix(self):
        code1 = Bit_Code('11111111', 2, 2, 2, 2, None) 
        code1_matrix = np.array([[.5, .5], [.5, .5]])
        code2 = Bit_Code(None, 2, 2, 2, 2, None) # Completely ambiguous code by default

        self.assertTrue(np.allclose(code1.code_matrix, code1_matrix))
        self.assertTrue(np.allclose(code1.effective_matrix, code1_matrix))
        
        self.assertTrue(np.allclose(code2.code_matrix, code1_matrix))
        self.assertTrue(np.allclose(code2.effective_matrix, code1_matrix))

    def test_equals_operator(self):
        case1 = Bit_Code('10101010', 2, 2, 2, 2, None)
        case2 = Bit_Code('10101010', 2, 2, 2, 2, None)
        case3 = Bit_Code('11101011', 2, 2, 2, 2, None)

        self.assertTrue(case1 == case2)
        self.assertFalse(case1 == case3)

    def test_misreading():
        self.assertTrue(False, msg="I need to make this still.")


if __name__ == '__main__':
    unittest.main()
