from abc import ABCMeta, abstractmethod
import numpy as np
from code import Bit_Code

class Code_Mutator(object):
    __metaclass__ = ABCMeta

    def __init__(self, initial_code, mu):
        self._initial_code = initial_code
        self._mu = mu

    @abstractmethod
    def get_all_code_mutants(self):
        pass

    @abstractmethod
    def get_mutation_probability(self, start_code, end_code):
        pass

    @abstractmethod
    def build_code(self, genetic_code):
        pass

    def get_mutation_probabilities(self):
        return {code:get_mutation_probability(self._inital_code, code) for code in get_all_code_mutants()}

class Bit_Code_Mutator(Code_Mutator):
    def __init__(self, initial_code, mu):
        super(Bit_Code_Mutator, self).__init__(initial_code, mu)

    def get_all_code_mutants(self):
        codes = []
        code_length = len(initial_code.tRNAs[0].sequence) * len(initial_code.tRNAs) + \
                      len(initial_code.AARSs[0].sequence) * len(initial_code.AARSs)

        for numerical_code in xrange(2**code_length):
            codes.append(bin(numberical_code)[2:])

        return tuple(codes)

    def get_mutation_probability(self, start_code, end_code):
        code_length = len(start_code)
        differences = sum([i == j for i, j in zip(start_code, end_code)])

        # mutation rate = (1-statewise mutation rate)^(# of similarities) * (statewise mutation rate)^(# of difference)
        # == (# of similarities)np.log(1-statewise mutation rate) * (# of differences)np.log(statewise mutation rate)
        # == np.log((# of similarities)np.log(1-statewise mutation rate)) + np.log((# of differences)np.log(statewise mutation rate))
        # np.exp(np.exp(answer) = mutation rate from start to end

        mu_rate = 0
        if differences > 0:
            mu_rate = np.log(differences * self._mu)
        if (code_length - differences) > 0:
            mu_rate += np.log((code_length - differences) * (1 - self._mu))

        mu_rate = np.exp(np.exp(mu_rate))

        print mu_rate

        assert 0 <= mu_rate <= 1, "Error calculating mutation rate, not in 0 to 1 range."

        return mu_rate  

    def build_code(self, genetic_code):
        return Bit_Code(genetic_code, len(self._initial_code.tRNAs), len(self._initial_code.tRNAs[0].sequence),
                        len(self._initial_code.AARSs), len(self._initial_code.AARSs[0].sequence),
                        self._initial_code.misreading_matrix)
        
