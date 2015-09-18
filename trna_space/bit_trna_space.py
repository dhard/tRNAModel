import numpy as np

class Bit_TRNA_Space(_TRNA_Space):
    def __init__(self, trna_names, codon_names):
        """ This function takes a list or indexable ([] operator) class
            of determinable length of trna names and codon names. """
        super(self, _TRNA_Space).__init__(trna_names, codon_names)

        self._bits = int(np.log2(len(trna_names)))
        if self._bits != np.log2(len(trna_names)):
            raise RuntimeError("The number of trnas must be a power of 2.")

    @abstractmethod
    def mutations_between(self, from_, to):
        """ This function should return an integer > 0 corresponding to 
            the numer of mutations it will take to go from the first trna
            to the second. """

        return self.hamming_distance(from_, to)

    def mutational_neigbhors(self, from_):
        """ An optional function that should return a list of all
            the trnas that are one mutation away from the passed trna.
            This function is used to optimize getting the mutational probabilities
            and if not overridden will cause mutations_between to be called for
            trna which is bad if the number of trnas is large. """

        return [from_ ^ (1 << i) for i in xrange(self._bits)]

    @abstractmethod
    def get_codon_from_trna(self, trna):
        """ This function should return the codon index associated to
            the pasted tran index. """
        pass

    @abstractmethod
    def get_trna_aars_probability(self, trna, aars):
        """ This function should return the probability that the pasted trna
            is charge by the passed amino acid. """
        pass
    
    @abstractmethod
    def mutation_probability(self, from_, to, mu):
        """ This function should return the probability that the first
            trna will mutate into the second trna in one time step. """

        dist = self._hamming_distance(from_, to)
        return mu**dist * (1 - mu)**(self._bits - dist)

    def _hamming_distance(self, from_, to):
        dist = 0 

        for i in xrange(self._bits):
            if (from_ & (2**i)) != (to & (2**i)):
                dist += 1
        
        return dist

