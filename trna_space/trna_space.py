import numpy as np

class _TRNA_Space(object):
    def __init__(self, trna_names, codon_names):
        """ This function takes a list or indexable ([] operator) class
            of determinable length of trna names and codon names. """
        self._trnas = trna_names
        self._codons = codon_names

    @abstractmethod
    def mutations_between(self, from_, to):
        """ This function should return an integer > 0 corresponding to 
            the numer of mutations it will take to go from the first trna
            to the second. """
        pass

    def mutational_neigbhors(self, from_):
        """ An optional function that should return a list of all
            the trnas that are one mutation away from the passed trna.
            This function is used to optimize getting the mutational probabilities
            and if not overridden will cause mutations_between to be called for
            trna which is bad if the number of trnas is large. """
        pass

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
        pass


    def get_trna_aars_map(self, trnas, aarss):
        """ This function return a probability matrix giving the probability
            a trna is charged by an aars. """
        trna_aars_mapping = np.array([[self.get_trna_aars_probability(trna, aars)
                                       for aars in aarss]
                                      for trna in trnas], dtype=np.float)

        row_sums = np.array([[row_sum] if row_sum > 0 else [1.]
                             for row_sum in trna_aars_mapping.sum(axis=1)])

        return trna_aars_mapping / row_sums

    def get_codon_trna_map(self, trnas):
        """ This function returns a probability matrix giving the probability 
            that a codon is read by a trna. """
        codon_trna_mapping = np.zeros((self.codons, len(trnas)), dtype=np.float)

        for i, trna in enumerate(trnas):
            codon = self.get_codon_from_trna(trna)
            codon_trna_mapping[codon, i] = 1.0

        # Avoid divide by zero
        row_sums = np.array([[row_sum] if row_sum > 0 else [1.]
                             for row_sum in codon_trna_mapping.sum(axis=1)])

        return codon_trna_mapping / row_sums

    @property
    def names(self):
        return self._trnas

    @property
    def codon_names(self):
        return self._codons

    @property
    def count(self):
        return len(self._trnas)

    @property
    def codons(self):
        return len(self._codons)
