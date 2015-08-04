import numpy as np

class Messages(object):
    def __init__(self, coding_matrix, fitness_matrix, message_mutation_matrix):
        if coding_matrix.shape[1] != fitness_matrix.shape[1]:
            raise ValueError("There must be as many amino acids in the coding matrix as there are "
                             "in the fitness matrix.")

        if message_mutation_matrix.shape[0] != coding_matrix.shape[0] or \
           message_mutation_matrix.shape[1] != coding_matrix.shape[0]:
            raise ValueError("The message mutation matrix must be square and have as many"
                             "codons as the coding matrix.")

        if not np.allclose(1.0, coding_matrix.sum(axis=0)):
            raise ValueError("All rows in the coding matrix must sum to 1.")

        if not np.allclose(1.0, message_mutation_matrix.sum(axis=0)):
            raise ValueError("All rows in the message mutation matrix must sum to 1.")

        self._code = coding_matrix
        self._fitness_matrix = fitness_matrix
        self._message_mutation_matrix = message_mutation_matrix

        self._last_codon_usage = None
        self._last_fitness_contributions = None

    @property
    def codon_usage(self):
        if self._last_codon_usage is None:
            self.calculate_at_equilibrium()

        return self._last_codon_usage

    @property 
    def fitness_contributions(self):
        if self._last_fitness_contributions is None:
            self.calculate_at_equilibrium()

        return self._last_fitness_contributions

    def calculate_at_equilibrium(self):
        num_sites = self._fitness_matrix.shape[0]
        num_codons = self._code_matrix.shape[0]
        self._last_codon_usage = np.zeros((num_sites, num_codons))
        self._last_fitness_contributions = np.zeros((num_sites, 1))

        for site in xrange(num_sites):
            mutation_selection_matrix = self._mutation_selection_matrix(site)

            eig_values, eig_vectors = np.linalg.eig(mutation_selection_matrix)
            maxi = np.argmax(eig_values)
            eig_vector = eig_vectors[:,maxi]
            eig_vector /= eig_vector.sum()
            
            self._last_codon_usage[site] = eig_vector
            self._last_fitness_contributions[site] = eig_values[maxi]

        self._last_codon_usage.setflags(write=False)
        self._last_fitness_contributions.setflags(write=False)

    def calculate_at_codon_usage(self, codon_usage):
        num_sites = self._fitness_matrix.shape[0]
        num_codons = self._code_matrix.shape[0]

        if codon_usage.shape != (num_sites, num_codons):
            raise ValueError("The codon usage matrix does not match the coding matrix "
                             "or fitness matrix.")

        self._last_codon_usage = codon_usage
        self._last_fitness_contributions = np.zeros((num_sites, 1))

        for site in xrange(num_sites):
            contribution = 0.0
            for beta in xrange(num_sites):
                for j in xrange(num_sites):
                    contribution += fitness_matrix[site, beta] * code[j, beta] * codon_usage[site, j]
            fitness_contributions[site] = contribution

        self._last_codon_usage.setflags(write=False)
        self._last_fitness_contributions.setflags(write=False)

    def _mutation_selection_matrix(self, site):
        """ Eq4 of SellaArdell06 """
        selection_matrix = np.diag(self._code.dot(self._fitness_matrix[site]))
        return self._message_mutation_matrix.dot(selection_matrix)
