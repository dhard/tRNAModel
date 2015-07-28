import numpy as np

class Messages(object):
    def __init__(self, code, site_types, message_mutation_matrix):
        self._site_types = site_types

        self._code = code
        self._message_mutation_matrix = message_mutation_matrix
        
        self._equilibrium_codon_usage = None
        self._equilibrium_fitness_contributions = None

        self._num_sites = len(site_types.sites)
        self._num_codons = code.code_matrix.shape[0]

    def get_equilibrium_codon_usage(self):
        if self._equilibrium_codon_usage is None:
            self._calculate_equilibrium_codon_usage()

        return self._equilibrium_codon_usage.copy()

    def get_equilibrium_fitness_contributions(self):
        if self._equilibrium_fitness_contributions is None:
            self._calculate_equilibrium_codon_usage()

        return self._equilibrium_fitness_contributions.copy()

    def get_mutation_selection_matrix(self, site):
        fitness_matrix = self._site_types.fitness_matrix
        code = self._code.effective_matrix

        # Eq4 of SellaArdell06
        selection_matrix = np.diag(code.dot(fitness_matrix[site]))
        return self._message_mutation_matrix.dot(selection_matrix)

    def _calculate_equilibrium_codon_usage(self):
        effective_code = self._code.effective_matrix
        num_codons = effective_code.shape[1]

        self._equilibrium_codon_usage = np.zeros((self._num_sites, self._num_codons))
        self._equilibrium_fitness_contributions = np.zeros((self._num_sites, 1))

        for site in xrange(self._num_sites):
            # Gets the eigenvector u_{c} and eigenvalue lambda of the site given by the index
            # from Eq3 of SellaArdell06 and stores them in the equilibrium variables
            mutation_selection_matrix = self.get_mutation_selection_matrix(site)

            eig_values, eig_vectors = np.linalg.eig(mutation_selection_matrix)
            maxi = np.argmax(eig_values)
            eig_vector = eig_vectors[:,maxi]
            eig_vector /= eig_vector.sum()
            
            self._equilibrium_codon_usage[site] = eig_vector
            self._equilibrium_fitness_contributions[site] = eig_values[maxi]
        
    def fitness_contributions_at_codon_usage(self, codon_usage):
        fitness_matrix = self._site_types.fitness_matrix
        code_matrix = self._code.effective_matrix
        
        fitness_contributions = np.zeroes(self._num_sites)
        for site in xrange(self._num_sites):
            contribution = 0.0
            for beta in xrange(self._num_sites):
                for j in xrange(self._num_sites):
                    contribution += fitness_matrix[site, beta] * code[j, beta] * codon_usage[site, j]
            fitness_contributions[site] = contribution

        return fitness_contributions
