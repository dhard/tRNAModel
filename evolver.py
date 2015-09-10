import numpy as np
from messages import Messages
from code_mutator import Code_Mutator

class Evolver(object):
    def __init__(self, initial_code, site_types, message_mutation_matrix, pop_size, rng):
        if pop_size < 1:
            raise ValueError("The population size must be greater than 0.")

        if initial_code._trna_space.codons != message_mutation_matrix.shape[0] or \
           initial_code._trna_space.codons != message_mutation_matrix.shape[1]:
            raise ValueError("The message mutation matrix must be square with {} "
                             "codons per side. It is currently {}".format(initial_code._trna_space.codons,
                                                                          message_mutation_matrix.shape))
            
        
        self._current_code = initial_code
        self._message_mutation_matrix = message_mutation_matrix
        self._fitness_matrix = site_types.fitness_matrix
        self._site_weights = site_types.weights
        self._pop_size = pop_size
        self._rng = rng

        # All these values depend on the current code
        # and are created/destroyed as needed
        self._last_messages = None
        self._current_code_fitness = None
        self._code_mutator = None

        self._frozen = False

    def _code_fitness(self, code, codon_usage = None):
        """ Gets the overall fitness of a code at the given codon usage,
            if no codon usage is given the codon usage of the 
            currently fixed code at equilibrium is used. """
        
        fitness_contributions = None

        if codon_usage is None:
            if self._last_messages is None:
                if code == self._current_code:
                    msgs = Messages(self._current_code.effective_code_matrix,
                                    self._fitness_matrix,
                                    self._message_mutation_matrix)
                    msgs.calculate_at_equilibrium()
                    self._last_messages = msgs
                    fitness_contributions = msgs.fitness_contributions
                else:
                    self._code_fitness(self._current_code)

        if fitness_contributions is None:
            msgs = Messages(code.effective_code_matrix,
                            self._fitness_matrix,
                            self._message_mutation_matrix)
            msgs.calculate_at_codon_usage(self._last_messages.codon_usage)
            fitness_contributions = msgs.fitness_contributions

        overall_fitness = self._overall_fitness(fitness_contributions)

        if codon_usage is None and code == self._current_code:
            self._current_code_fitness = overall_fitness

        return overall_fitness

    def _overall_fitness(self, fitness_contributions):
        """ Eq5 SellaArdell06, uses the fitness contribution of each site
            and the accompanying weights of each site to get an overall
            fitness of the code. """

        return np.prod([fitness_contributions[site]**self._site_weights[site]
                        for site in xrange(len(self._site_weights))])

    def _transition_probability(self, to, mutation_probability):
        """ Returns the probability a code mutates from the current code
            to the other code. This is just Eq3 Sella09. """

        if self._current_code == to: # This case must be handled when the transition matrix is already built
            # It is really 1 - summation(transition_probability(from_, to), 
            #                            for all potential values of end)
            return 0

        if mutation_probability == 0.0:
            return 0

        fix_probability = self._fixation_probability(to)
        if fix_probability == 0.0:
            return 0
        
        trans_probability = self._pop_size * fix_probability * \
                            mutation_probability
        #trans_probability = np.log(self._pop_size) + np.log(fix_probability) + \
        #                    np.log(mutation_probability)
        #trans_probability = np.exp(trans_probability)

        assert 0 <= trans_probability <= 1, "Transition probability was not in the range [0,1]."

        return trans_probability

    def _fixation_probability(self, to):
        """ Returns the proability of the next fixed code being the code passed
            Eq1 Sella09 - Fixation probability of a Moran birth-death process. """

        if self._last_messages is None:
            # sets self._last_messages
            # and self._current_code_fitness
            self._code_fitness(self._current_code)

        to_fitness = self._code_fitness(to, self._last_messages.codon_usage)

        if to_fitness == 0.0:
            # This is the limit of the function below
            # as the "to_fitness" becomes 0
            return 0 

        assert self._current_code_fitness > 0, \
            "The fitness of the current code should always be greater than 0."

        if np.isclose(self._current_code_fitness, to_fitness):
            # This is the limit of the function below
            # as the fitnesses become equal
            return 1.0 / self._pop_size 

        fitness_ratio = self._current_code_fitness / to_fitness

        # Check if calculation will overflow
        if (self._pop_size * np.log10(fitness_ratio)) > 308: # largest number floats can represent ~10^308
            return 0 # Same limit as above for to_fitness

        # Calculates the Moran fixation probability
        num = 1 - fitness_ratio
        dem = 1 - (fitness_ratio)**self._pop_size

        fix_prob = num / dem

        assert 0 <= fix_prob <= 1, "Fixation probability must be in the range [0,1]."
        
        return num / dem

    def get_transition_probabilities(self):
        """ Returns a dictionary of the probability of going from the current code
            to a mutant code in one time step. """

        # Keep things cached in case the code doesn't evolve over one step
        if self._code_mutator is None:
            self._code_mutator = Code_Mutator(self._current_code,
                                              self._current_code._aars_space.mutation_matrix,
                                              self._current_code._trna_space.mutation_matrix)

            mutant_prob = self._code_mutator.get_one_gene_mutation_probabilities()

            self._transition_probabilities = {to_code:self._transition_probability(to_code, mutant_prob[to_code])
                                              for to_code in mutant_prob}

            # Handles the corner case in _transition_probability()
            self._transition_probabilities[self._current_code] = 1.0 - sum(self._transition_probabilities.values())

        return self._transition_probabilities

    def step_time(self):
        """ This function steps one generation allowing the current code to
            potentially mutate into another code. """

        #if self.frozen:
        #    return
        
        transition_probabilities = self.get_transition_probabilities()

        range_bottom = 0
        random_num = self._rng.random()
        code_changed = False

        possible_transitions = 0

        for to_code in transition_probabilities:
            range_top = transition_probabilities[to_code] + range_bottom

            # Avoid infinite loop
            if transition_probabilities[to_code] > 0:
                possible_transitions += 1

            if range_bottom <= random_num <= range_top:
                code_changed = True
                if to_code != self._current_code:
                    self._current_code = to_code

                    # Invalidate all variables that depend on the current code
                    self._last_messages = None
                    self._current_code_fitness = None
                    self._code_mutator = None
                break

            range_bottom = range_top

        #if possible_transitions <= 1:
        #    self._frozen = True

    @property
    def current_code(self):
        return self._current_code

    @property
    def frozen(self):
        return self._frozen
