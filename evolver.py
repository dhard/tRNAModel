import numpy as np
from messages import Messages
from code_mutator import Code_Mutator

class Evolver(object):
    def __init__(self, initial_code, site_types, message_mutation_matrix, pop_size, rng):
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

    def _code_fitness(self, code, codon_usage = None):
        # codon usage being None means use the codon usages
        # from the current code
        if codon_usage is None:
            if self._last_messages is None:
                # it needs to be calculated
                self._code_fitness(self._current_code)

            codon_usage = self._last_messages.codon_usage

        msgs = Messages(code.effective_code_matrix, self._fitness_matrix,
                        self._message_mutation_matrix)

        if codon_usage is None and code == self._current_code:
            # Only calculate codon usage at equilibrium for the current code
            msgs.calculate_at_equilibrium()
            self._last_messages = msgs
        else:
            msgs.calculate_at_codon_usage(codon_usage)

        fitness_contributions = msgs.fitness_contributions

        overall_fitness = np.prod([fitness_contributions[site]**self._site_weights[site]
                                   for site in xrange(len(self._site_weights))])

        if codon_usage is None and code == self._current_code:
            self._current_code_fitness = overall_fitness

        return overall_fitness

    def _transition_probability(self, from_, to, mutation_probability):
        if from_ == to: # This case must be handled when the transition matrix is already built
            # It is really 1 - summation(transition_probability(from_, to), 
            #                            for all potential values of end)
            return 0

        probability = log(self._pop_size) + log(self._fixation_probability(from_, to)) + \
                      log(mutation_probability)
        probability = exp(probability)

        assert 0 <= probability <= 1, "Transition probability was not in the range [0,1]."

        return probability

    def _fixation_probability(self, to):
        if self._last_messages is None:
            # sets self._last_messages
            # and self._current_code_fitness
            self._code_fitness(self._current_code)

        to_fitness = self.get_code_fitness(to, self._last_messages.codon_usage)

        if np.isclose(self._current_code_fitness, to_fitess):
            return 1.0 / pop_size # This is the limit of the function below

        # Calculates the Moran fixation probability
        num = 1 - (self._current_code_fitness / to_fitess) # TODO these calculations in log space
        dem = 1 - (self._current_code_fitness / to_fitess)**pop_size

        return num / dem

    def step_time(self):
        if self._code_mutator is None:
            self._code_mutator = Code_Mutator(self._current_code,
                                              self._current_code._aars_space.mutation_matrix,
                                              self._current_code._trna_space.mutation_matrix)

        mutant_prob = self._code_mutator.get_mutation_probabilities()
        
        transition_probabilities = {to_code:self._transition_probability(to_code, mutant_prob[to_code])
                                    for to_code in mutant_prob}

        # Handles the corner case in _transition_probability()
        transition_probabilities[self._current_code] = 1.0 - sum(transition_probabilities.values())
        
        range_bottom = 0
        random_num = self._rng.rand()
        code_changed = False

        for to_code in transition_probabilities:
            range_top = transition_probabilities[to_code] + range_bottom

            if range_bottom <= random_num <= range_top:
                code_changed = True
                if to_code != self._current_code:
                    self._current_code = to_code

                    # Invalidate all variables that depend on the current code
                    self._last_messages = None
                    self._current_code_fitness = None
                    self._code_mutator = None
                break

        assert code_changed, "You used bad logic here double check it!"
