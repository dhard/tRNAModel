from abc import ABCMeta, abstractmethod
from .messages import Messages
from .code import Bit_Code
from random import Random


class Evolver(object):
    def __init__(self, inital_code, site_types, message_mutation_matrix, code_mutator, pop_size, rng):
        self._current_code = inital_code
        self._message_mutation_matrix
        self._site_types = site_types
        self._code_mutator = code_mutator
        self._pop_size = pop_size
        self._rng = rng

        message = Messages(initial_code, site_types, message_mutation_matrix)
        self._codon_usage = message.get_equilibrium_codon_usage()

    def get_code_fitness(self, code, codon_usage = None):
        fitness_contributions = None
        messages = Messages(code, self._site_types, self._message_mutation_matrix)
        
        if code_usage == None:
            fitness_contributions = messages.get_equilibrium_fitness_contributions()
        else:
            fitness_contributions = messages.fitness_contributions_at_codon_usage(codon_usage)

        site_weights = self._site_types.site_weights
        return np.prod([fitness_contributions[site]**site_weights[site] for site in xrange(len(site_weights))])
    
    def get_transition_probability(self, start, end, mutation_probability):
        if start == end: # This case must be handled when the transition matrix is already built
            return 0 # It is really 1 - summation(transition_probability(start, end), for all potential values of end)

        probability = log(self._pop_size) + log(self.get_fixation_probability(start, end)) + \ 
                      log(mutation_probability)
        probability = exp(probability)

        assert 0 <= probability <= 1, "Error calculating transition probability, not in 0 to 1 range."

        return probability

    def get_fixation_probability(self, start, end):
        start_codon_usage = None if start == self._current_code else self._codon_usage
        end_codon_usage = None if end == self._current_code else self._codon_usage

        start_fitness = self.get_code_fitness(start, start_codon_usage)
        end_fitness = self.get_code_fitness(end, end_codon_usage)

        if start_fitness == end_fitess:
            return 1.0 / pop_size # This is the limit of the function below

        # Calculates the Moran fixation probability
        num = 1 - (start_fitness / end_fitess) # TODO these calculations in log space
        dem = 1 - (start_fitness / end_fitess)**pop_size

        return num / dem
        
    def step_time(self):
        code_mutants = self._code_mutator.get_mutation_probabilities()
        transition_probabilities = {}
        
        for genetic_code in code_mutants:
            end_code = self._code_mutator.build_code(genetic_code)
            transition_probabilities[end_code] = self.get_transition_probability(self._current_code, end_code,
                                                                                 code_mutants[genetic_code])

        range_bottom = 0
        random_num = self._rng.rand()
        code_changed = False

        assert np.close(sum(transition_probabilities.values()), 1), \
            "Transitions probabilities do not sum to 1."

        for code in transition_probabilities:
            range_top = transition_probabilities[code] + range_bottom

            if range_bottom < random_num <= range_top:
                code_changed = True
                if code != self._current_code:
                    self._current_code = code
                    message = Messages(code, self._site_types, self._message_mutation_matrix)
                    self._codon_usage = message.get_equilibrium_codon_usage()

                break

        assert code_changed, "You used bad logic here double check it!"

            
            
