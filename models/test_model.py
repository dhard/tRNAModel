import numpy as np

from models import Model

from aars_space.id_aars_space import Id_AARS_Space
from trna_space.bit_trna_space import Bit_TRNA_Space
from site_types.ring_site_types import Ring_Site_Types
from code_mutator import Code_Mutator
from code import Code
from evolver import Evolver

def get_ring_mutation_matrix(size, mu):
    mut_matrix = np.diag([1.0 - (2 * mu)] * size)

    for i in xrange(size):
        mut_matrix[i][i - 1] = mu
        mut_matrix[i][i + 1 if (i + 1) < size else 0 ] = mu # Handles wraping around to 0

    return mut_matrix

class Test_Model(Model):
    def __init__(self, args, rng):
        super(Test_Model, self).__init__(args, rng)

    def get_summary(self):
        print "Running test model"

    def get_output(self):
        print self.evolver.current_code
        print "Model done"

    def run(self):
        trnas = 8
        aarss = 8
        phi = 0.9
        mu = .01
        msg_mu = .1
        pop = 100

        trna_names = ["tra_" + str(i) for i in xrange(trnas)]
        codon_names = ["codon_" + str(i) for i in xrange(trnas)]

        aars_names = ["aars_" + str(i) for i in xrange(aarss)]
        aa_names = ["aa_" + str(i) for i in xrange(aarss)]

        site_names = ["site_" + str(i) for i in xrange(aarss)]

        trna_space = Bit_TRNA_Space(trna_names, codon_names)
        aars_space = Id_AARS_Space(aars_names, aa_names)

        trna = range(trnas)
        trna[1] = 0
        trna[2] = 0
        code = Code(trna, range(aarss), trna_space, aars_space)
        
        site_types = Ring_Site_Types(phi, zip(aa_names, np.arange(0, 1, .125)),
                                     zip(site_names, np.arange(0, 1, .125)), [1] * len(site_names))

        
        self.evolver = Evolver(code, site_types, mu,
                          get_ring_mutation_matrix(len(codon_names), msg_mu), pop, self.rng)

        mut = Code_Mutator(code, mu)

        trna[0] = 1
        code2 = Code(trna, range(aarss), trna_space, aars_space)
        
        #print mut.mutation_probability(code2)
        #print [trna_space.mutation_probability(0, i, mu) for i in xrange(8)]
        # mut = Code_Mutator(code, mu)
        # code_dict = mut.get_one_gene_mutation_probabilities()
        #for code in code_dict:
        #    print code, code_dict[code]
        for i in xrange(10000):
            self.evolver.step_time()
