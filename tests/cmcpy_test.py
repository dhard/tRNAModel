import numpy as np


class Test_Space():
    def get_aars_aa_map(self, aarss):
        return np.diag([1] * len(aarss))

    def get_trna_aars_map(self, trnas, aarss):
        pass

    def get_codon_trna_map(self, trnas):
        return np.diag([1] * len(trnas))
