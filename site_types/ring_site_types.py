import _Site_Types

class Ring_Site_Types(_Site_Types):
    def __init__(self, phi, amino_acids, sites, weights):
        super(Ring_Site_Types, self).__init__(phi, amino_acids, sites, weights)

    def distance(self, site, amino_acid):
        return min(abs(site - amino_acid), (1 - abs(site - amino_acid))) # Ring space from cmcpy
