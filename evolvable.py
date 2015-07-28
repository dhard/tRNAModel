
class Evolvable(object):
    def __init__(self, sequence):
        self._sequence = sequence

    @property
    def sequence(self):
        return self._sequence

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __hash__(self):
        return hash(sequence)
