from abc import ABCMeta, abstractmethod

class Model(object):
    __metaclass__ = ABCMeta

    def __init__(self, args, rng):
        """ This function takes a namespace of arguments and a
            random number generator. """
        self.args = args
        self.rng = rng

    @abstractmethod
    def get_summary(self):
        """ This function should return a summary of the model
            describing information like what parameters used. """
        pass

    @abstractmethod
    def get_output(self):
        """ This function should return information about the 
            model after it has finished running. """
        pass

    @abstractmethod
    def run(self):
        """ This function should run the model. """
        pass
