from abc import ABCMeta, abstractmethod

class _Model(object):
    __metaclass__ = ABCMeta

    def __init__(self, args, rng):
        """ This function takes a namespace of arguments and a
            random number generator. """
        pass

    @abstractmethod
    def get_model_summary():
        """ This function should return a summary of the model
            describing information like what parameters used. """
        pass

    @abstractmethod
    def get_model_output():
        """ This function should return information about the 
            model after it has finished running. """
        pass

    @abstractmethod
    def run_model():
        """ This function should run the model. """
        pass
