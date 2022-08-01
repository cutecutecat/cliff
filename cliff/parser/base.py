import abc


class Scenery:
    sequence: list[str]
    fitness: list[float]


class Parser(metaclass=abc.ABCMeta):
    """
    Abstract API for parser of any type
    """
    @classmethod
    @abc.abstractmethod
    def parse(cls) -> Scenery:
        """
        generate scenery from source
        """
        pass
