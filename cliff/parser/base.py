import abc
from typing import List


class Scenery:
    sequence: List[str]
    fitness: List[float]


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
