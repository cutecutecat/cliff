"""abstract interface for parser"""
import abc
from typing import List


class Scenery:
    """data struct for mutation dataset"""
    sequence: List[str]
    fitness: List[float]


class Parser(metaclass=abc.ABCMeta):
    """
    abstract API for parser of any type
    """
    @classmethod
    @abc.abstractmethod
    def parse(cls) -> Scenery:
        """
        generate scenery from source
        """
        pass
