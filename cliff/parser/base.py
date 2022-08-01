"""abstract interface for parser"""
import abc
from typing import Any, List, Union

import pandas as pd


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
    def parse(cls, data: Union[str, pd.DataFrame], args: Any) -> Scenery:
        """
        generate scenery from source
        """
