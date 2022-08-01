"""parser for `mutation` dataset"""
from typing import Union, cast


import pandas as pd

from .base import Parser, Scenery


class SeqArgs:
    """arguments for `sequence` parser"""
    sequence_label: str
    fitness_label: str


class SeqParser(Parser):
    """parser for `mutation` dataset"""

    @classmethod
    def parse(cls, data: Union[str, pd.DataFrame], args: SeqArgs) -> Scenery:
        if isinstance(data, str):
            file = pd.read_csv(data)
            return cls.parse(file, args)

        data = cast(pd.DataFrame, data)
        assert (
            args.sequence_label in data.columns and args.fitness_label in data.columns
        )

        sce = Scenery()
        sce.sequence = data[args.sequence_label].to_list()
        sce.fitness = data[args.fitness_label].to_list()
        return sce
