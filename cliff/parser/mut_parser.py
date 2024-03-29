"""parser for `mutation` dataset"""
from typing import Union, cast

import pandas as pd

from .base import Parser, Scenery


class MutArgs:
    """arguments for `mutation` parser"""
    mutation_label: str
    fitness_label: str
    wile_type: str
    # 0 for index range [1 -> num], 1 for index range [0 -> num-1]
    vt_offset: int


class MutParser(Parser):
    """parser for `mutation` dataset"""

    @staticmethod
    def generate_mut_seq(mut_line: str, wild_type: str, vt_offset: int) -> str:
        """generate mutation sequence of a sequence"""
        if len(mut_line) == 0:
            return wild_type
        muts = mut_line.split(":")
        now = wild_type
        for mut in muts:
            index = int(mut[1:-1]) + vt_offset - 1
            assert (
                now[index] == mut[0]
            ), f"mismatch between mutation {mut[0]}{mut[1:-1]} and wild-type {now[index]}{index}"
            now = now[:index] + mut[-1] + now[index + 1:]
        return now

    @classmethod
    def parse(cls, data: Union[str, pd.DataFrame], args: MutArgs) -> Scenery:
        if isinstance(data, str):
            file = pd.read_csv(
                data, dtype={args.mutation_label: str, args.fitness_label: float})
            return cls.parse(file, args)

        data = cast(pd.DataFrame, data)
        assert (
            args.mutation_label in data.columns and args.fitness_label in data.columns
        )

        sce = Scenery()
        data[args.mutation_label].fillna('', inplace=True)
        sce.sequence = (
            data[args.mutation_label]
            .apply(cls.generate_mut_seq, args=(args.wile_type, args.vt_offset))
            .to_list()
        )
        sce.fitness = data[args.fitness_label].to_list()
        return sce
