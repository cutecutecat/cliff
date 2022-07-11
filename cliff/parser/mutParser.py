from typing import Union, cast

from .base import Parser, Scenery
import pandas as pd


class MutArgs:
    mutation_label: str
    fitness_label: str
    wile_type: str
    # 0 for index range [1 -> num], 1 for index range [0 -> num-1]
    vt_offset: int


class MutParser(Parser):
    @staticmethod
    def generate_mut_seq(mut_line: str, wt: str, vt_offset: int) -> str:
        if len(mut_line) == 0:
            return wt
        muts = mut_line.split(":")
        now = wt
        for mut in muts:
            index = int(mut[1:-1]) + vt_offset - 1
            assert (
                now[index] == mut[0]
            ), f"mismatch between mutation {mut[0]}{mut[1:-1]} and wild-type {now[index]}{index}"
            now = now[:index] + mut[-1] + now[index + 1 :]
        return now

    @classmethod
    def parse(cls, data: Union[str, pd.DataFrame], args: MutArgs) -> Scenery:
        if type(data) == str:
            file = pd.read_csv(data)
            return cls.parse(file, args)

        data = cast(pd.DataFrame, data)
        assert (
            args.mutation_label in data.columns and args.fitness_label in data.columns
        )

        sce = Scenery()
        sce.sequence = (
            data[args.mutation_label]
            .apply(cls.generate_mut_seq, args=(args.wile_type, args.vt_offset))
            .to_list()
        )
        sce.fitness = data[args.fitness_label].to_list()
        return sce
