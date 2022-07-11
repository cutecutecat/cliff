from itertools import product

import numpy as np

from cliff.metadata import MetaData

class Ruggness:
    def __init__(self, meta: MetaData) -> None:
        self.meta = meta
        if len(meta.neighbour) == 0:
            self.meta.get_neighbour(save=True)

        self.sequence_length = self.meta.sequence_length
        self.neighbour = self.meta.neighbour
        self.variables = self.meta.variables
        self.fitness = self.meta.fitness

    @staticmethod
    def dif_char(A: str, B: str):
        for i, (c1, c2) in enumerate(zip(A, B)):
            if c1 != c2:
                return "{:d}{:s}{:s}".format(i, c1, c2)
        raise ValueError

    def calculate(self) -> float:
        mutation_label = [
            "{}{}{}".format(m, n, k)
            for m, n, k in product(
                range(self.sequence_length), self.variables, self.variables
            )
            if n != k
        ]
        mutation_derivation: dict[str, list[float]] = dict(
            zip(mutation_label, [[] for _ in range(len(mutation_label))])
        )

        for i in range(self.sequence_length):
            neighbour_of_one = self.neighbour[i]
            items = filter(lambda x: x.target > i, neighbour_of_one)

            for item in items:
                diff_value = self.fitness[item.target] - self.fitness[i]
                mutation_derivation["{}{}".format(item.index, item.diff)].append(
                    diff_value
                )

        mutation_derivation = {
            key: value for key, value in mutation_derivation.items() if len(value) > 0
        }
        diff_recenter = []
        for key in mutation_derivation.keys():
            mean = np.mean(mutation_derivation[key])
            diff_recenter.extend([v - mean for v in mutation_derivation[key]])
        return np.var(diff_recenter)
