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

    def calculate(self) -> float:
        """
        calculate ruggness of a scenery

        Returns
        -------
        ruggness : float
            ruggness of scenery
        """
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
                index = item.index[0]
                diff_value = self.fitness[item.target] - self.fitness[i]
                mutation_derivation["{}{}".format(index, item.diff)].append(
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
