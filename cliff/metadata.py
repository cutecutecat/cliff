from __future__ import annotations
from os import getcwd
from os.path import isfile, join
import pickle
from itertools import product
from typing import cast, Union

from tqdm import tqdm

from cliff.parser.base import Scenery

MULTI_RESIDUE = tuple[int]
SEQ = tuple[str]

class NeighbourItem:
    # ABCD -> ABCE
    # ABCE
    target: int
    diff: str
    # 3(index from 0)
    index: MULTI_RESIDUE


class Neighbourhood:
    def __init__(
        self,
        use_residues: tuple[MULTI_RESIDUE],
        sequence: list[str],
        variables: set[str],
        neighbour_path: str,
        tqdm_enable: bool,
    ) -> None:
        self.sequence: list[str] = sequence
        self.neighbour_path: str = neighbour_path
        self.tqdm_enable = tqdm_enable

        # inferred attributes
        self.sequence_length: int = len(sequence[0])
        self.sequence_num: int = len(self.sequence)
        self.seq_to_index: dict[str, int] = {
            seq: index for index, seq in enumerate(sequence)
        }
        self.res_variables: dict[MULTI_RESIDUE, tuple[SEQ]] = {
            res: tuple(product(variables, repeat=2)) for res in use_residues
        }
        # used for inner calculation
        # [((0, 1), (A, B)), (2, A)]
        self.SUBSTITUDE_TUPLE: list[tuple[MULTI_RESIDUE, SEQ]] = []
        for res in use_residues:
            self.SUBSTITUDE_TUPLE.extend(
                [(res, seq) for seq in self.res_variables[res]]
            )

    @staticmethod
    def substitude(target: str, index: MULTI_RESIDUE, sub: SEQ) -> str:
        ret = []
        last = 0
        for i, res in enumerate(index):
            ret.append(target[last:res])
            last = len(target) if i == len(index) - 1 else index[i + 1]
            ret.append(sub[i])
            ret.append(target[res + 1 : last])
        return "".join(ret)

    @staticmethod
    def select(target: str, index: MULTI_RESIDUE) -> str:
        return "".join([target[i] for i in index])

    def prefetch_neighbour(self) -> dict[int, tuple[NeighbourItem]]:
        neighbour: dict[int, list[NeighbourItem]] = {
            i: [] for i in range(self.sequence_num)
        }
        iter = product(range(self.sequence_num), self.SUBSTITUDE_TUPLE)
        if self.tqdm_enable:
            iter = tqdm(
                iter,
                desc="creating {}".format(self.neighbour_path),
                total=self.sequence_num * len(self.SUBSTITUDE_TUPLE),
            )
        for seq_index, (sub_index, sub_char) in iter:
            seq_index, sub_index, sub_char = cast(
                tuple[int, MULTI_RESIDUE, SEQ], (seq_index, sub_index, sub_char),
            )
            new_seq = self.substitude(self.sequence[seq_index], sub_index, sub_char)
            if new_seq not in self.seq_to_index:
                continue
            item = NeighbourItem()
            item.target = self.seq_to_index[new_seq]
            item.diff = "{}{}".format(
                self.select(self.sequence[seq_index], sub_index),
                self.select(new_seq, sub_index),
            )
            item.index = sub_index

            neighbour[seq_index].append(item)
        return {key: tuple(value) for key, value in neighbour.items()}

    def load_neighbour(self) -> dict[int, tuple[NeighbourItem]]:
        if not isfile(self.neighbour_path):
            raise RuntimeError("path not exist")

        with open(self.neighbour_path, "rb") as f:
            neighbour: dict[int, tuple[NeighbourItem]] = pickle.load(f)

        return neighbour

    def prefetch(self) -> None:
        if isfile(self.neighbour_path):
            raise RuntimeError("path exist")

        neighbour = self.prefetch_neighbour()

        with open(self.neighbour_path, "wb") as f:
            pickle.dump(neighbour, f)

    def get(self) -> dict[int, tuple[NeighbourItem]]:
        try:
            neighbour = self.load_neighbour()
        except Exception as e:
            neighbour = self.prefetch_neighbour()
        return neighbour


class Dictionary:
    def __init__(self) -> None:
        self.chars: set[str]

    @classmethod
    def from_factory(cls, src: Union[list[str], str]) -> Dictionary:
        if type(src) == str:
            src = cast(str, src)
            with open(src, "r") as f:
                dictionary_str = f.read()
            chars = dictionary_str.rstrip("\n").lstrip("\n").split("\n")
            return cls.from_factory(chars)

        src = cast(list[str], src)
        dic = cls()
        dic.chars = set(src)
        return dic


class MetaData:
    def __init__(
        self,
        scenery: Scenery,
        chars: Union[list[str], str],
        neighbour_path: str = "neighbour.pkl",
    ) -> None:

        self.dictionary = Dictionary.from_factory(chars)

        self.variables: set[str] = set(self.dictionary.chars)

        # set attributes
        self.sequence_length: int = len(scenery.sequence[0])
        self.sequence: list[str] = scenery.sequence
        self.neighbour_path: str = join(getcwd(), neighbour_path)

        self.fitness = scenery.fitness

        # inferred attributes
        self.sequence_num: int = len(self.sequence)
        self.seq_index: dict[str, int] = {
            seq: index for index, seq in enumerate(scenery.sequence)
        }
        assert len(self.sequence) == len(self.seq_index)

        # lazy attributes
        self.neighbour: dict[int, tuple[NeighbourItem]] = {}

    def get_neighbour(
        self, use_keys: tuple[MULTI_RESIDUE] = tuple(), save=False, tqdm_enable=True
    ) -> None:
        use_keys = tuple((i,) for i in range(self.sequence_length))
        if save:
            Neighbourhood(
                use_keys,
                self.sequence,
                self.variables,
                self.neighbour_path,
                tqdm_enable,
            ).prefetch()
        self.neighbour: dict[int, tuple[NeighbourItem]] = Neighbourhood(
            use_keys, self.sequence, self.variables, self.neighbour_path, tqdm_enable
        ).get()

    def prefetch(
        self, use_keys: tuple[MULTI_RESIDUE] = tuple(), tqdm_enable=True
    ) -> None:
        use_keys = tuple((i,) for i in range(self.sequence_length))
        Neighbourhood(
            use_keys, self.sequence, self.variables, self.neighbour_path, tqdm_enable
        ).prefetch()
