"""metadata contains data struct of mutation dataset"""
from __future__ import annotations
from itertools import product
from typing import cast, Union, Tuple, List, Dict, Set

from tqdm import tqdm

from cliff.parser.base import Scenery

MultiResidue = Tuple[int]
Seq = Tuple[str]


class NeighbourItem:
    """a `NeighbourItem` means the difference of two Neighboring sequences"""
    # ABCD -> ABCE
    # ABCE
    target: int
    diff: str
    # 3(index from 0)
    index: MultiResidue


class Neighbourhood:
    """`Neighbourhood` means the adjacency list of sequences"""

    def __init__(
        self,
        use_residues: Tuple[MultiResidue],
        sequence: List[str],
        variables: Set[str],
        tqdm_enable: bool,
    ) -> None:
        self.sequence: List[str] = sequence
        self.tqdm_enable = tqdm_enable

        # inferred attributes
        self.sequence_length: int = len(sequence[0])
        self.sequence_num: int = len(self.sequence)
        self.seq_to_index: Dict[str, int] = {
            seq: index for index, seq in enumerate(sequence)
        }
        self.res_variables: Dict[MultiResidue, Tuple[Seq]] = {
            res: tuple(product(variables, repeat=2)) for res in use_residues
        }
        # used for inner calculation
        # [((0, 1), (A, B)), (2, A)]
        self.substitude_tuple: List[Tuple[MultiResidue, Seq]] = []
        for res in use_residues:
            self.substitude_tuple.extend(
                [(res, seq) for seq in self.res_variables[res]]
            )

    @staticmethod
    def substitude(target: str, index: MultiResidue, sub: Seq) -> str:
        """generate a mutation of a sequence"""
        ret = []
        last = 0
        for i, res in enumerate(index):
            ret.append(target[last:res])
            last = len(target) if i == len(index) - 1 else index[i + 1]
            ret.append(sub[i])
            ret.append(target[res + 1: last])
        return "".join(ret)

    @staticmethod
    def select(target: str, index: MultiResidue) -> str:
        """select some part of a str"""
        return "".join([target[i] for i in index])

    def prefetch_neighbour(self) -> Dict[int, Tuple[NeighbourItem]]:
        """calculate and store the neighbour"""
        neighbour: Dict[int, List[NeighbourItem]] = {
            i: [] for i in range(self.sequence_num)
        }
        merge_iter = product(range(self.sequence_num), self.substitude_tuple)
        if self.tqdm_enable:
            merge_iter = tqdm(
                merge_iter,
                desc="creating neighbour",
                total=self.sequence_num * len(self.substitude_tuple),
            )
        for seq_index, (sub_index, sub_char) in merge_iter:
            seq_index, sub_index, sub_char = cast(
                Tuple[int, MultiResidue,
                      Seq], (seq_index, sub_index, sub_char),
            )
            new_seq = self.substitude(
                self.sequence[seq_index], sub_index, sub_char)
            if new_seq not in self.seq_to_index:
                continue
            item = NeighbourItem()
            item.target = self.seq_to_index[new_seq]

            diff_before = self.select(self.sequence[seq_index], sub_index)
            diff_after = self.select(new_seq, sub_index)
            item.diff = f"{diff_before}{diff_after}"
            item.index = sub_index

            neighbour[seq_index].append(item)
        return {key: tuple(value) for key, value in neighbour.items()}

    def get(self) -> Dict[int, Tuple[NeighbourItem]]:
        """store and return the neighbour"""
        neighbour = self.prefetch_neighbour()
        return neighbour


class Dictionary:
    """used variables of sequence"""

    def __init__(self) -> None:
        self.chars: Set[str]

    @classmethod
    def from_factory(cls, src: Union[List[str], str]) -> Dictionary:
        """make a distionary"""
        dic = cls()
        dic.chars = set(src)
        return dic


class MetaData:
    """a model for mutation dataset"""

    def __init__(
        self,
        scenery: Scenery,
        chars: Union[List[str], str],
    ) -> None:

        self.dictionary = Dictionary.from_factory(chars)

        self.variables: Set[str] = set(self.dictionary.chars)

        # set attributes
        self.sequence_length: int = len(scenery.sequence[0])
        self.sequence: List[str] = scenery.sequence

        self.fitness = scenery.fitness

        # inferred attributes
        self.sequence_num: int = len(self.sequence)
        self.seq_index: Dict[str, int] = {
            seq: index for index, seq in enumerate(scenery.sequence)
        }
        assert len(self.sequence) == len(self.seq_index)

        # lazy attributes
        self.neighbour: Dict[int, Tuple[NeighbourItem]] = {}

    def get_neighbour(
        self, use_keys: Tuple[MultiResidue] = tuple(), tqdm_enable=True
    ) -> None:
        """fetch the neighbour adjacency list"""
        use_keys = tuple((i,) for i in range(self.sequence_length))
        self.neighbour: Dict[int, Tuple[NeighbourItem]] = Neighbourhood(
            use_keys, self.sequence, self.variables, tqdm_enable
        ).get()
