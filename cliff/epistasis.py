"""Cauculation of dataset Epistasis"""
from functools import cmp_to_key
from itertools import combinations, product, zip_longest
import logging
import sys
from typing import Union, Set, Tuple, Dict, List

from matplotlib import colors, gridspec
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
import numpy as np
from joblib import Parallel, delayed

from cliff.metadata import MetaData, NeighbourItem
from cliff.parser.base import Scenery
from cliff.epi_utils import (select_substr,
                             mk_combine_subset,
                             mk_combine,
                             mk_combine_with_k,
                             get_epi_from_diff,
                             fetch_lower_select,
                             MultiResidue,
                             EpiResidue,
                             EpiNet,
                             SeqDiff,
                             Seq)


logging.basicConfig(level=logging.INFO)

LARGE = sys.float_info.max


class Epi2Show:
    """module for plot Epistasis"""

    @staticmethod
    def sort_key(one: Tuple[int], two: Tuple[int]):
        """sorted two tuple by length then elements"""
        first_cmp = (len(one) > len(two)) - (len(one) < len(two))
        if first_cmp != 0:
            return first_cmp
        return (one > two) - (one < two)

    def __init__(self, varible: Tuple[str], possible_keys: Set[MultiResidue],
                 epi: Dict[MultiResidue, EpiResidue]):
        self.max_keys_num = max(max(b) for b in possible_keys) + 1
        self.varibles = varible
        self.varibles_index = {key: i for i, key in enumerate(varible)}

        self.possible_keys: List[MultiResidue] = sorted(
            list(possible_keys), key=cmp_to_key(self.sort_key))
        self.epi = epi
        self.orders = [len(base) for base in self.possible_keys]
        self.cols = sum(len(epi[base]) for base in self.possible_keys)

        # Prepare an cycle of colors
        max_order = max(self.orders)
        prop_cycle = plt.rcParams['axes.prop_cycle']
        self.color_cycle = prop_cycle.by_key()['color']
        color_scalar = int(max_order / len(self.color_cycle)) + 1
        self.color_cycle *= color_scalar

        select = []
        for base in self.possible_keys:
            select.extend([len(base)] * len(epi[base]))
        self.colors_for_graph = np.array([colors.colorConverter.to_rgba(
            self.color_cycle[i - 1]) for i in select])
        # Prepare figure
        self.fig = plt.figure()
        grid_spec = gridspec.GridSpec(3, 1,
                                      height_ratios=[1, 1, 0.3],
                                      hspace=0.00)

        all_axis = [plt.subplot(grid_spec[0])]
        all_axis.append(plt.subplot(grid_spec[1], sharex=all_axis[0]))
        all_axis.append(plt.subplot(grid_spec[2], sharex=all_axis[0]))
        self.bar_axis = all_axis[0]
        self.residue_axis = all_axis[1]
        self.chars_table = all_axis[2]
        # Set tick invisible
        for axis in [self.bar_axis.xaxis,
                     self.residue_axis.xaxis, self.residue_axis.yaxis,
                     self.chars_table.xaxis, self.chars_table.yaxis]:
            self.invisible(axis)

    @staticmethod
    def invisible(axis):
        """make axis invisible"""
        for tick in axis.get_major_ticks():
            tick.tick1line.set_visible(False)
            tick.tick2line.set_visible(False)
            tick.label1.set_visible(False)
            tick.label2.set_visible(False)

    def plot(self) -> Figure:
        """draw a graph of Epistasis"""
        values = []
        for bases in self.possible_keys:
            inner_dict = self.epi[bases]
            chars = sorted(inner_dict.keys())
            values.extend([inner_dict[c] for c in chars])

        self.bar_axis.bar(np.arange(len(values)), values,
                          width=0.9,
                          color=self.colors_for_graph,
                          linewidth=1)

        now_at = 0
        residue_corr = np.zeros((self.max_keys_num, self.cols, 4))

        chars_table = []
        for bases in self.possible_keys:
            for chars, _ in self.epi[bases].items():
                residue_corr[bases, now_at, :] = self.colors_for_graph[now_at]

                chars_table.append(list(chars))
                now_at += 1
        chars_col = list(zip_longest(*chars_table, fillvalue=''))
        self.chars_table.table(cellText=chars_col)
        self.chars_table.axis('off')
        self.residue_axis.imshow(residue_corr)

        self.residue_axis.set_xticks(np.arange(-.5, self.cols, 1), minor=True)
        self.residue_axis.set_yticks(
            np.arange(-.5, self.max_keys_num, 1), minor=True)
        self.residue_axis.grid(which='minor', color='black',
                               linestyle='-', linewidth=1)
        return self.fig


class Epistasis:
    """Cauculation of dataset Ruggness"""

    def __init__(
        self, scenery: Scenery, max_order: int, variables: Union[List[str], str]
    ) -> None:

        self.variables = variables

        self.scenery = scenery
        self.sequence_length = len(scenery.sequence[0])
        self.fitness = scenery.fitness
        self.sequence = scenery.sequence

        assert 1 <= max_order <= self.sequence_length
        self.max_order = max_order

        # inner calculator varibles
        self.epi_net: EpiNet = {}
        self.possible_keys: Set[MultiResidue] = set()

    def to_draw(self, epi: Dict[MultiResidue, EpiResidue]) -> Epi2Show:
        """plot Epistasis"""
        return Epi2Show(self.variables, self.possible_keys, epi)

    def cal_epi_link(
        self, target: NeighbourItem
    ) -> Dict[MultiResidue, List[Tuple[int, NeighbourItem]]]:
        """calculate neighbour of every residue combinations"""
        epi_link: Dict[MultiResidue, List[Tuple[int, NeighbourItem]]] = {
            key: [] for key in mk_combine(self.sequence_length, self.max_order)
        }
        for seq_i_index, nearbys in target.items():
            for neighbour in nearbys:
                epi_keys = mk_combine_with_k(
                    neighbour.index, self.sequence_length, self.max_order
                )
                for key in epi_keys:
                    epi_link[key].append((seq_i_index, neighbour))
        return epi_link

    def cal_order(
        self,
        sorted_at_key: MultiResidue,
        neighbour=None,
    ) -> EpiResidue:
        """calculate Epistasis of a residue combinations"""
        used_base = sorted_at_key
        if neighbour is None:
            meta = MetaData(self.scenery, self.variables)
            meta.get_neighbour(used_base, tqdm_enable=False)
            neighbour = meta.neighbour
        epi_link = self.cal_epi_link(neighbour)
        diff_group: Dict[SeqDiff, List[float]] = {}
        used_neighbours = epi_link[sorted_at_key]

        for i, neighbor in used_neighbours:
            src_key = select_substr(self.sequence[i], sorted_at_key)
            tgt_key = select_substr(
                self.sequence[neighbor.target], sorted_at_key)
            diff_group.setdefault((src_key, tgt_key), []).append(
                self.fitness[i] - self.fitness[neighbor.target]
            )

        diff_group = {key: value for key,
                      value in diff_group.items() if len(value) > 0}
        keys_set: Set[Seq] = set()
        for key, _ in diff_group.items():
            src, tgt = key
            if src not in keys_set:
                keys_set.add(src)
            if tgt not in keys_set:
                keys_set.add(tgt)

        possiable_keys = list(keys_set)
        diff = {key: sum(value) / len(value)
                for key, value in diff_group.items()}

        epi_values = get_epi_from_diff(diff, possiable_keys)

        return possiable_keys, epi_values

    def sub(self, epi_value: EpiResidue, possiable_keys: List[Seq],
            sorted_at_key: MultiResidue):
        """substitude lower-order contribution"""
        lower_base_comb = mk_combine_subset(sorted_at_key)
        for lower_base, seq in product(lower_base_comb, possiable_keys):
            lower_index = fetch_lower_select(lower_base, sorted_at_key)
            lower_seq = select_substr(seq, lower_index)
            epi_value[seq] -= self.epi_net[lower_base][lower_seq]
        self.epi_net[sorted_at_key] = epi_value

    def calculate(self) -> Dict[MultiResidue, EpiResidue]:
        """
        calculate epistasis of a scenery

        Returns
        -------
        epistasis : Dict[MultiResidue, EpiResidue]
            epistasis of scenery
        """
        epi_order_keys: List[MultiResidue] = []
        for i in range(1, self.max_order + 1):
            epi_order_keys.extend(
                list(combinations(range(self.sequence_length), i)))
            if i < self.max_order + 1:
                self.possible_keys.update(set(epi_order_keys))
        all_ans = Parallel(n_jobs=len(epi_order_keys))(delayed(self.cal_order)(
            sorted_at_key) for sorted_at_key in epi_order_keys)

        for (possiable_keys, epi_value), sorted_at_key in zip(all_ans, epi_order_keys):
            self.sub(epi_value, possiable_keys, sorted_at_key)
        return self.epi_net
