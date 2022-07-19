from functools import cmp_to_key
from itertools import combinations, product, zip_longest
import logging
from math import isnan
import sys
from typing import Union, cast
from bisect import insort, bisect_left
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.figure import Figure

import numpy as np
from scipy.special import softmax
from tqdm import tqdm

from cliff.metadata import MetaData, NeighbourItem
from cliff.parser.base import Scenery
from cliff.epi_utils import *


logging.basicConfig(level=logging.INFO)

LARGE = sys.float_info.max


class Epi2Prob:
    @staticmethod
    def sort_key(a: tuple[int], b: tuple[int]):
        first_cmp = (len(a) > len(b)) - (len(a) < len(b))
        if first_cmp != 0:
            return first_cmp
        else:
            return (a > b) - (a < b)

    def __init__(self, possible_bases: set[MULTI_RESIDUE], epi: dict[MULTI_RESIDUE, VAL_EACH_RESIDUE]):
        self.max_keys_num = max(max(b) for b in possible_bases) + 1

        self.possible_bases: list[MULTI_RESIDUE] = sorted(
            list(possible_bases), key=cmp_to_key(self.sort_key))
        self.possible_bases = list(
            filter(lambda b: len(b) > 1, self.possible_bases))
        self.epi = epi
        self.orders = [len(base) for base in self.possible_bases]
        self.base_prob: dict[MULTI_RESIDUE, float] = None

        # Prepare an cycle of colors
        max_order = max(self.orders)
        prop_cycle = plt.rcParams['axes.prop_cycle']
        self.color_cycle = prop_cycle.by_key()['color']
        color_scalar = int(max_order / len(self.color_cycle)) + 1
        self.color_cycle *= color_scalar

    def _calculate_prob(self):
        base_all_prob: dict[MULTI_RESIDUE, list[float]] = {
            key: [] for key in self.possible_bases}
        # softmax化
        for _, val_res in self.epi.items():
            val_res_keys = val_res.keys()

            val_res_softmax = softmax([1/val_res[k] for k in val_res_keys])
            for i, bases in enumerate(val_res_keys):
                for b in bases:
                    if b in base_all_prob:
                        base_all_prob[b].append(val_res_softmax[i])
        self.base_prob = {key: np.mean(
            base_all_prob[key]) for key in self.possible_bases}

    def plot(self) -> Figure:
        if self.base_prob == None:
            self._calculate_prob()

        # Create a plot with an upper and lower panel, sharing the x-axis
        fig = plt.figure()
        gs = mpl.gridspec.GridSpec(2, 1,
                                   height_ratios=[1, 1],
                                   hspace=0.00)
        values = [self.base_prob[k] for k in self.possible_bases]
        ax = [plt.subplot(gs[0])]
        ax.append(plt.subplot(gs[1], sharex=ax[0]))
        bar_axis = ax[0]
        grid_axis = ax[1]

        colors_for_bar = np.array([mpl.colors.colorConverter.to_rgba(
            self.color_cycle[i - 1]) for i in self.orders])
        bar_axis.bar(np.arange(len(values)), values,
                     width=0.9,
                     color=colors_for_bar,
                     linewidth=1)

        grid_corr = np.zeros((self.max_keys_num, len(self.orders), 4))
        for i, bases in enumerate(self.possible_bases):
            for b in bases:
                grid_corr[b, i] = colors_for_bar[i]
        grid_axis.imshow(grid_corr)

        grid_axis.set_xticks(np.arange(-.5, len(self.orders), 1), minor=True)
        grid_axis.set_yticks(np.arange(-.5, self.max_keys_num, 1), minor=True)
        for tick in grid_axis.xaxis.get_major_ticks():
            tick.tick1line.set_visible(False)
            tick.tick2line.set_visible(False)
            tick.label1.set_visible(False)
            tick.label2.set_visible(False)
        for tick in grid_axis.yaxis.get_major_ticks():
            tick.tick1line.set_visible(False)
            tick.tick2line.set_visible(False)
            tick.label1.set_visible(False)
            tick.label2.set_visible(False)
        for tick in bar_axis.yaxis.get_major_ticks():
            tick.tick1line.set_visible(False)
            tick.tick2line.set_visible(False)
            tick.label1.set_visible(False)
            tick.label2.set_visible(False)
        for tick in bar_axis.xaxis.get_major_ticks():
            tick.tick1line.set_visible(False)
            tick.tick2line.set_visible(False)
            tick.label1.set_visible(False)
            tick.label2.set_visible(False)
        grid_axis.grid(which='minor', color='black',
                       linestyle='-', linewidth=1)
        return fig


class Connection:
    def __init__(
        self, scenery: Scenery, max_order: int, variables: Union[list[str], str]
    ) -> None:

        self.variables = variables

        self.scenery = scenery
        self.sequence_length = len(scenery.sequence[0])
        self.fitness = scenery.fitness
        self.sequence = scenery.sequence

        assert 1 <= max_order <= self.sequence_length
        self.max_order = max_order

        # inner calculator varibles
        self.epi_net_base: EPI_NET_BASE = {}
        self.possible_bases: set[MULTI_RESIDUE] = set()

    def to_prob(self, epi: dict[MULTI_RESIDUE, VAL_EACH_RESIDUE]) -> Epi2Prob:
        return Epi2Prob(self.possible_bases, epi)

    def cal_epi_link(
        self, neighbour
    ) -> dict[MULTI_RESIDUE, list[tuple[int, NeighbourItem]]]:
        epi_link: dict[MULTI_RESIDUE, list[tuple[int, NeighbourItem]]] = {
            key: [] for key in mk_combine(self.sequence_length, self.max_order)
        }
        # 此处可做并行优化
        for seq_i_index, nearbys in neighbour.items():
            for neighbour in nearbys:
                epi_keys = mk_combine_with_k(
                    neighbour.index, self.sequence_length, self.max_order
                )
                for key in epi_keys:
                    epi_link[key].append((seq_i_index, neighbour))
        return epi_link

    def calculate_order_base(
        self,
        sorted_at_key: MULTI_RESIDUE,
        used_base: tuple[MULTI_RESIDUE],
        neighbour=None,
    ) -> float:
        used_base = tuple(sorted(used_base))
        if neighbour == None:
            meta = MetaData(self.scenery, self.variables)
            meta.get_neighbour(used_base, tqdm_enable=False)
            neighbour = meta.neighbour
        epi_link = self.cal_epi_link(neighbour)
        diff_group: dict[SEQ_DIFF, list[float]] = {}
        used_neighbours = epi_link[sorted_at_key]

        for i, neighbor in used_neighbours:
            src_key = select_substr(self.sequence[i], sorted_at_key)
            tgt_key = select_substr(
                self.sequence[neighbor.target], sorted_at_key)
            diff_group.setdefault((src_key, tgt_key), []).append(
                self.fitness[i] - self.fitness[neighbor.target]
            )

        # 清理diff_group/keys
        diff_group = {key: value for key,
                      value in diff_group.items() if len(value) > 0}
        keys_set: set[SEQ] = set()
        for key, _ in diff_group.items():
            src, tgt = key
            if src not in keys_set:
                keys_set.add(src)
            if tgt not in keys_set:
                keys_set.add(tgt)
        # 根据diff_group计算平均差
        possiable_keys = list(keys_set)
        diff = {key: sum(value) / len(value)
                for key, value in diff_group.items()}

        epi_values = get_epi_from_diff(diff, possiable_keys)
        # 此处存在差异：SUB操作，计算减去低阶量
        # (0, 1, 2) - (0, 1) - (0, 2) - (1, 2) - (0) - (1) - (2)
        lower_base_comb = mk_combine_subset(used_base)
        for lower_base, seq in product(lower_base_comb, possiable_keys):
            # 有可能低阶数据不存在，这时给一个nan指标
            lower_multi_res_raw: list[int] = []
            for multi_residue in lower_base:
                lower_multi_res_raw.extend([i for i in multi_residue])
            lower_multi_res: MULTI_RESIDUE = tuple(lower_multi_res_raw)

            lower_index = fetch_lower_select(lower_multi_res, sorted_at_key)
            lower_seq = select_substr(seq, lower_index)
            if (
                lower_base not in self.epi_net_base
                or lower_multi_res not in self.epi_net_base[lower_base]
            ):
                self.calculate_order_base(
                    lower_multi_res, lower_base, neighbour)
            if lower_seq in self.epi_net_base[lower_base][lower_multi_res]:
                epi_values[seq] -= self.epi_net_base[lower_base][lower_multi_res][
                    lower_seq
                ]
            else:
                epi_values[seq] = float("nan")

        self.epi_net_base.setdefault(used_base, {})[sorted_at_key] = epi_values
        var: float = np.var(
            [val for _, val in epi_values.items() if not isnan(val)])
        print(var)
        return var

    def calculate_order(self, sorted_at_key: MULTI_RESIDUE) -> VAL_EACH_RESIDUE:
        # 计算某个位点组合的所有氨基酸组合的上位
        # 计算（0，1，2）其中计差组（1，2） 变化组（0）A->E 其余组(3, 4)
        # 两组样例数据
        # ABCDE -> EBCDE = ABC - EBC
        # ABCPQ -> EBCPQ = ABC - EBC

        # 计算实时局部差值邻接表
        # (0, 1, 2) = （0， 1）+（2） 实时计算它的局部差值
        # 计算可用基
        bases = gen_lower_base(sorted_at_key)
        epi_ans: VAL_EACH_RESIDUE = {}
        for used_base in tqdm(bases):
            epi_ans[used_base] = self.calculate_order_base(
                sorted_at_key, used_base)
        return epi_ans

    def calculate(self) -> dict[MULTI_RESIDUE, VAL_EACH_RESIDUE]:
        """
        ret: {(0, 1):{("A","B"): 0.1, ("A","C"): 1.2, ("B","C"): 0.5}}
        """
        epi_val: dict[MULTI_RESIDUE, VAL_EACH_RESIDUE] = {}
        # 计算上位效应
        for i in range(1, self.max_order + 1):
            logging.info("begin calculate epistasis net -- order {}".format(i))
            epi_order_keys: list[MULTI_RESIDUE] = list(
                combinations(range(self.sequence_length), i)
            )
            if i < self.max_order + 1:
                self.possible_bases.update(set(epi_order_keys))
            # 此处可做并行优化
            for sorted_at_key in epi_order_keys:
                logging.info("epistasis net -- order {}".format(sorted_at_key))
                sorted_at_key = cast(MULTI_RESIDUE, sorted_at_key)
                epi_order = self.calculate_order(sorted_at_key)
                epi_val[sorted_at_key] = epi_order
        return epi_val
