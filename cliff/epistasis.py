from itertools import combinations, product, zip_longest
import logging
from math import isnan
import sys
from typing import Union, cast
from bisect import insort, bisect_left

import numpy as np
from prettytable import PrettyTable
from tqdm import tqdm

from cliff.metadata import MetaData, NeighbourItem
from cliff.parser.base import Scenery
from cliff.epi_utils import *


logging.basicConfig(level=logging.INFO)

LARGE = sys.float_info.max

class Epi2Prob:
    def __init__(self, possible_bases : set[MULTI_RESIDUE], epi: dict[MULTI_RESIDUE, VAL_EACH_RESIDUE]):
        self.possible_bases : set[MULTI_RESIDUE] = possible_bases
        self.epi = epi
        self.base_probability: dict[MULTI_RESIDUE, float] = None

    def _calculate_prob(self):
        all_possiable_base_exist: dict[MULTI_RESIDUE, int] = {key: 0 for key in self.possible_bases}
        all_actual_base_exist: dict[MULTI_RESIDUE, int] = {key: 0 for key in self.possible_bases}
        for key, (bases, _) in self.epi.items():
            if key in bases:
                continue
            all_possiable_bases = gen_lower_all_base(key)
            for possible_base in all_possiable_bases:
                all_possiable_base_exist[possible_base] += 1
            for actual_base in bases:
                if actual_base in all_actual_base_exist:
                    all_actual_base_exist[actual_base] += 1
        self.base_probability = {key: all_actual_base_exist[key] / all_possiable_base_exist[key] for key in self.possible_bases}
    
    def table_show(self):
        if self.base_probability ==None:
            self._calculate_prob()

        table = PrettyTable()
        max_len = max(len(key) for key in self.base_probability.keys())
        splits_keys = [[] for _ in range(max_len)]
        splits_vals = [[] for _ in range(max_len)]
        for key, val in self.base_probability.items():
            splits_keys[len(key) - 1].append(key)
            splits_vals[len(key) - 1].append(val)
        
        for i in range(max_len):
            sorted_ind = argsort(splits_keys[i])
            splits_keys[i] = [splits_keys[i][j] for j in sorted_ind]
            splits_vals[i] = [splits_vals[i][j] for j in sorted_ind]

        splits_keys = list(zip(*list(zip_longest(*splits_keys, fillvalue=""))))
        splits_vals = list(zip(*list(zip_longest(*splits_vals, fillvalue=""))))

        for i in range(max_len):
            table.add_column("res_{}".format(i+1), splits_keys[i])
            table.add_column("prob_{}".format(i+1), splits_vals[i])
        print(table)

class Epistasis:
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
        self.possible_bases : set[MULTI_RESIDUE] = set()

    def to_prob(self, epi: dict[MULTI_RESIDUE, VAL_EACH_RESIDUE]) -> Epi2Prob:
        return Epi2Prob(self.possible_bases, epi)

    def gen_lower_base(
        self, sorted_at_key: MULTI_RESIDUE,
    ) -> list[tuple[MULTI_RESIDUE]]:

        ret: list[tuple[MULTI_RESIDUE]] = [(sorted_at_key,)]
        if len(sorted_at_key) == 1:
            return ret
        sorted_at_key = cast(MULTI_RESIDUE, sorted_at_key)
        # 获得基准组合
        ret.append(tuple((i,) for i in sorted_at_key))
        if len(sorted_at_key) == 2:
            return ret
        # 获得低一阶的全部组合
        # (1, 2, 3) = (1) + (2, 3) / (2) + (1, 3) / (3) + (1, 2)
        for i in sorted_at_key:
            tmp = list(sorted_at_key)
            tmp.__delitem__(bisect_left(tmp, i))
            one: tuple[MULTI_RESIDUE] = tuple([tuple(tmp), (i,)])
            ret.append(one)
        return ret

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
            tgt_key = select_substr(self.sequence[neighbor.target], sorted_at_key)
            diff_group.setdefault((src_key, tgt_key), []).append(
                self.fitness[i] - self.fitness[neighbor.target]
            )

        # 清理diff_group/keys
        diff_group = {key: value for key, value in diff_group.items() if len(value) > 0}
        keys_set: set[SEQ] = set()
        for key, _ in diff_group.items():
            src, tgt = key
            if src not in keys_set:
                keys_set.add(src)
            if tgt not in keys_set:
                keys_set.add(tgt)
        # 根据diff_group计算平均差
        possiable_keys = list(keys_set)
        diff = {key: sum(value) / len(value) for key, value in diff_group.items()}

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
                self.calculate_order_base(lower_multi_res, lower_base, neighbour)
            if lower_seq in self.epi_net_base[lower_base][lower_multi_res]:
                epi_values[seq] -= self.epi_net_base[lower_base][lower_multi_res][
                    lower_seq
                ]

            else:
                epi_values[seq] = float("nan")

        self.epi_net_base.setdefault(used_base, {})[sorted_at_key] = epi_values
        var: float = np.var([val for _, val in epi_values.items() if not isnan(val)])
        return var

    def calculate_order(
        self, sorted_at_key: MULTI_RESIDUE, early_stop=True,
    ) -> VAL_EACH_RESIDUE:
        # 计算某个位点组合的所有氨基酸组合的上位
        # 计算（0，1，2）其中计差组（1，2） 变化组（0）A->E 其余组(3, 4)
        # 两组样例数据
        # ABCDE -> EBCDE = ABC - EBC
        # ABCPQ -> EBCPQ = ABC - EBC

        # 计算实时局部差值邻接表
        # (0, 1, 2) = （0， 1）+（2） 实时计算它的局部差值
        # 计算可用基
        bases = self.gen_lower_base(sorted_at_key)
        epi_ans: VAL_EACH_RESIDUE = (tuple(), LARGE)
        for used_base in tqdm(bases):
            var = self.calculate_order_base(sorted_at_key, used_base)
            if var < epi_ans[1]:
                epi_ans = (used_base, var)
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
            if i < self.max_order:
                self.possible_bases.update(set(epi_order_keys))
            # 此处可做并行优化
            for sorted_at_key in epi_order_keys:
                logging.info("epistasis net -- order {}".format(sorted_at_key))
                sorted_at_key = cast(MULTI_RESIDUE, sorted_at_key)
                epi_order = self.calculate_order(sorted_at_key, True)
                epi_val[sorted_at_key] = epi_order

        return epi_val
