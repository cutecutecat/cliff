from bisect import insort
from itertools import combinations
from typing import Generator, Union, cast

import numpy as np
from sklearn.linear_model import LinearRegression

from cliff.metadata import SEQ, MULTI_RESIDUE

SEQ_DIFF = tuple[SEQ, SEQ]
VAL_EACH_RESIDUE = dict[tuple[MULTI_RESIDUE], float]
EPI_EACH_RESIDUE = dict[SEQ, float]
EPI_NET = dict[MULTI_RESIDUE, EPI_EACH_RESIDUE]
EPI_NET_BASE = dict[tuple[MULTI_RESIDUE], EPI_NET]


def select_substr(
    src: Union[str, SEQ, list[str]], select: Union[MULTI_RESIDUE, list[int]]
) -> SEQ:
    """
    utility for select sequence subset from string

    Parameters
    ----------
    src: str
        target string to be select

    select: list[int]
        list of index to select

    Returns
    -------
    substr : list[str]
        selected sub-string

    Examples
    --------
    >> assert(select_substr("ABCDE", [0, 2, 3]) == ('A', 'C', 'E'))
    """
    if type(select) == set:
        select = sorted(list(select))
    select = cast(list[int], select)
    return tuple(src[s] for s in select)


def mk_combine_subset(
    used_base: tuple[MULTI_RESIDUE],
) -> list[tuple[MULTI_RESIDUE]]:
    ret: list[tuple[MULTI_RESIDUE]] = []
    for order in range(1, len(used_base)):
        once: list[tuple[MULTI_RESIDUE]] = sorted(
            list(combinations(used_base, order))
        )
        ret.extend(once)
    return ret


def mk_combine(n: int, max_order: int) -> list[MULTI_RESIDUE]:
    ret: list[MULTI_RESIDUE] = []
    for order in range(1, max_order + 1):
        once: list[MULTI_RESIDUE] = list(
            combinations([j for j in range(n)], order))
        ret.extend(once)
    return ret


def mk_combine_with_k(
    k: MULTI_RESIDUE, n: int, max_order: int
) -> list[MULTI_RESIDUE]:
    """
    connect the combination of (1, 2, ...,n) from order 1 to max_order, which must contains k.

    mk_combine_with_k(1, 2, 3) --> (1), (1, 0), (1, 2), (0, 1, 2)
    """
    k_set = set(k)
    ret: list[MULTI_RESIDUE] = []
    for order in range(0, max_order):
        other: list[list[int]] = [
            list(res)
            for res in combinations([j for j in range(n) if j not in k_set], order)
        ]
        for res in other:
            for add in k:
                insort(res, add)
        once: list[MULTI_RESIDUE] = [tuple(res) for res in other]
        ret.extend(once)
    return ret


def gen_lower_base(
    sorted_at_key: MULTI_RESIDUE,
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


def get_epi_from_diff(
    diff: dict[SEQ_DIFF, float], possiable_keys: list[SEQ],
) -> EPI_EACH_RESIDUE:
    """
    calaulate averaging epistasis value from epistasis delta
    by Linear Regression

    Parameters
    ----------
    diff: dict[SEQ_DIFF, float]
        epistasis delta of each variance combination

    possiable_keys: list[SEQ]
        list of variance combination

    Returns
    -------
    epi_values : EPI_EACH_RESIDUE
        averaging epistasis value of all variance combination in each residue
    """
    # 研究位点作为点，邻接组作为边，构图
    # 直接做最小二乘法计算
    # D = [0, d1, d2, d3, ...]
    # V = [v0, v1, v2, ...]
    # H = [1, 1, 1, ..., 1]
    # H =  [..., 1, -1, ...]
    diff_keys = list(diff.keys())
    diff_vals = [diff[key] for key in diff_keys]
    # 结果矩阵
    Y = np.hstack([[0], diff_vals])
    # 系数矩阵
    keys_index = {key: i for i, key in enumerate(possiable_keys)}
    edges_positive = [[i+1, keys_index[diff[0]]]
                      for i, diff in enumerate(diff_keys)]
    edges_negetive = [[i+1, keys_index[diff[1]]]
                      for i, diff in enumerate(diff_keys)]

    D = np.zeros((len(diff) + 1, len(possiable_keys)))
    D[0, :] = 1.0
    D[tuple(zip(*edges_positive))] = 1.0
    D[tuple(zip(*edges_negetive))] = -1.0

    LR = LinearRegression()
    # 训练模型
    # sample_weight = [1.0, 0.1, ...]
    sample_weight = 0.1 * np.ones_like(Y)
    sample_weight[0] = 1.0
    LR.fit(D, Y, sample_weight)
    X = LR.coef_
    epi_values = {possiable_keys[ind]: val for ind, val in enumerate(X)}
    return epi_values


def fetch_lower_select(
    lower_multi_res: MULTI_RESIDUE, sorted_at_key: MULTI_RESIDUE
) -> MULTI_RESIDUE:
    index = argsort(sorted_at_key)
    res_to_ind = dict(zip(sorted_at_key, index))
    return tuple(res_to_ind[i] for i in lower_multi_res)


def argsort(seq):
    # http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
    return sorted(range(len(seq)), key=seq.__getitem__)
