"""utils for cauculating epistasis"""
from bisect import insort
from itertools import combinations
from typing import Generator, Union, cast, List, Tuple, Dict, Set

import networkx as nx

from cliff.metadata import Seq, MultiResidue

SeqDiff = Tuple[Seq, Seq]
EpiResidue = Dict[Seq, float]
EpiNet = Dict[MultiResidue, EpiResidue]


def select_substr(
    src: Union[str, Seq, List[str]], select: Union[MultiResidue, List[int]]
) -> Seq:
    """
    select sequence subset from string

    Parameters
    ----------
    src: str
        target string to be select

    select: List[int]
        list of index to select

    Returns
    -------
    substr : List[str]
        selected sub-string

    Examples
    --------
    >> assert(select_substr("ABCDE", [0, 2, 3]) == ('A', 'C', 'E'))
    """
    if isinstance(select, set):
        select = sorted(List(select))
    select = cast(List[int], select)
    return tuple(src[s] for s in select)


def mk_combine_subset(
    all_res: Tuple[MultiResidue],
) -> List[Tuple[MultiResidue]]:
    """
    generate combination of multiply residues

    Parameters
    ----------
    all_res: Tuple[MultiResidue]
        several residues of comb input

    Returns
    -------
    comb : List[Tuple[MultiResidue]]
        result of combinations

    Examples
    --------
    >> assert(mk_combine_subset(((1, 2), (3, 4))) == (((1, 2),), ((3, 4),), ((1, 2), (3, 4))))
    """
    ret: List[Tuple[MultiResidue]] = []
    for order in range(1, len(all_res)):
        once: List[Tuple[MultiResidue]] = sorted(
            list(combinations(all_res, order))
        )
        ret.extend(once)
    return ret


def mk_combine(num: int, max_order: int) -> List[MultiResidue]:
    """
    generate combination of order 1 to max_order 
    for [1, 2, ..., n] input

    Parameters
    ----------
    num: int
        input elements range

    max_order: int
        maximum input for combinations

    Returns
    -------
    comb : List[MultiResidue]
        result of combinations

    Examples
    --------
    >> assert(mk_combine(3, 2) == [(0, ), (1, ), (2, ), (0, 1), (0, 2), (1, 2)])
    """
    ret: List[MultiResidue] = []
    for order in range(1, max_order + 1):
        once: List[MultiResidue] = list(
            combinations(list(range(num)), order))
        ret.extend(once)
    return ret


def mk_combine_with_k(
    k: MultiResidue, num: int, max_order: int
) -> List[MultiResidue]:
    """
    connect the combination of (1, 2, ...,n) from order 1 to max_order, which must contains k.

    Parameters
    ----------
    k: int
        must contain element

    n: int
        input elements range

    max_order: int
        maximum input for combinations

    Returns
    -------
    comb : List[MultiResidue]
        result of combinations

    Examples
    --------
    >> assert(mk_combine(3, 2) == [(0, ), (1, ), (2, ), (0, 1), (0, 2), (1, 2)])
    """
    k_set = set(k)
    ret: List[MultiResidue] = []
    for order in range(0, max_order):
        other: List[List[int]] = [
            list(res)
            for res in combinations([j for j in range(num) if j not in k_set], order)
        ]
        for res in other:
            for add in k:
                insort(res, add)
        once: List[MultiResidue] = [tuple(res) for res in other]
        ret.extend(once)
    return ret


def get_epi_from_diff(
    diff: Dict[SeqDiff, float], possiable_keys: List[Seq],
) -> EpiResidue:
    """
    calaulate averaging epistasis value from epistasis delta
    by Graph

    Parameters
    ----------
    diff: Dict[SeqDiff, float]
        epistasis delta of each variance combination

    possiable_keys: List[Seq]
        list of variance combination

    Returns
    -------
    epi_values : EpiResidue
        averaging epistasis value of all variance combination in each residue
    """
    graph = nx.Graph()
    keys_index = {key: i for i, key in enumerate(possiable_keys)}
    epi_values = {key: 0.0 for key in possiable_keys}

    graph.add_edges_from(
        [(keys_index[diff[0]], keys_index[diff[1]]) for diff in list(diff.keys())])
    rings: Generator[Set[int], None, None] = nx.connected_components(graph)

    for ring in rings:
        edges = list(nx.bfs_edges(graph.subgraph(ring), 0))
        for edge in edges:
            edge = cast(Tuple[int, int], edge)
            src_seq, tgt_seq = possiable_keys[edge[0]], possiable_keys[edge[1]]
            # one of (src_seq, tgt_seq) / (tgt_seq, src_seq) must in diff
            delta = (
                diff[(src_seq, tgt_seq)]
                if (src_seq, tgt_seq) in diff
                else -diff[(tgt_seq, src_seq)]
            )
            epi_values[tgt_seq] = epi_values[src_seq] + delta
        select_values = [epi_values[possiable_keys[v]] for v in ring]
        value_mean = sum(select_values) / len(select_values)
        for point in ring:
            epi_values[possiable_keys[point]] -= value_mean
    return epi_values


def fetch_lower_select(
    lower_multi_res: MultiResidue, sorted_at_key: MultiResidue
) -> MultiResidue:
    """
    pick lower_multi_res by sequence of sorted_at_key

    Parameters
    ----------
    lower_multi_res: MultiResidue
        pick source of residues

    sorted_at_key: MultiResidue
        key to be sorted and fetch index

    Returns
    -------
    epi_values : EpiResidue
        averaging epistasis value of all variance combination in each residue
    """
    index = argsort(sorted_at_key)
    res_to_ind = dict(zip(sorted_at_key, index))
    return tuple(res_to_ind[i] for i in lower_multi_res)


def argsort(seq):
    """sort and return index"""
    # http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
    return sorted(range(len(seq)), key=seq.__getitem__)
