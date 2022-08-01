from bisect import insort
from itertools import combinations
from typing import Generator, Union, cast, List, Tuple, Dict, Set

import networkx as nx

from cliff.metadata import SEQ, MULTI_RESIDUE

SEQ_DIFF = Tuple[SEQ, SEQ]
VAL_EACH_RESIDUE = Dict[Tuple[MULTI_RESIDUE], float]
EPI_EACH_RESIDUE = Dict[SEQ, float]
EPI_NET = Dict[MULTI_RESIDUE, EPI_EACH_RESIDUE]
EPI_NET_BASE = Dict[Tuple[MULTI_RESIDUE], EPI_NET]


def select_substr(
    src: Union[str, SEQ, List[str]], select: Union[MULTI_RESIDUE, List[int]]
) -> SEQ:
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
    if type(select) == set:
        select = sorted(List(select))
    select = cast(List[int], select)
    return tuple(src[s] for s in select)


def mk_combine_subset(
    all_res: Tuple[MULTI_RESIDUE],
) -> List[Tuple[MULTI_RESIDUE]]:
    """
    generate combination of multiply residues

    Parameters
    ----------
    all_res: Tuple[MULTI_RESIDUE]
        several residues of comb input

    Returns
    -------
    comb : List[Tuple[MULTI_RESIDUE]]
        result of combinations

    Examples
    --------
    >> assert(mk_combine_subset(((1, 2), (3, 4))) == (((1, 2),), ((3, 4),), ((1, 2), (3, 4))))
    """
    ret: List[Tuple[MULTI_RESIDUE]] = []
    for order in range(1, len(all_res)):
        once: List[Tuple[MULTI_RESIDUE]] = sorted(
            list(combinations(all_res, order))
        )
        ret.extend(once)
    return ret


def mk_combine(n: int, max_order: int) -> List[MULTI_RESIDUE]:
    """
    generate combination of order 1 to max_order 
    for [1, 2, ..., n] input

    Parameters
    ----------
    n: int
        input elements range

    max_order: int
        maximum input for combinations

    Returns
    -------
    comb : List[MULTI_RESIDUE]
        result of combinations

    Examples
    --------
    >> assert(mk_combine(3, 2) == [(0, ), (1, ), (2, ), (0, 1), (0, 2), (1, 2)])
    """
    ret: List[MULTI_RESIDUE] = []
    for order in range(1, max_order + 1):
        once: List[MULTI_RESIDUE] = list(
            combinations([j for j in range(n)], order))
        ret.extend(once)
    return ret


def mk_combine_with_k(
    k: MULTI_RESIDUE, n: int, max_order: int
) -> List[MULTI_RESIDUE]:
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
    comb : List[MULTI_RESIDUE]
        result of combinations

    Examples
    --------
    >> assert(mk_combine(3, 2) == [(0, ), (1, ), (2, ), (0, 1), (0, 2), (1, 2)])
    """
    k_set = set(k)
    ret: List[MULTI_RESIDUE] = []
    for order in range(0, max_order):
        other: List[List[int]] = [
            list(res)
            for res in combinations([j for j in range(n) if j not in k_set], order)
        ]
        for res in other:
            for add in k:
                insort(res, add)
        once: List[MULTI_RESIDUE] = [tuple(res) for res in other]
        ret.extend(once)
    return ret


def get_epi_from_diff(
    diff: Dict[SEQ_DIFF, float], possiable_keys: List[SEQ],
) -> EPI_EACH_RESIDUE:
    """
    calaulate averaging epistasis value from epistasis delta
    by Graph

    Parameters
    ----------
    diff: Dict[SEQ_DIFF, float]
        epistasis delta of each variance combination

    possiable_keys: List[SEQ]
        list of variance combination

    Returns
    -------
    epi_values : EPI_EACH_RESIDUE
        averaging epistasis value of all variance combination in each residue
    """
    G = nx.Graph()
    diff_keys = list(diff.keys())
    keys_index = {key: i for i, key in enumerate(possiable_keys)}
    epi_values = {key: 0.0 for key in possiable_keys}

    edges = [(keys_index[diff[0]], keys_index[diff[1]]) for diff in diff_keys]
    G.add_edges_from(edges)
    rings: Generator[Set[int], None, None] = nx.connected_components(G)

    for ring in rings:
        G_subset = G.subgraph(ring)
        edges = nx.bfs_edges(G_subset, 0)
        edges = list(edges)
        for edge in edges:
            edge = cast(Tuple[int, int], edge)
            src, tgt = edge
            src_seq, tgt_seq = possiable_keys[src], possiable_keys[tgt]
            # one of (src_seq, tgt_seq) / (tgt_seq, src_seq) must in diff
            delta = (
                diff[(src_seq, tgt_seq)]
                if (src_seq, tgt_seq) in diff
                else -diff[(tgt_seq, src_seq)]
            )
            epi_values[tgt_seq] = epi_values[src_seq] + delta
        select_values = [epi_values[possiable_keys[v]] for v in ring]
        value_mean = sum(select_values) / len(select_values)
        for v in ring:
            epi_values[possiable_keys[v]] -= value_mean
    return epi_values


def fetch_lower_select(
    lower_multi_res: MULTI_RESIDUE, sorted_at_key: MULTI_RESIDUE
) -> MULTI_RESIDUE:
    """
    pick lower_multi_res by sequence of sorted_at_key

    Parameters
    ----------
    lower_multi_res: MULTI_RESIDUE
        pick source of residues

    sorted_at_key: MULTI_RESIDUE
        key to be sorted and fetch index

    Returns
    -------
    epi_values : EPI_EACH_RESIDUE
        averaging epistasis value of all variance combination in each residue
    """
    index = argsort(sorted_at_key)
    res_to_ind = dict(zip(sorted_at_key, index))
    return tuple(res_to_ind[i] for i in lower_multi_res)


def argsort(seq):
    # http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
    return sorted(range(len(seq)), key=seq.__getitem__)
