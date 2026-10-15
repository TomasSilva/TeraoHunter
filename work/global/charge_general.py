"""Planar-graph prediction for jumps that are critically fixed maps, any number of points.
A multiarrangement of s points with sum(k_i - 1) = 2D - 2 jumps in degree D exactly at the critically fixed
maps of degree D with critical points the k_i >= 2 (local degree k_i) and the k_i = 1 points among the
non-critical fixed points (TC_GLOBAL.md, Theorem B and section 5, step 2)."""
import sys
from math import prod
from charge_graphs import multigraphs, plane_maps

def predicted_general(k, D):
    val = [x - 1 for x in k if x >= 2]
    c = len(val)
    assert sum(val) == 2 * D - 2
    nmaps = sum(plane_maps(val, g) for g in multigraphs(val))
    free = D + 1 - c
    need = len(k) - c
    if need > free:
        return 0
    choices = prod(range(free - need + 1, free + 1)) if need > 0 else 1
    if c == 2:
        return nmaps * choices // (D - 1)
    return nmaps * choices

if __name__ == '__main__':
    from itertools import combinations_with_replacement
    # 6-point lines of the 17-line window, jump to (6, 10): D = 6
    vecs = sorted({tuple(sorted(p, reverse=True)) for p in combinations_with_replacement(range(1, 6), 6) if sum(p) == 16},
                  reverse=True)
    for k in vecs:
        print(list(k), 'predicted (6,10) positions:', predicted_general(list(k), 6), flush=True)
