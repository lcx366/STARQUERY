"""
Some slightly modified subroutines of package ASTROALIGN from Martin Beroiz.
"""

import numpy as _np
from scipy.spatial import KDTree
from itertools import combinations
from functools import partial

# The number of nearest neighbors of a given star(including itself) to construct the triangle invariants.
NUM_NEAREST_NEIGHBORS = 6

def _invariantfeatures(x1, x2, x3):
    """
    Given 3 points x1, x2, x3, return the invariant features for the set.
    """
    sides = _np.sort(
        [
            _np.linalg.norm(x1 - x2),
            _np.linalg.norm(x2 - x3),
            _np.linalg.norm(x1 - x3),
        ]
    )
    return [sides[2] / sides[1], sides[1] / sides[0]]

def _arrangetriplet(sources, vertex_indices):
    """
    Return vertex_indices ordered in an (a, b, c) form where:
      a is the vertex defined by L1 & L2
      b is the vertex defined by L2 & L3
      c is the vertex defined by L3 & L1
    and L1 < L2 < L3 are the sides of the triangle defined by vertex_indices.
    """
    ind1, ind2, ind3 = vertex_indices
    x1, x2, x3 = sources[vertex_indices]

    side_ind = _np.array([(ind1, ind2), (ind2, ind3), (ind3, ind1)])
    side_lengths = list(map(_np.linalg.norm, (x1 - x2, x2 - x3, x3 - x1)))
    l1_ind, l2_ind, l3_ind = _np.argsort(side_lengths)

    # the most common vertex in the list of vertices for two sides is the point at which they meet.
    from collections import Counter

    count = Counter(side_ind[[l1_ind, l2_ind]].flatten())
    a = count.most_common(1)[0][0]
    count = Counter(side_ind[[l2_ind, l3_ind]].flatten())
    b = count.most_common(1)[0][0]
    count = Counter(side_ind[[l3_ind, l1_ind]].flatten())
    c = count.most_common(1)[0][0]

    return _np.array([a, b, c])

def _generate_invariants(sources):
    """
    Return an array of (unique) invariants derived from the array `sources`.
    Return an array of the indices of `sources` that correspond to each
    invariant, arranged as described in _arrangetriplet.
    """
    arrange = partial(_arrangetriplet, sources=sources)

    inv = []
    triang_vrtx = []
    coordtree = KDTree(sources)
    # The number of nearest neighbors to request (to work with few sources)
    knn = min(len(sources), NUM_NEAREST_NEIGHBORS)
    for asrc in sources:
        __, indx = coordtree.query(asrc, knn)

        # Generate all possible triangles with the 6 indx provided, and store
        # them with the order (a, b, c) defined in _arrangetriplet
        all_asterism_triang = [
            arrange(vertex_indices=list(cmb)) for cmb in combinations(indx, 3)
        ]
        triang_vrtx.extend(all_asterism_triang)

        inv.extend(
            [
                _invariantfeatures(*sources[triplet])
                for triplet in all_asterism_triang
            ]
        )

    # Remove here all possible duplicate triangles
    uniq_ind = [
        pos for (pos, elem) in enumerate(inv) if elem not in inv[pos + 1 :]
    ]
    inv_uniq = _np.array(inv)[uniq_ind]
    triang_vrtx_uniq = _np.array(triang_vrtx)[uniq_ind]

    return inv_uniq, triang_vrtx_uniq