"""
Slightly modified subroutines from package ASTROALIGN by Martin Beroiz.
"""

import numpy as np
from scipy.spatial import KDTree
from itertools import combinations
from functools import partial
from collections import Counter

# The number of nearest neighbors of a given star(including itself) to construct the triangle invariants.
NUM_NEAREST_NEIGHBORS = 6

def _invariantfeatures(x1, x2, x3):
    """
    Compute invariant features for a triangle formed by three points.

    Usage:
        >>> ratios = _invariantfeatures(x1, x2, x3)
    Inputs:
        x1, x2, x3 -> [array-like] Coordinates of the triangle's vertices.
    Outputs:
        ratios -> [list] Two ratios of the sorted side lengths of the triangle.
    """
    # Calculate the lengths of the triangle's sides and sort them
    sides = np.sort([
        np.linalg.norm(x1 - x2),
        np.linalg.norm(x2 - x3),
        np.linalg.norm(x1 - x3),
    ])

    # Return the ratios of the longest side to the middle and the middle to the shortest side
    return [sides[2] / sides[1], sides[1] / sides[0]]

def _arrangetriplet(sources, vertex_indices):
    """
    Arrange vertices of a triangle in a specific order based on side lengths.

    Usage:
        >>> vertex_indices_ordered = _arrangetriplet(sources, vertex_indices)
    Inputs:
        sources -> [2D array-like] Source points.
        vertex_indices_ordered -> [array-like] Indices of the vertices in the 'sources'.
    Outputs:
        vertex_indices -> [array-like] Ordered vertex indices in form of (a, b, c), where
            a is the vertex defined by L1 & L2
            b is the vertex defined by L2 & L3
            c is the vertex defined by L3 & L1
            and L1 < L2 < L3 are the sides of the triangle defined by vertex_indices.
    """
    # Extract the vertex indices and corresponding source points
    ind1, ind2, ind3 = vertex_indices
    x1, x2, x3 = sources[vertex_indices]

    # Calculate and sort the sides of the triangle
    side_lengths = list(map(np.linalg.norm, (x1 - x2, x2 - x3, x3 - x1)))
    l1_ind, l2_ind, l3_ind = np.argsort(side_lengths)

    # Find the common vertex in the shortest and middle sides (a)
    # and in the middle and longest sides (b), and in the longest and shortest sides (c)
    side_ind = np.array([(ind1, ind2), (ind2, ind3), (ind3, ind1)])
    count = Counter(side_ind[[l1_ind, l2_ind]].flatten())
    a = count.most_common(1)[0][0]
    count = Counter(side_ind[[l2_ind, l3_ind]].flatten())
    b = count.most_common(1)[0][0]
    count = Counter(side_ind[[l3_ind, l1_ind]].flatten())
    c = count.most_common(1)[0][0]

    return np.array([a, b, c])

def _generate_invariants(sources):
    """
    Generate unique invariants from a set of source points.

    Usage:
        >>> inv_uniq, triang_vrtx_uniq = _generate_invariants(sources)
    Inputs:
        sources -> [array-like] Array of source points.
    Outputs:
        inv_uniq -> [array-like] Unique invariants
        triang_vrtx_uniq-> [array-like] Corresponding triangle vertex indices.
    """
    # Prepare a partial function for arranging triplets
    arrange = partial(_arrangetriplet, sources=sources)

    # Initialize lists for invariants and triangle vertices
    inv = []
    triang_vrtx = []

    # Create a KDTree for efficient nearest neighbor search
    coordtree = KDTree(sources)
    knn = min(len(sources), NUM_NEAREST_NEIGHBORS)  # Number of nearest neighbors to consider

    # Iterate over each source point
    for asrc in sources:
        # Find nearest neighbors and generate all possible triangles
        __, indx = coordtree.query(asrc, knn)
        all_asterism_triang = [arrange(vertex_indices=list(cmb)) for cmb in combinations(indx, 3)]
        triang_vrtx.extend(all_asterism_triang)

        # Generate invariant features for each triangle
        inv.extend([_invariantfeatures(*sources[triplet]) for triplet in all_asterism_triang])

    # Remove duplicate triangles
    uniq_ind = [pos for (pos, elem) in enumerate(inv) if elem not in inv[pos + 1:]]
    inv_uniq = np.array(inv)[uniq_ind]
    triang_vrtx_uniq = np.array(triang_vrtx)[uniq_ind]

    return inv_uniq, triang_vrtx_uniq