"""
Slightly modified subroutines from package ASTROALIGN by Martin Beroiz.
"""

import numpy as np
from numpy.linalg import norm
from scipy.spatial import KDTree
from itertools import combinations

def unique_triangles(points,NUM_NEAREST_NEIGHBORS=5):
    """
    Generates unique triangles from a set of points by finding the nearest neighbors for each point,
    then calculating the ratios of the sides of these triangles. This helps in identifying invariant features
    of triangles that can be used for matching purposes.

    Inputs:
        points -> [array-like] A numpy array of shape (n, 2), representing n points in 2D space.
        NUM_NEAREST_NEIGHBORS -> [int,optional,default=5] Number of nearest neighbors to consider for triangle formation.
    Outputs:
        inv_uniq -> [array-like] A numpy array containing the ratios [L3/L2, L2/L1] for each unique triangle, where L3 is the
      longest side, L2 is the middle, and L1 is the shortest side of the triangle.
        triang_vrtx_uniq -> [array-like] A numpy array of vertex indices [A, B, C] for each triangle, sorted such that A is
      opposite L3, B is opposite L1, and C is opposite L2.
    """
    # Construct a KDTree for efficient nearest neighbor search.
    tree = KDTree(points)

    # Use a set to track unique triangle combinations based on their vertex indices.
    triangles = set()

    # For each point, query the KDTree to find its k nearest neighbors, including itself.
    for i in range(len(points)):
        # Generate all possible combinations of 3 points (triangles) from these neighbors.
        _, indices = tree.query(points[i], k=NUM_NEAREST_NEIGHBORS+1)
        for combo in combinations(indices, 3):
            # Add each triangle to the set, ensuring uniqueness.
            triangles.add(tuple(sorted(combo)))

    # Lists to hold the output ratios and vertex indices.
    inv_uniq = []
    triang_vrtx_uniq = []

    # Process each unique triangle.
    for tri in triangles:
        # Extract the coordinates of the triangle's vertices.
        tri_points = points[list(tri)]
        # Calculate the lengths of each side of the triangle and associate them with the
        # index of the vertex opposite each side. This creates a mapping of side lengths to vertices.
        sides = [
            (norm(tri_points[0] - tri_points[1]), tri[2]),  # Side opposite to vertex C
            (norm(tri_points[0] - tri_points[2]), tri[1]),  # Side opposite to vertex B
            (norm(tri_points[1] - tri_points[2]), tri[0])   # Side opposite to vertex A
        ]
        # Sort the sides by their lengths to identify L1, L2, L3, along with their corresponding opposite vertices.
        sides.sort(key=lambda x: x[0])

        # Calculate the ratios of the side lengths: L3/L2 and L2/L1.
        ratios = [sides[2][0] / sides[1][0], sides[1][0] / sides[0][0]]  # L3/L2, L2/L1
        inv_uniq.append(ratios)
        
        # Sort the vertices according to the problem specification: A opposite L3, B opposite L1, and C opposite L2.
        # This is done by retrieving the indices of the vertices in the order of the side lengths they are opposite to.
        indices_sorted_by_length = [sides[2][1], sides[0][1], sides[1][1]]
        triang_vrtx_uniq.append(indices_sorted_by_length)

    return np.array(inv_uniq), np.array(triang_vrtx_uniq)      