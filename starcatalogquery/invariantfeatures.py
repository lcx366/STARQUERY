import numpy as np
from numpy.linalg import norm
from scipy.spatial import KDTree
from itertools import combinations

def unique_quads(points, NUM_NEAREST_NEIGHBORS=10):
    """
    Generates unique quads from a set of points by finding the nearest neighbors for each point,
    then calculating invariant features of quads.

    Usage:
        >>> points = np.random.rand(100, 2)  # Generate 100 random 2D points
        >>> invariants, quads = unique_quads(points)
    Inputs:
        points -> [array-like] A numpy array of shape (n, 2), representing n points in 2D space.
        NUM_NEAREST_NEIGHBORS -> [int, optional, default=8] Number of nearest neighbors to consider for quad formation.
    Outputs:
        inv_uniq -> [array-like] A numpy array containing the invariant features of each unique quad.
        quad_vrtx_uniq -> [array-like] A numpy array of vertex indices for each quad.
    """
    # Construct a KDTree for efficient nearest neighbor search.
    tree = KDTree(points)

    # Use a set to track unique quad combinations based on their vertex indices.
    quads = set()

    n = len(points)
    if n < 4:
        raise ValueError('Number of points must be >= 4 for quads.')
    elif n < NUM_NEAREST_NEIGHBORS+1:
        NUM_NEAREST_NEIGHBORS = n - 1

    # For each point, query the KDTree to find its k nearest neighbors, including itself.
    for i in range(n):
        # Generate all possible combinations of 4 points (quads) from these neighbors.
        _, indices = tree.query(points[i], NUM_NEAREST_NEIGHBORS + 1)
        for combo in combinations(indices, 4):
            # Add each quad to the set, ensuring uniqueness.
            quads.add(tuple(sorted(combo)))

    # Lists to hold the output quads and vertex indices.
    inv_uniq = []
    quad_vrtx_uniq = []

    # Process each unique quad.
    for quad in quads:
        quad_indices = np.array(quad)
        # Extract the coordinates of the quad's vertices.
        quad_points = points[quad_indices]

        # Calculate the distances between any two points.
        diff = quad_points[:, np.newaxis, :] - quad_points[np.newaxis, :, :]
        distances = norm(diff, axis=-1)

        # Find the two points with the longest distance.
        A_index, B_index = np.unravel_index(np.argmax(distances, axis=None), distances.shape)
        C_index, D_index = list(set(range(4)) - {A_index, B_index})
        A, B, C, D = quad_points[[A_index, B_index, C_index, D_index]]

        # Establish a coordinate system with A as the origin.
        origin = A
        x_axis = B - A
        length_scale = norm(x_axis)

        x_axis = x_axis / length_scale  # Unit vector

        # Rotate clockwise by 45 degrees.
        x_axis = np.array([[1, 1], [-1, 1]]) @ x_axis * np.sqrt(2) / 2

        # Compute the y-axis perpendicular to the x-axis.
        y_axis = np.array([-x_axis[1], x_axis[0]])

        # Normalize coordinates of C and D.
        C = (C - A) / length_scale * np.sqrt(2)
        D = (D - A) / length_scale * np.sqrt(2)

        base_matrix = np.vstack([x_axis, y_axis])

        C = base_matrix @ C
        D = base_matrix @ D

        if C[0] + D[0] > 1:
            A_index, B_index = B_index, A_index
            C, D = 1 - C, 1 - D

        if C[0] > D[0]:
            C_index, D_index = D_index, C_index
            C, D = D, C

        inv_uniq.append(np.hstack([C, D]))
        quad_indices = quad_indices[[A_index, B_index, C_index, D_index]]
        quad_vrtx_uniq.append(quad_indices)

    return np.array(inv_uniq), np.array(quad_vrtx_uniq)

def unique_triangles(points,NUM_NEAREST_NEIGHBORS=10):
    """
    Generates unique triangles from a set of points by finding the nearest neighbors for each point,
    then calculating the ratios of the sides of these triangles.

    Usage:
        >>> points = np.random.rand(100, 2)  # Generate 100 random 2D points
        >>> invariants, triangles = unique_triangles(points)
    Inputs:
        points -> [array-like] A numpy array of shape (n, 2), representing n points in 2D space.
        NUM_NEAREST_NEIGHBORS -> [int, optional, default=5] Number of nearest neighbors to consider for triangle formation.
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

    n = len(points)
    if n < 3:
        raise ValueError('Number of points must be >= 3 for triangles.')
    elif n < NUM_NEAREST_NEIGHBORS+1:
        NUM_NEAREST_NEIGHBORS = n - 1

    # For each point, query the KDTree to find its k nearest neighbors, including itself.
    for i in range(len(points)):
        # Generate all possible combinations of 3 points (triangles) from these neighbors.
        _, indices = tree.query(points[i], NUM_NEAREST_NEIGHBORS+1)
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

def calculate_invariantfeatures(sources,mode_invariants):
    """
    Computes geometric invariant features for a set of source points by generating unique triangles or quads.
    These features are then used to build a KDTree for efficient matching.

    Inputs:
        sources -> [array-like] A 2D numpy array representing the pixel coordinates of source points.
        mode_invariants -> [str] Mode of geometric invariants to use. Available options are 'triangles' or 'quads'.
    Outputs:
        invariants -> [array-like] A 2D numpy array containing the invariant ratios (L2/L1, L1/L0) for each triangle, where L2 >= L1 >= L0 are the sides of the triangle.
        asterisms -> [array-like] A 2D numpy array containing the indices of source points that form each triangle.
        kdtree -> An instance of scipy.spatial.KDTree built from the invariants for quick nearest-neighbor lookup.
    """
    # Derive unique geometric invariants
    if mode_invariants == 'triangles':
        invariants,asterisms = unique_triangles(sources)
    elif mode_invariants == 'quads':
        invariants,asterisms = unique_quads(sources)
    else:
        raise ValueError('mode of invariant features must be either "triangles" or "quads"')
    # Construct a KDTree structure using the unique invariants
    kdtree = KDTree(invariants)
    return invariants,asterisms,kdtree