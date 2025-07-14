import numpy as np
from numpy.linalg import norm

TWOPI = 2 * np.pi

def unit_vector(vector):
    """
    Normalize a vector or an array of vectors to unit vectors.
    Supports 1D or ND array, normalizes along the last axis.

    Usage:
        >>> print(unit_vector([1, 2, 3]))
        [0.26726124 0.53452248 0.80178373]
        >>> print(unit_vector([[1, 2, 3], [4, 5, 6]]))
        [[0.26726124 0.53452248 0.80178373]
         [0.45584231 0.56980288 0.68376346]]

    Inputs:
        vector -> [array-like] The input vector(s).
    Outputs:
        unit_vector -> [array-like] The normalized unit vector(s).
    """
    norm_vector = norm(vector, axis=-1, keepdims=True)
    return vector / norm_vector

def matrix_dot_vector(matrix, vector):
    """
    Computes the dot product of matrices and vectors.

    Usage:
        >>> import numpy as np
        >>> matrix = np.arange(1296).reshape(8, 18, 3, 3)
        >>> vector = np.arange(432).reshape(8, 18, 3)
        >>> vector_trans = matrix_dot_vector(matrix, vector)

    Inputs:
        matrix -> [array-like] Multi-dimensional matrix.
        vector -> [array-like] Multi-dimensional vector.
    Outputs:
        vector_trans -> [array-like] Transformed vectors.
    """
    vector = np.asarray(vector)
    return np.matmul(matrix, vector[..., None]).squeeze(-1)

def spherical_to_cartesian(ra, dec, r, degrees=True):
    """
    Convert spherical coordinates (Right Ascension, Declination, Range) to Cartesian coordinates.

    Usage:
        >>> spherical_to_cartesian(45, 45, 1)
        array([0.5      , 0.5      , 0.70710678])
        >>> spherical_to_cartesian([45, 90], [45, 45], [1, 1])
        array([[5.00000000e-01, 5.00000000e-01, 7.07106781e-01],
               [3.74939946e-33, 7.07106781e-01, 7.07106781e-01]])
    Inputs:
        ra -> [array-like] Right Ascension.
        dec -> [array-like] Declination.
        r -> [array-like] Distance.
        degrees -> [bool, optional, default=True] Specifies if RA and Dec are in degrees or radians.
    Outputs:
        xyz -> [array-like] Cartesian coordinates as an array of shape (3,) or (N, 3).
    """
    if degrees:
        ra = np.radians(ra)
        dec = np.radians(dec)

    x = r * np.cos(dec) * np.cos(ra)
    y = r * np.cos(dec) * np.sin(ra)
    z = r * np.sin(dec)

    return np.stack([x, y, z], axis=-1)


def cartesian_to_spherical(x, y, z, degrees=True):
    """
    Convert Cartesian coordinates to spherical coordinates (Right Ascension, Declination, Distance).

    Usage:
        >>> cartesian_to_spherical(0.5, 0.5, 0.70710678)
        array([ 45.,  45.,   1.])
        >>> cartesian_to_spherical([0.5, 0], [0.5, 0], [0.70710678, 1])
        array([[ 45.,  45.,   1.],
               [  0.,  90.,   1.]])
    Inputs:
        x -> [float or array-like] X coordinate.
        y -> [float or array-like] Y coordinate.
        z -> [float or array-like] Z coordinate.
        degrees -> [bool, optional, default=True] Specifies if RA and Dec should be returned in degrees or radians.
    Outputs:
        ra_dec_r -> [array-like] Spherical coordinates as an array of shape (3,) or (N, 3).
    """
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    ra = np.arctan2(y, x) % TWOPI # normalize RA to [0, 2pi)
    dec = np.arcsin(z / r)

    if degrees:
        ra = np.degrees(ra)
        dec = np.degrees(dec)

    return np.stack([ra, dec, r], axis=-1)

def separation(ra1, dec1, ra2, dec2):
    """
    Calculate the cosine of the spherical angular separation between two points.

    Usage:
        >>> separation(np.radians(0), np.radians(0), np.radians(90), np.radians(0))
        6.123233995736766e-17
        >>> separation(np.radians([0, 45]), np.radians([0, 45]), np.radians([90, 45]), np.radians([0, 45]))
        array([6.12323400e-17, 5.00000000e-01])
    Inputs:
        ra1, dec1 -> [float or array-like] Right Ascension and Declination of the first point (in radians).
        ra2, dec2 -> [float or array-like] Right Ascension and Declination of the second point (in radians).
    Outputs:
        angular_distance_cos -> [float or array-like] Cosine of the angular separation.
    """
    angular_distance_cos = (
            np.sin(dec1) * np.sin(dec2) +
            np.cos(dec1) * np.cos(dec2) *
            np.cos(ra1 - ra2)
    )

    return angular_distance_cos


def axes_to_antisymmetric_matrices(axes):
    """
    Convert rotation axes to corresponding antisymmetric matrices.
    When a vector rotates around an axis by a small angle θ in a clockwise direction, the rotation matrix can be approximated as
    I - θ * K, where I is the identity matrix and K is the antisymmetric matrix corresponding to the rotation axis.
    The antisymmetric matrix K for a rotation axis [ax, ay, az] is defined as:
        K = [[  0, -az,  ay],
             [ az,   0, -ax],
             [-ay,  ax,   0]]

    Usage:
        >>> axes = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        >>> axes_to_antisymmetric_matrices(axes)
        array([[[ 0.,  0.,  0.],
                [ 0.,  0., -1.],
                [ 0.,  1.,  0.]],

               [[ 0.,  0.,  1.],
                [ 0.,  0.,  0.],
                [-1.,  0.,  0.]],

               [[ 0., -1.,  0.],
                [ 1.,  0.,  0.],
                [ 0.,  0.,  0.]]])

    Inputs:
        axes -> [numpy.ndarray] A 2D array of shape (n, 3), where each row is a rotation axis [ax, ay, az].
    Outputs:
        antisymmetric_matrices -> [numpy.ndarray] A 3D array of shape (n, 3, 3), where each item is a 3x3 antisymmetric matrix.
    """
    n = axes.shape[0]
    antisymmetric_matrices = np.zeros((n, 3, 3))

    antisymmetric_matrices[:, 0, 1] = -axes[:, 2]
    antisymmetric_matrices[:, 0, 2] = axes[:, 1]
    antisymmetric_matrices[:, 1, 0] = axes[:, 2]
    antisymmetric_matrices[:, 1, 2] = -axes[:, 0]
    antisymmetric_matrices[:, 2, 0] = -axes[:, 1]
    antisymmetric_matrices[:, 2, 1] = axes[:, 0]

    return antisymmetric_matrices
