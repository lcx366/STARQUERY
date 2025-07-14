import numpy as np
from scipy.spatial.transform import Rotation
from astropy.wcs import WCS

def Rot(seq, angles, degrees=True):
    """
    Generates a rotation matrix for a sequence of rotations around specified axes.

    Usage:
        >>> seq = 'XY'
        >>> angles = [60,30]
        >>> rotation_matrix = Rot(seq,angles)
    Inputs:
        seq -> [str] Sequence of axes ('X', 'Y', 'Z') for rotation, e.g., 'ZY'.
        angles -> [float or list] Rotation angle(s) in radians (default) or degrees.
        degrees -> [bool,optional,default=True] If True, angles are interpreted as degrees, otherwise radians.
    Returns:
        rotation_matrix -> [2D/3D array-like] Rotation matrix (3x3 or array of 3x3 matrices).
    Note:
        - The function is only suitable for a right-handed reference frame.
        - If multiple angles are provided, they should match the length of 'seq'.
    """
    # Check if 'angles' is a single value and not an array
    if np.isscalar(angles):
        # Generate a single rotation matrix for the given angle and axes
        # 'as_matrix().T' converts the rotation object to a 3x3 matrix and transposes it
        rotation_matrix = Rotation.from_euler(seq, angles, degrees).as_matrix().T
    else:
        # Convert 'angles' to a NumPy array for processing
        angles = np.array(angles)

        # Check if the sequence length is more than one and handle accordingly
        if len(seq) > 1:
            if angles.ndim == 1:
                # For a single set of angles, generate a single 3x3 rotation matrix
                rotation_matrix = Rotation.from_euler(seq, angles, degrees).as_matrix().T
            elif angles.ndim == 2:
                # For multiple sets of angles, generate a series of 3x3 rotation matrices
                rotation_matrix = Rotation.from_euler(seq, angles, degrees).as_matrix().transpose(0, 2, 1)
        else:
            # Handle the case for a single-axis rotation with multiple angle values
            rotation_matrix = Rotation.from_euler(seq, angles, degrees).as_matrix().transpose(0, 2, 1)

    # Return the rotation matrix or matrices
    return rotation_matrix

def wcs_trans(pixel_width, fp_radec):
    """
    Defines the WCS(World Coordinate System) transformation between celestial and pixel coordinates.

    Usage:
        >>> wcs = wcs_trans(0.01,[20,30]):
    Inputs: 
        pixel_width -> [list or tuple of int] Pixel width in degrees. If a single float, it's used for both axes.
        fp_radec -> [list or tuple of int] Fiducial point (center pointing) of the camera in RA and Dec (in degrees).
    Returns:
        wcs -> [astropy.wcs.WCS] WCS object representing the transformation.
    """
    # Handling single or tuple pixel width inputs
    if np.isscalar(pixel_width):
        # If the pixel width is a scalar, use the same value for both axes
        pixel_width_axis1 = pixel_width_axis2 = pixel_width
    else:
        # If the pixel_width is a tuple, assign values to each axis separately
        pixel_width_axis1, pixel_width_axis2 = pixel_width

    # Extracting fiducial point coordinates (RA and DEC) from the input
    fp_ra, fp_dec = fp_radec

    # Setting up the WCS transformation parameters in a dictionary
    wcs_input_dict = {
        'CTYPE1': 'RA---TAN',  # Specifies the first(RA) axis type and projection (tangential)
        'CUNIT1': 'deg',       # Unit for the first coordinate axis (degrees)
        'CDELT1': pixel_width_axis1,  # Pixel scale (width) for the first coordinate axis in degrees
        'CRPIX1': 1,           # Reference pixel location for the first coordinate axis
        'CRVAL1': fp_ra,        # RA value at the reference pixel
        'CTYPE2': 'DEC--TAN',  # Specifies the second(DEC) axis type and projection
        'CUNIT2': 'deg',       # Unit for the second coordinate axis
        'CDELT2': pixel_width_axis2,  # Pixel scale (height) for the second coordinate axis in degrees
        'CRPIX2': 1,           # Reference pixel location for the second coordinate axis
        'CRVAL2': fp_dec       # DEC value at the reference pixel
    }

    # Creating the WCS object with the specified transformation parameters
    wcs = WCS(wcs_input_dict)

    # Returning the WCS object
    return wcs  

def xy_catalog(fp_radec,radec,pixel_width,theta=0):
    """
    Calculates the pixel coordinates of stars given their celestial coordinates.

    Usage:
        >>> x,y = xy_catalog([10,20],[[11,15],[22,-4]],0.01)
    Inputs:
        fp_radec -> [list or tuple of int] Fiducial point of the camera in RA and Dec (in degrees).
        radec -> [2d array of float] Celestial coordinates of stars in form of (RA, Dec) in degrees.
        pixel_width -> [float] Pixel width in degrees.
        theta -> [float,optional,default=0] Rotation angle from WCS frame(equivalent to ENU) to image frame in radians.
    Returns:
        x -> [array of float] x pixel coordinates of the stars      
        y -> [array of float] y pixel coordinates of the stars   
        wcs -> [astropy.wcs.WCS] WCS object representing the transformation.
    """
    # Define the WCS transformation based on pixel width and fiducial point
    wcs = wcs_trans(pixel_width,fp_radec)

    # Directly convert world coordinates (RA/DEC) to pixel coordinates
    x, y = wcs.world_to_pixel_values(radec[:, 0], radec[:, 1])

    # Check if rotation is needed
    if theta != 0:
        # Apply rotation on pixel coordinates
        # 'Rot' function is called to get the rotation matrix
        x,y = Rot('Z',theta,degrees=False)[:-1,:-1] @ np.stack([x,y])

    return x,y,wcs   