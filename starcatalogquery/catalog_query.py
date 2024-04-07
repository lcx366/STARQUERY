import os
import math
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord

def from_cap(grid_size, clat, clon, theta):
    """
    Determines which points in a grid are inside a cone on the surface of a sphere.

    Usage:
        >>> mask_cap = from_cap(grid_size, clat, clon, theta)
    Inputs:
        rid_size -> [int] The grid size in degrees.
        clat -> [float] The latitude of the cap center in degrees.
        clon -> [float] The longitude of the cap center in degrees.
        theta -> [float] Angular radius of the spherical cap in degrees.
    Outputs:
        mask_cap -> [2D array-like] A 2D boolean array representing the grid. 
        Each element is True if the corresponding grid point is inside the cone, and False otherwise.
    """
    # Generate grid points. Latitude ranges from 90 to -90, and longitude from 0 to 360.
    lat_points = np.arange(90, -90, -grid_size)
    lon_points = np.arange(0, 360, grid_size)
    grid_lat, grid_lon = np.meshgrid(lat_points, lon_points, indexing='ij')

    # Convert the cap center coordinates and grid points to radians for calculation.
    clat_rad = math.radians(clat)
    clon_rad = math.radians(clon)
    
    # Calculate the angular distance between the cap center and each grid point.
    angular_distance_cos = \
    np.sin(clat_rad) * np.sin(np.radians(grid_lat)) + \
    np.cos(clat_rad) * np.cos(np.radians(grid_lat)) * \
    np.cos(np.radians(grid_lon) - clon_rad)

    # Determine whether each grid point is inside the cap.
    inside_cone = angular_distance_cos >= np.cos(np.deg2rad(theta))

    return inside_cone    

def seq2radec(seq, tile_size):
    """
    Converts tile file indices to spherical rectangles (RA and DEC coordinates).

    Usage:
        >>> radec_box,box_center = seq2radec(100,3)
        >>> # radec_box,box_center = seq2radec([120,130],5)
    Inputs:
        seq -> [int,array-like] Index or indices of the star catalog tile file.
        tile_size -> [int] Size of the tile in degrees.
    Outputs:
        radec_box -> Spherical rectangles in form of [ra_min, dec_min, ra_max, dec_max] in degrees.
        box_center -> Center of spherical rectangle in form of [ra, dec] in degrees.
    """
    # Ensure seq is a NumPy array for vectorized operations
    if type(seq) is list: seq = np.array(seq)

    # Calculate the number of tiles along the longitude
    count_lon = 360 // tile_size

    # Compute RA and DEC boundaries for each tile
    ra_min = (seq % count_lon) * tile_size
    ra_max = ra_min + tile_size
    codec_min = (seq // count_lon) * tile_size
    codec_max = codec_min + tile_size
    dec_min = 90 - codec_max
    dec_max = 90 - codec_min

    # Stack coordinates to form the RA and DEC boxes
    radec_box = np.stack([ra_min, dec_min, ra_max, dec_max]).T

    # Calculate the center of each box
    box_center = np.stack([(ra_min + ra_max) / 2, (dec_min + dec_max) / 2]).T

    return radec_box, box_center

def box2seqs(radec_box, tile_size):
    """
    Calculate indices of star catalog tiles covering a rectangular search area.

    Usage:
        >>> seqs = box2seqs([350.6,77.7,373.1,85.4],2)
    Inputs:
        radec_box -> [list or array-like] Rectangular search area in form of [ra_min, dec_min, ra_max, dec_max] in degrees.
        tile_size -> [int] Size of each tile in degrees.
    Outputs:
        seqs -> [int, array-like] Array of indices for the star catalog tile files covering the search area.
    """
    # Boundary tolerance to handle edge cases
    tol = 1e-8

    # Extract coordinates from the search box
    ra_min, dec_min, ra_max, dec_max = radec_box

    # Validate the coordinates to ensure proper rectangular area
    if ra_min >= ra_max or dec_min >= dec_max:
        raise ValueError('Rectangular area must start at the lower left corner and end at the upper right corner.')

    # Calculate the CO-DEC
    codec_min = 90 - dec_max
    codec_max = 90 - dec_min

    # Calculate RA indices covering the search area
    ra_seqs = np.arange(int(ra_min // tile_size), int((ra_max - tol) // tile_size) + 1)

    # Calculate DEC indices covering the search area
    codec_seqs = np.arange(int(codec_min // tile_size), int((codec_max - tol) // tile_size) + 1)

    # Calculate the total number of tiles along the longitude
    count_lon = 360 // tile_size

    # Compute the sequence indices for all tiles within the search area
    seqs = codec_seqs[:, None] * count_lon + ra_seqs

    # Adjust for crossing the primary meridian
    seqs += (ra_seqs < 0) * count_lon - (ra_seqs >= count_lon) * count_lon

    # Flatten the array to get a list of indices
    return seqs.flatten()

def cone2seqs(ra_c, dec_c, radius, tile_size):
    """
    Calculates the indices of star catalog tile files covering a conical search area.

    Usage:
        >>> seqs = cone2seqs(100,20,10,2)
    Inputs:
        ra_c -> [float] Right Ascension of the cone's center in degrees.
        dec_c -> [float] Declination of the cone's center in degrees.
        radius -> [float] Angular radius of the cone in degrees.
        tile_size -> [int] Size of each tile in degrees.
    Outputs:
        seqs -> [int, array-like] Indices of the star catalog tile files covering the cone search area.
    """
    # Calculate the number of tiles along longitude
    count_lon = 360 // tile_size

    # Increase the radius slightly to ensure complete coverage
    radius += tile_size * 1.2

    if 180 % tile_size:
        raise Exception('180/tile_size must be integer.')

    # Calculate mask for the spherical cap
    mask_cap = from_cap(tile_size, dec_c, ra_c, radius)

    # Find indices within the masked area
    dec_index, ra_index = np.where(mask_cap)

    # Compute the sequence indices for the masked tiles
    seqs = dec_index * count_lon + ra_index

    return seqs  

def _load_files(sc_indices, sc_path, sc_name, _mode, max_num_per_tile=None):
    """
    Generator function for loading multiple star catalog tile files.

    Usage:
        >>> _load_files([10,15,178,3430,8009,10002],'starcatalogs/raw/hygv3.7/res5/','hygv3。7,'raw')
    Inputs:
        sc_indices -> [int,array-like] Indices of the star catalog tile files to load.
        sc_path -> [str] Path of star catalog tile files, e.g., 'starcatalogs/raw/hygv3.7/res5/'.
        sc_name -> [str] Name of the star catalog, e.g., 'hygv3.7'.
        _mode -> [str] Type of star catalogs ('raw', 'reduced', 'simplified').
    Yields:
        pandas.DataFrame: Dataframe loaded from each tile file.
    """
    # Set the number of rows to skip based on the catalog type
    skiprows = 1 if _mode == 'raw' else 0

    # Loop through each index in the star catalog
    for sc_index in sc_indices:
        # Construct the filename for each tile
        filename = f'{sc_path}{sc_name}-{sc_index}.csv'

        # Check if file is large enough to be non-empty (arbitrary threshold set at 25 bytes)
        if os.path.getsize(filename) > 25:
            # Yield a DataFrame for each non-empty file
            df = pd.read_csv(filename, skiprows=skiprows, dtype=str)
            if max_num_per_tile is None:
                yield df
            else:    
                df['mag'] = pd.to_numeric(df['mag'])
                df_sorted = df.sort_values(by='mag').head(max_num_per_tile)
                yield df_sorted

def search_box_magpm(radec_box, sc_path, sc_name, tile_size, mag_threshold, t_pm,max_num_per_tile=None):
    """
    Performs a rectangular search on the reduced star catalog considering magnitude and proper motion.

    Usage:
        >>> df = search_box_magpm([20,30,30,40],'starcatalogs/reduced/hygv3.7/res5/','hygv3.7',5,8,2023.5)
    Inputs:
        radec_box -> [array-like] Rectangular search area in form of [ra_min, dec_min, ra_max, dec_max] in degrees.
        sc_path -> [str] Path of star catalog tile files.
        sc_name -> [str] Name of the star catalog.
        tile_size -> [int] Size of the tile in degrees.
        mag_threshold -> [float] Apparent magnitude limit of the detector.
        t_pm -> [float] Epoch when the search was performed.
    Outputs:
        df -> [DataFrame] Dataframe containing stars within the search area.
    """
    # Extract the boundary from the search box
    ra_min, dec_min, ra_max, dec_max = radec_box

    # Convert the search box to catalog tile file indices
    sc_indices = box2seqs(radec_box, tile_size)

    # Concatenate data from relevant tile files
    df = pd.concat(_load_files(sc_indices, sc_path, sc_name,'reduced',max_num_per_tile))
    
    # Remove duplicate entries based on RA and DEC
    df.drop_duplicates(subset=['ra', 'dec'], inplace=True)

    # Convert relevant columns to numeric for calculations
    if {'pmra', 'pmdec'}.issubset(df.columns): 
        df[['ra', 'dec', 'pmra', 'pmdec', 'mag', 'epoch']] = df[['ra', 'dec', 'pmra', 'pmdec', 'mag', 'epoch']].apply(pd.to_numeric)
    else:    
        df[['ra', 'dec', 'mag', 'epoch']] = df[['ra', 'dec', 'mag', 'epoch']].apply(pd.to_numeric)

    # Filter stars based on the magnitude threshold
    mag_flag = df['mag'] < mag_threshold
    df = df[mag_flag].sort_values(by=['mag'])

    # Calculate the time difference from the epoch
    dt = t_pm - df['epoch']

    # Apply proper motion correction if data is available
    if {'pmra', 'pmdec'}.issubset(df.columns):    
        df['ra'] += df['pmra'] / 3.6e6 * dt   
        df['dec'] += df['pmdec'] / 3.6e6 * dt
    else:
        warnings.warn(f'Proper motion data for stars in catalog {sc_name} are not found.')

    # Filter stars within the search area
    ra, dec = df['ra'], df['dec']
    ra_flag = np.abs(ra - (ra_min + ra_max) / 2) < (ra_max - ra_min) / 2
    dec_flag = np.abs(dec - (dec_min + dec_max) / 2) < (dec_max - dec_min) / 2
    flag = ra_flag & dec_flag 
    df = df[flag]

    # Update the epoch to the search epoch
    df['epoch'] = t_pm
    
    return df.reset_index(drop=True)

def search_cone_magpm(center, radius, sc_path, sc_name, tile_size, mag_threshold, t_pm,max_num_per_tile=None):
    """
    Performs a conical search of stars on the reduced star catalog considering magnitude and proper motion.

    Usage:
        >>> df = search_cone([20,30],10,'starcatalogs/reduced/hygv3.7/res5/','hygv3.7',5,8,2023.5)
    Inputs:
        center -> [tuple] Center of the cap in form of [Ra, Dec] in degrees.
        radius -> [float] Angular radius of the cap in degrees.
        sc_path -> [str] Path of star catalog tile files.
        sc_name -> [str] Name of the star catalog.
        tile_size -> [int] Tile size in degrees.
        mag_threshold -> [float] Apparent magnitude limit.
        t_pm -> [float] Epoch when the search was performed.
    Outputs:
        df -> [DataFrame] Dataframe of stars within the cone search area.
    """
    # Extract center coordinates for the cone search
    ra_c, dec_c = center

    # Get tile file indices for the cone search area
    sc_indices = cone2seqs(ra_c, dec_c, radius, tile_size)

    # Load and concatenate data from relevant tile files
    df = pd.concat(_load_files(sc_indices, sc_path, sc_name,'reduced',max_num_per_tile))

    # Remove duplicate entries based on RA and DEC
    df.drop_duplicates(subset=['ra', 'dec'], inplace=True)

    # Convert relevant columns to numeric for calculations
    if {'pmra', 'pmdec'}.issubset(df.columns): 
         df[['ra','dec','pmra','pmdec','mag','epoch']] = df[['ra', 'dec','pmra','pmdec','mag','epoch']].apply(pd.to_numeric)
    else:    
        df[['ra','dec','mag','epoch']] = df[['ra','dec','mag','epoch']].apply(pd.to_numeric)

    # Filter stars based on magnitude threshold
    df = df[df['mag'] < mag_threshold].sort_values(by='mag')

    # Proper motion adjustment, if available
    dt = t_pm - df['epoch']
    if 'pmra' in df.columns and 'pmdec' in df.columns:
        df['ra'] += df['pmra'] / 3.6e6 * dt
        df['dec'] += df['pmdec'] / 3.6e6 * dt
    else:
        warnings.warn(f'Proper motion data for stars in catalog {sc_name} not found.')

    # Filter stars within the cone area
    c1 = SkyCoord(df['ra'], df['dec'], unit='deg')
    c2 = SkyCoord(ra_c, dec_c, unit='deg')
    sep = c1.separation(c2).deg
    df = df[sep < radius]

    # Update the epoch for the search results
    df['epoch'] = t_pm

    return df.reset_index(drop=True)  

def search_box(radec_box, sc_path, sc_name, tile_size,max_num_per_tile=None):
    """
    Performs a rectangular search on the simplified star catalog without considering magnitude and proper motion.

    Usage:
        >>> df = search_box([20,30,30,40],'starcatalogs/simplified/hygv3.7/res5/mag8.0/epoch2023.0/','hygv3.7',5)
    Inputs:
        radec_box -> [array-like] Rectangular search area in form of [ra_min, dec_min, ra_max, dec_max] in degrees.
        sc_path -> [str] Path of star catalog tile files.
        sc_name -> [str] Name of the star catalog.
        tile_size -> [int] Size of the tile in degrees.
    Outputs:
        df -> [DataFrame] Dataframe of stars within the rectangular search area.
    """
    # Convert search area to catalog tile indices
    sc_indices = box2seqs(radec_box, tile_size)

    # Load and concatenate data from tile files
    df = pd.concat(_load_files(sc_indices, sc_path, sc_name,'simplified',max_num_per_tile))

    # Remove duplicates and convert coordinates to numeric
    df.drop_duplicates(subset=['ra', 'dec'], inplace=True)
    df[['ra', 'dec', 'mag']] = df[['ra', 'dec', 'mag']].apply(pd.to_numeric)

    # Filter stars within the search box
    ra_min, dec_min, ra_max, dec_max = radec_box
    ra_flag = np.abs(df['ra'] - (ra_min + ra_max) / 2) < (ra_max - ra_min) / 2
    dec_flag = np.abs(df['dec'] - (dec_min + dec_max) / 2) < (dec_max - dec_min) / 2
    df = df[ra_flag & dec_flag].sort_values(by='mag')

    return df.reset_index(drop=True)

def search_cone(center, radius, sc_path, sc_name, tile_size,max_num_per_tile=None):
    """
    Performs a cone search on the simplified star catalog without considering magnitude and proper motion.

    Usage:
        >>> df = search_cone([20,30],10,'starcatalogs/reduced/hygv3.7/res5/','hygv3.7',5)
    Inputs:
        center -> [tuple] Center of the cap in form of [Ra, Dec] in degrees.
        radius -> [float] Angular radius of the cap in degrees.
        sc_path -> [str] Path of star catalog tile files.
        sc_name -> [str] Name of the star catalog.
        tile_size -> [int] Tile size in degrees.
    Outputs:
        df -> [DataFrame] Dataframe of stars within the cone search area.
    """
    # Extract center coordinates and get tile indices
    ra_c, dec_c = center
    sc_indices = cone2seqs(ra_c, dec_c, radius, tile_size)

    # Load and concatenate data from relevant tile files
    df = pd.concat(_load_files(sc_indices, sc_path, sc_name,'simplified',max_num_per_tile))
    df.drop_duplicates(subset=['ra', 'dec'], inplace=True)
    df[['ra', 'dec', 'mag']] = df[['ra', 'dec', 'mag']].apply(pd.to_numeric)

    # Filter stars within the cone search area
    c1 = SkyCoord(df['ra'], df['dec'], unit='deg')
    c2 = SkyCoord(ra_c, dec_c, unit='deg')
    sep = c1.separation(c2).deg
    df = df[sep < radius].sort_values(by='mag')

    return df.reset_index(drop=True)     