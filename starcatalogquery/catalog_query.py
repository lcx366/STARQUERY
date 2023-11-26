import os
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord

def from_cap(theta, clat, clon, lmax):
    """
    Creates a mask for a spherical cap on a grid.

    Usage:
        >>> mask_cap = from_cap(theta,clat,clon,lmax)
    Inputs:
        theta -> [float] Angular radius of the spherical cap in degrees.
        clat -> [float] Latitude of the center of the spherical cap in degrees.
        clon -> [float] Longitude of the center of the spherical cap in degrees.
        lmax -> [int] Maximum spherical harmonic degree resolvable by the grid.
    Outputs:
        mask_cap -> [2D array-like] A grid mask with the spherical cap area masked.
    """

    # Calculate the grid size based on lmax
    step = 90 / (lmax + 1)
    
    # Generate arrays of latitudes and longitudes based on the grid size
    lats = np.deg2rad(np.arange(90, -90, -step))
    lons = np.deg2rad(np.arange(0, 360, step))
    
    # Initialize an array to represent the grid mask
    nlat, nlon = len(lats), len(lons)
    mask_cap = np.zeros((nlat, nlon))

    # Convert angular parameters to radians
    theta = np.deg2rad(theta)
    clat = np.deg2rad(clat)
    clon = np.deg2rad(clon)

    # Identify the index range for latitudes within the cap
    imin, imax = np.inf, 0
    for i, lat in enumerate(lats):
        if lat <= clat + theta:
            imin = min(imin, i)
        if lat >= clat - theta:
            imax = max(imax, i)

    # Calculate cartesian coordinates for the cap center
    x = np.cos(clat) * np.cos(clon)
    y = np.cos(clat) * np.sin(clon)
    z = np.sin(clat)

    # Precompute cosine and sine for longitudes
    coslon = np.cos(lons)
    sinlon = np.sin(lons)
    costheta = np.cos(theta)

    # Iterate over the grid and set values within the cap
    for i in range(imin, imax + 1):
        coslat = np.cos(lats[i])
        sinlat = np.sin(lats[i])
        for j in range(nlon):
            # Calculate the angular distance to the cap center
            dist = coslat * (x * coslon[j] + y * sinlon[j]) + z * sinlat
            # Mark points within the cap
            if dist >= costheta: mask_cap[i, j] = 1

    return mask_cap

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
    radius += tile_size * 1.5

    # Determine the grid resolution
    N = int(180 / tile_size)
    if N % 2:
        raise Exception('tile_size must ensure 180/tile_size is even.')

    # Calculate mask for the spherical cap
    L = int(N / 2) - 1
    mask_cap = from_cap(radius, dec_c, ra_c, L)

    # Find indices within the masked area
    dec_index, ra_index = np.where(mask_cap)

    # Compute the sequence indices for the masked tiles
    seqs = dec_index * count_lon + ra_index

    return seqs  

def _load_files(sc_indices, sc_path, sc_name, _mode):
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
            yield pd.read_csv(filename, skiprows=skiprows, dtype=str)

def search_box_magpm(radec_box, sc_path, sc_name, tile_size, mag_threshold, t_pm):
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
    df = pd.concat(_load_files(sc_indices, sc_path, sc_name,'reduced'))
    
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

def search_cone_magpm(center, radius, sc_path, sc_name, tile_size, mag_threshold, t_pm):
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
    df = pd.concat(_load_files(sc_indices, sc_path, sc_name, 'reduced'))

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

def search_box(radec_box, sc_path, sc_name, tile_size):
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
    df = pd.concat(_load_files(sc_indices, sc_path, sc_name, 'simplified'))

    # Remove duplicates and convert coordinates to numeric
    df.drop_duplicates(subset=['ra', 'dec'], inplace=True)
    df[['ra', 'dec', 'mag']] = df[['ra', 'dec', 'mag']].apply(pd.to_numeric)

    # Filter stars within the search box
    ra_min, dec_min, ra_max, dec_max = radec_box
    ra_flag = np.abs(df['ra'] - (ra_min + ra_max) / 2) < (ra_max - ra_min) / 2
    dec_flag = np.abs(df['dec'] - (dec_min + dec_max) / 2) < (dec_max - dec_min) / 2
    df = df[ra_flag & dec_flag].sort_values(by='mag')

    return df.reset_index(drop=True)

def search_cone(center, radius, sc_path, sc_name, tile_size):
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
    df = pd.concat(_load_files(sc_indices, sc_path, sc_name, 'simplified'))
    df.drop_duplicates(subset=['ra', 'dec'], inplace=True)
    df[['ra', 'dec', 'mag']] = df[['ra', 'dec', 'mag']].apply(pd.to_numeric)

    # Filter stars within the cone search area
    c1 = SkyCoord(df['ra'], df['dec'], unit='deg')
    c2 = SkyCoord(ra_c, dec_c, unit='deg')
    sep = c1.separation(c2).deg
    df = df[sep < radius].sort_values(by='mag')

    return df.reset_index(drop=True)     