import os
import numpy as np
import pandas as pd
import healpy as hp
from astropy.time import Time

from .catalog_index import find_healpix_level,index_query_sql
from .utils.math import separation
from .astrometry_corrections import apply_astrometry_corrections

def create_empty_df(file_path, skip_initial_description=False):
    """
    Creates an empty DataFrame using column headers read from a specified file.

    Inputs:
        file_path -> [str] Path to the CSV file from which to read the column headers.
        skip_initial_description -> [bool] Indicates whether to skip an initial description line (e.g., for 'raw' mode).
    Outputs:
        empty_df -> [pd.DataFrame] An empty DataFrame with column names set as per the file's header.
    Raises:
        FileNotFoundError: If the specified file does not exist.
        ValueError: If the file is empty or headers cannot be determined.
    """
    try:
        # Open the file and read the first line to get the headers
        with open(file_path, 'r') as file:
            if skip_initial_description:
                description = next(file)  # Skip the first line if it's a description
            header_line = next(file).strip()  # Read the next line as the header
            headers = header_line.split(',')
    except FileNotFoundError:
        raise FileNotFoundError(f"The file at {file_path} does not exist.")
    except StopIteration:
        raise ValueError(f"The file at {file_path} is empty or does not have headers.")

    # Create and return an empty DataFrame with these headers
    return pd.DataFrame(columns=headers)

def _load_files(K5_indices_list, dir_sc, sc_name, _mode, max_num_per_tile=None):
    """
    Generator function that yields star catalog data from specified files and rows based on HEALPix indices.

    Inputs:
        K5_indices_list -> [list] Each sublist contains tuples of (file index, [row indices]) specifying which rows to load from each file group.
        dir_sc -> [str] Directory where the star catalog files are stored.
        sc_name -> [str] Base name of the star catalog files, expected to be formatted as "<sc_name>-<index>.csv".
        _mode -> [str] Mode that determines the number of header rows to skip; 'raw' mode skips two rows, other modes skip one row.
        max_num_per_tile -> [int, optional, default=None] Maximum number of rows (stars) to keep for each tile after sorting by magnitude.
    Outputs:
        pd.DataFrame -> DataFrames containing the specified rows from each file group after sorting and limiting by magnitude.
    """

    if not K5_indices_list:
        # If no data was generated throughout the loop, yield an empty DataFrame
        file_path = os.path.join(dir_sc, f"{sc_name}-0.csv")
        yield create_empty_df(file_path, initial_skip)

    # Setup initial conditions based on _mode before processing any files
    initial_skip = 1 if _mode == 'raw' else 0  # Number of initial rows to skip due to mode
    # Process each file group specified in K5_indices_list
    for file_group in K5_indices_list:

        all_dfs = []  # List to hold all dataframes from current group before sorting and limiting
        for file_index, row_indices in file_group:
            file_path = os.path.join(dir_sc, f"{sc_name}-{file_index}.csv")

            if os.path.getsize(file_path) > 25:
                # Determine header row dynamically if the mode requires it
                with open(file_path, 'r') as file:
                    if _mode == 'raw':
                        description = next(file)  # Skip description line only for 'raw' mode
                    header = next(file).strip().split(',')

                # Prepare adjusted row indices considering the initial skip
                adjusted_row_indices = [i + initial_skip + 1 for i in row_indices]  # Adjust for the actual data rows

                # Define the rows to be loaded using adjusted indices
                skip_rows = lambda x: x not in adjusted_row_indices

                # Read the file with proper headers and adjusted skipping
                df = pd.read_csv(file_path, skiprows=skip_rows, header=None, names=header)

                # Add dataframe to list
                all_dfs.append(df)

        # Combine all dataframes from current group
        combined_df = pd.concat(all_dfs, ignore_index=True)

        # Process 'mag' column if necessary
        if max_num_per_tile is not None:
            combined_df['mag'] = pd.to_numeric(combined_df['mag'])
            combined_df = combined_df.sort_values(by='mag').head(max_num_per_tile)

        yield combined_df

def search_box_raw(radec_box, dir_sc, sc_name, _mode, tb_name, catalog_indices_db, mag_threshold, t_pm, fov_min, max_num_per_tile=None):
    """
    This function performs a rectangular search of stars in specified star catalogs within a given RA/Dec box.
    It applies magnitude filtering and proper motion correction based on the input parameters.

    Inputs:
        radec_box -> [list] Rectangular search area in form of [ra_min, dec_min, ra_max, dec_max] in degrees.
        dir_sc -> [str] Directory where the star catalog files are stored.
        sc_name -> [str] Name of the star catalog.
        _mode -> [str] Mode that determines the number of header rows to skip; 'raw' mode skips two rows, other modes skip one row.
        tb_name -> [str] Name of the star catalog table.
        catalog_indices_db -> [str] Path to the database containing the star catalog indices.
        mag_threshold -> [float] Apparent magnitude limit.
        t_pm -> [float] Epoch to which the stars are unified.
        fov_min -> [float] Field of view parameters in degrees. It determines the hierarchical division of the sky region in HEALPix,
        ensuring that each tile's size is greater than one-quarter of the FOV and less than or equal to half of the FOV.
        max_num_per_tile -> [int, optional, default=None] Maximum number of stars to include in each tile, sorted by brightness.
        If None, includes all stars meeting the magnitude criteria.
    Outputs:
        df -> [pd.DataFrame] DataFrame containing the search results.
        level -> [str] HEALPix level used for the search.
        nside -> [int] NSIDE parameter of the HEALPix grid used for the search.
        ids -> [list] List of HEALPix pixel IDs covering the search box.
        pixel_size -> [float] Approximate size of each HEALPix pixel in degrees.
        fov_min -> [float] Field of view parameters in degrees.
    """
    # Extract the boundary from the search box
    ra_min, dec_min, ra_max, dec_max = radec_box
    dec_c,ra_c = np.mean([dec_max,dec_min]),np.mean([ra_max,ra_min])
    ra_range,dec_range = ra_max-ra_min,dec_max-dec_min
    if fov_min is None:
        fov_min = min(ra_range*np.cos(np.radians(dec_c)),dec_range)

    vertices = np.array([
        [ra_min, dec_min],  # Bottom left
        [ra_max, dec_min],  # Bottom right
        [ra_max, dec_max],  # Top right
        [ra_min, dec_max]   # Top left
    ])

    vertices_uec = hp.ang2vec(vertices[:,0],vertices[:,1],lonlat=True)

    # Find the HEALPix parameters for a given field of view (FOV) in degrees
    level, nside, npix, pixel_size = find_healpix_level(fov_min)

    # Query the HEALPix pixels covering the search box
    ids = hp.query_polygon(nside, vertices_uec, inclusive=True,fact=64)

    # Query the database to retrieve specific data columns for pixel numbers at a given level from a specified star catalog table. 
    K5_indices_list = index_query_sql(catalog_indices_db, tb_name, level, ids)

    # Load and concatenate data from relevant tile files
    dfs = _load_files(K5_indices_list,dir_sc,sc_name,_mode,max_num_per_tile)
    
    df = pd.concat(dfs,ignore_index=True)
    if df.empty: 
        return df,level,nside,ids,pixel_size,fov_min

    # Specific handling for different catalogs
    if sc_name in ['hyg37','at-hyg24']:
        df['ra'] = (df['ra'].astype(float)*15).round(8) # Convert hourangle to deg
        df['epoch'] = 2000.0
        columns_dict = {'pmra':'pm_ra', 'pmdec':'pm_dec'}
        df.rename(columns=columns_dict, inplace=True)
    elif sc_name == 'gsc30':
        columns_dict = {'rapm':'pm_ra', 'decpm':'pm_dec'}
        df.rename(columns=columns_dict, inplace=True)
    elif sc_name == 'gaiadr3':
        columns_dict = {'pmra': 'pm_ra', 'pmdec': 'pm_dec'}
        df.rename(columns=columns_dict, inplace=True)
    elif sc_name == 'ucac5':
        columns_dict = {'pmur':'pm_ra', 'pmud':'pm_dec','epu':'epoch'}
        df.rename(columns=columns_dict, inplace=True) 
    elif sc_name == 'usnob':
        columns_dict = {'pmRA':'pm_ra', 'pmDEC':'pm_dec','Epoch':'epoch'}
        df.rename(columns=columns_dict, inplace=True)
    elif sc_name == '2mass':
        df['epoch'] = Time(df['jdate'].astype(float), format='jd').jyear.round(3)

    # Ensure required columns are numeric
    required_columns = ['ra', 'dec', 'mag', 'epoch']
    pm_columns = ['pm_ra', 'pm_dec']

    if set(pm_columns).issubset(df.columns):
        required_columns.extend(pm_columns)

    df[required_columns] = df[required_columns].apply(pd.to_numeric, errors='coerce')

    # Filter by magnitude threshold
    mag_flag = (df['mag'] < mag_threshold)
    df = df[mag_flag].sort_values(by=['mag'])

    # Correct the proper motion
    dt = float(t_pm) - df['epoch']
    if {'pm_ra', 'pm_dec'}.issubset(df.columns):
        df['ra'] +=  df['pm_ra']/3.6e6 * dt
        df['dec'] += df['pm_dec']/3.6e6 * dt
        df['epoch'] = t_pm
    else:
        warnings.warn(f'Proper motion data for stars in catalog {sc_name} are not found.')

    # Filter stars within the search box
    ra_flag = np.abs(df['ra'] - (ra_min + ra_max)/2) < (ra_max - ra_min)/2
    dec_flag = np.abs(df['dec']- (dec_min + dec_max)/2) < (dec_max - dec_min)/2
    df = df[ra_flag & dec_flag]
    df.reset_index(drop=True,inplace=True)  

    return df,level,nside,ids,pixel_size,fov_min

def search_cone_raw(center, radius, dir_sc, sc_name, _mode, tb_name, catalog_indices_db, mag_threshold,t_pm,fov_min,max_num_per_tile=None):
    """
    This function performs a cone search of stars in specified star catalogs within a given RA/Dec center and radius.
    It applies magnitude filtering and proper motion correction based on the input parameters.

    Inputs:
        center -> [list] Center of the cone in form of [ra_c, dec_c] in degrees.
        radius -> [float] Angular radius of the cone.
        dir_sc -> [str] Directory where the star catalog files are stored.
        sc_name -> [str] Name of the star catalog.
        _mode -> [str] Mode that determines the number of header rows to skip; 'raw' mode skips two rows, other modes skip one row.
        tb_name -> [str] Name of the star catalog table.
        catalog_indices_db -> [str] Path to the database containing the star catalog indices.
        mag_threshold -> [float] Apparent magnitude limit.
        t_pm -> [float] Epoch to which the stars are unified.
        fov_min -> [float] Field of view parameters in degrees. It determines the hierarchical division of the sky region in HEALPix,
        ensuring that each tile's size is greater than one-quarter of the FOV and less than or equal to half of the FOV.
        max_num_per_tile -> [int, optional, default=None] Maximum number of stars to include in each tile, sorted by brightness.
        If None, includes all stars meeting the magnitude criteria.
    Outputs:
        df -> [pd.DataFrame] DataFrame containing the search results.
        level -> [str] HEALPix level used for the search.
        nside -> [int] NSIDE parameter of the HEALPix grid used for the search.
        ids -> [list] List of HEALPix pixel IDs covered by the search cone.
        pixel_size -> [float] Approximate size of each HEALPix pixel in degrees.
        fov_min -> [float] Field of view parameters in degrees.
    """
    # Extract center coordinates for the cone search
    ra_c, dec_c = center

    # Find the HEALPix parameters for a given field of view (FOV) in degrees
    if fov_min is None: fov_min = radius*2

    level, nside, npix, pixel_size = find_healpix_level(fov_min)

    # Compute the unit vector for the cone center
    uec = hp.ang2vec(ra_c, dec_c, lonlat=True)

    # Query the HEALPix pixels covering the cone
    ids = hp.query_disc(nside, uec, np.radians(radius), inclusive=True,fact=64)

    # Query the database to retrieve specific data columns for pixel numbers at a given level from a specified star catalog table. 
    K5_indices_list = index_query_sql(catalog_indices_db, tb_name, level, ids)

    # Load and concatenate data from relevant tile files
    dfs = _load_files(K5_indices_list,dir_sc,sc_name,_mode,max_num_per_tile)
    
    df = pd.concat(dfs,ignore_index=True)
    if df.empty: 
        return df,level,nside,ids,pixel_size,fov_min

    # Specific handling for different catalogs
    if sc_name in ['hyg37','at-hyg24']:
        df['ra'] = (df['ra'].astype(float)*15).round(8) # Convert hourangle to deg
        df['epoch'] = 2000.0
        columns_dict = {'pmra':'pm_ra', 'pmdec':'pm_dec'}
        df.rename(columns=columns_dict, inplace=True)
    elif sc_name == 'gsc30':
        columns_dict = {'rapm':'pm_ra', 'decpm':'pm_dec'}
        df.rename(columns=columns_dict, inplace=True)
    elif sc_name == 'gaiadr3':
        columns_dict = {'pmra': 'pm_ra', 'pmdec': 'pm_dec'}
        df.rename(columns=columns_dict, inplace=True)
    elif sc_name == 'ucac5':
        columns_dict = {'pmur':'pm_ra', 'pmud':'pm_dec','epu':'epoch'}
        df.rename(columns=columns_dict, inplace=True) 
    elif sc_name == 'usnob':
        columns_dict = {'pmRA':'pm_ra', 'pmDEC':'pm_dec','Epoch':'epoch'}
        df.rename(columns=columns_dict, inplace=True)
    elif sc_name == '2mass':
        df['epoch'] = Time(df['jdate'].astype(float), format='jd').jyear.round(3)

    # Ensure required columns are numeric
    required_columns = ['ra', 'dec', 'mag', 'epoch']
    pm_columns = ['pm_ra', 'pm_dec']

    if set(pm_columns).issubset(df.columns):
        required_columns.extend(pm_columns)

    df[required_columns] = df[required_columns].apply(pd.to_numeric, errors='coerce')

    # Filter by magnitude threshold
    mag_flag = (df['mag'] < mag_threshold)
    df = df[mag_flag].sort_values(by=['mag'])

    # Correct the proper motion
    dt = float(t_pm) - df['epoch']
    if {'pm_ra', 'pm_dec'}.issubset(df.columns):
        df['ra'] +=  df['pm_ra']/3.6e6 * dt
        df['dec'] += df['pm_dec']/3.6e6 * dt
        df['epoch'] = t_pm
    else:
        warnings.warn('Proper motion data for stars in catalog {:s} are not found.'.format(sc_name))

    # Calculate the angular distance between the cone center and each grid point.
    ra_c_rad,dec_c_rad = np.radians([ra_c,dec_c])
    ra_rad,dec_rad = np.radians(df[['ra', 'dec']].values).T

    angular_distance_cos = separation(ra_c_rad,dec_c_rad,ra_rad,dec_rad)

    # Determine whether each grid point is inside the cone.
    inside_cone = angular_distance_cos >= np.cos(np.deg2rad(radius))

    df = df[inside_cone]
    df.reset_index(drop=True,inplace=True)   

    return df,level,nside,ids,pixel_size,fov_min

def search_box_reduced(radec_box, dir_sc, sc_name, _mode, tb_name, catalog_indices_db, mag_threshold, t_pm, fov_min,max_num_per_tile=None):
    """
    This function performs a rectangular search of stars in specified reduced star catalogs within a given RA/Dec box.
    It applies magnitude filtering and proper motion correction based on the input parameters.

    Inputs:
        radec_box -> [list] Rectangular search area in form of [ra_min, dec_min, ra_max, dec_max] in degrees.
        dir_sc -> [str] Path of star catalog tile files.
        sc_name -> [str] Name of the star catalog.
        _mode -> [str] Mode that determines the number of header rows to skip; 'raw' mode skips two rows, other modes skip one row.
        tb_name -> [str] Name of the star catalog table.
        catalog_indices_db -> [str] Path to the database containing the star catalog indices.
        mag_threshold -> [float] Apparent magnitude limit.
        t_pm -> [float] Epoch to which the stars are unified.
        fov-min -> [float] Field of view parameters in degrees. It determines the hierarchical division of the sky region in HEALPix,
        ensuring that each tile's size is greater than one-quarter of the FOV and less than or equal to half of the FOV.
        max_num_per_tile -> [int, optional, default=None] Maximum number of stars to include in each tile, sorted by brightness. If None, includes all stars meeting the magnitude criteria.
    Outputs:
        df -> [pd.DataFrame] DataFrame containing stars within the search area.
        level -> [str] HEALPix level used for the search.
        nside -> [int] NSIDE parameter of the HEALPix grid used for the search.
        ids -> [list] List of HEALPix pixel IDs covered by the search box.
        pixel_size -> [float] Approximate size of each HEALPix pixel in degrees.
        fov_min -> [float] Field of view parameters in degrees.
    """
    # Extract the boundary from the search box
    ra_min, dec_min, ra_max, dec_max = radec_box
    dec_c,ra_c = np.mean([dec_max,dec_min]),np.mean([ra_max,ra_min])
    ra_range,dec_range = ra_max-ra_min,dec_max-dec_min
    if fov_min is None: fov_min = min(ra_range*np.cos(np.radians(dec_c)),dec_range)

    vertices = np.array([
        [ra_min, dec_min],  # Bottom left
        [ra_max, dec_min],  # Bottom right
        [ra_max, dec_max],  # Top right
        [ra_min, dec_max]   # Top left
    ])

    vertices_uec = hp.ang2vec(vertices[:,0],vertices[:,1],lonlat=True)

    # Find the HEALPix parameters for a given field of view (FOV) in degrees
    level, nside, npix, pixel_size = find_healpix_level(fov_min)

    # Query the HEALPix pixels covering the search box
    ids = hp.query_polygon(nside, vertices_uec, inclusive=True,fact=64)

    # Query the database to retrieve specific data columns for pixel numbers at a given level from a specified star catalog table. 
    K5_indices_list = index_query_sql(catalog_indices_db, tb_name, level, ids)

    # Load and concatenate data from relevant tile files
    dfs = _load_files(K5_indices_list,dir_sc,sc_name,_mode,max_num_per_tile)
    
    df = pd.concat(dfs,ignore_index=True)
    if df.empty: 
        return df,level,nside,ids,pixel_size,fov_min

    # Ensure required columns are numeric
    required_columns = ['ra', 'dec', 'mag', 'epoch']
    pm_columns = ['pm_ra', 'pm_dec']

    if set(pm_columns).issubset(df.columns):
        required_columns.extend(pm_columns)

    df[required_columns] = df[required_columns].apply(pd.to_numeric, errors='coerce')

    # Filter stars based on the magnitude threshold
    mag_flag = df['mag'] < mag_threshold
    df = df[mag_flag].sort_values(by=['mag'])

    # Calculate the time difference from the epoch
    dt = t_pm - df['epoch']

    # Apply proper motion correction if data is available
    if {'pm_ra', 'pm_dec'}.issubset(df.columns):
        df['ra'] += df['pm_ra'] / 3.6e6 * dt
        df['dec'] += df['pm_dec'] / 3.6e6 * dt
        # Update the epoch to the search epoch
        df['epoch'] = t_pm
    else:
        warnings.warn(f'Proper motion data for stars in catalog {sc_name} are not found.')

    # Filter stars within the search area
    ra_flag = np.abs(df['ra'] - (ra_min + ra_max) / 2) < (ra_max - ra_min) / 2
    dec_flag = np.abs(df['dec'] - (dec_min + dec_max) / 2) < (dec_max - dec_min) / 2
    df = df[ra_flag & dec_flag]

    df.reset_index(drop=True,inplace=True)
    
    return df,level,nside,ids,pixel_size,fov_min


def search_cone_reduced(center, radius, dir_sc, sc_name, _mode, tb_name, catalog_indices_db, mag_threshold,t_pm,fov_min,max_num_per_tile=None):
    """
    This function performs a conical search of stars in specified reduced star catalogs within a given RA/Dec center and radius.
    It applies magnitude filtering and proper motion correction based on the input parameters.

    Inputs:
        center -> [list] Center of the cap in form of [Ra, Dec] in degrees.
        radius -> [float] Angular radius of the cap in degrees.
        dir_sc -> [str] Path of star catalog tile files.
        sc_name -> [str] Name of the star catalog.
        _mode -> [str] Mode that determines the number of header rows to skip; 'raw' mode skips two rows, other modes skip one row.
        tb_name -> [str] Name of the star catalog table.
        catalog_indices_db -> [str] Path to the database containing the star catalog indices.
        mag_threshold -> [float] Apparent magnitude limit.
        t_pm -> [float] Epoch to which the stars are unified.
        fov_min -> [float] Field of view parameters in degrees. It determines the hierarchical division of the sky region in HEALPix,
        ensuring that each tile's size is greater than one-quarter of the FOV and less than or equal to half of the FOV.
        max_num_per_tile -> [int, optional, default=None] Maximum number of stars to include in each tile, sorted by brightness. If None, includes all stars meeting the magnitude criteria.
    Outputs:
        df -> [pd.DataFrame] DataFrame containing stars within the cone search area.
        level -> [str] HEALPix level used for the search.
        nside -> [int] NSIDE parameter of the HEALPix grid used for the search.
        ids -> [list] List of HEALPix pixel IDs covered by the search cone.
        pixel_size -> [float] Approximate size of each HEALPix pixel in degrees.
        fov_min -> [float] Field of view parameters in degrees.
    """
    # Extract center coordinates for the cone search
    ra_c, dec_c = center

    # Find the HEALPix parameters for a given field of view (FOV) in degrees
    if fov_min is None: fov_min = radius*2

    level, nside, npix, pixel_size = find_healpix_level(fov_min)

    # Compute the unit vector for the cone center
    uec = hp.ang2vec(ra_c, dec_c, lonlat=True)

    # Query the HEALPix pixels covering the cone
    ids = hp.query_disc(nside, uec, np.radians(radius), inclusive=True,fact=64)

    # Query the database to retrieve specific data columns for pixel numbers at a given level from a specified star catalog table. 
    K5_indices_list = index_query_sql(catalog_indices_db, tb_name, level, ids)
 
    # Load and concatenate data from relevant tile files
    dfs = _load_files(K5_indices_list,dir_sc,sc_name,_mode,max_num_per_tile)
    
    df = pd.concat(dfs,ignore_index=True)
    if df.empty: 
        return df,level,nside,ids,pixel_size,fov_min

    # Ensure required columns are numeric
    required_columns = ['ra', 'dec', 'mag', 'epoch']
    pm_columns = ['pm_ra', 'pm_dec']

    if set(pm_columns).issubset(df.columns):
        required_columns.extend(pm_columns)

    df[required_columns] = df[required_columns].apply(pd.to_numeric, errors='coerce')

    # Filter stars based on magnitude threshold
    df = df[df['mag'] < mag_threshold].sort_values(by='mag')

    # Proper motion adjustment, if available
    dt = t_pm - df['epoch']
    if {'pm_ra', 'pm_dec'}.issubset(df.columns):
        df['ra'] += df['pm_ra'] / 3.6e6 * dt
        df['dec'] += df['pm_dec'] / 3.6e6 * dt
        # Update the epoch for the search results
        df['epoch'] = t_pm
    else:
        warnings.warn(f'Proper motion data for stars in catalog {sc_name} not found.')

    # Calculate the angular distance between the cap center and each grid point.
    ra_c_rad,dec_c_rad = np.radians([ra_c,dec_c])
    ra_rad,dec_rad = np.radians(df[['ra', 'dec']].values).T

    angular_distance_cos = separation(ra_c_rad,dec_c_rad,ra_rad,dec_rad)

    # Determine whether each grid point is inside the cone.
    inside_cone = angular_distance_cos >= np.cos(np.deg2rad(radius))

    df = df[inside_cone]
    df.reset_index(drop=True,inplace=True)

    return df,level,nside,ids,pixel_size,fov_min

def search_box_simplified(radec_box, dir_sc, sc_name, _mode, tb_name, catalog_indices_db, fov_min, max_num_per_tile=None,astrometry_corrections={}):
    """
    This function performs a rectangular search of stars in specified simplified star catalogs within a given RA/Dec box.
    It applies various astrometry corrections based on the input parameters.

    Inputs:
        radec_box -> [array-like] Rectangular search area in form of [ra_min, dec_min, ra_max, dec_max] in degrees.
        dir_sc -> [str] Path of star catalog tile files.
        sc_name -> [str] Name of the star catalog.
        _mode -> [str] Mode that determines the number of header rows to skip; 'raw' mode skips two rows, other modes skip one row.
        tb_name -> [str] Name of the star catalog table.
        catalog_indices_db -> [str] Path to the database containing the star catalog indices.
        fov_min -> [float] Field of view parameters in degrees. It determines the hierarchical division of the sky region in HEALPix,
        ensuring that each tile's size is greater than one-quarter of the FOV and less than or equal to half of the FOV.
        max_num_per_tile -> [int, optional, default=None] Maximum number of stars to include in each tile, sorted by brightness. If None, includes all stars meeting the criteria.
        astrometry_corrections -> [dict, optional, default={}] Dictionary specifying the types of astrometry corrections to apply.
            - 't' -> [str] Observation time in UTC, such as '2019-02-26T20:11:14.347'.
            - 'proper-motion' -> [None] If present, apply proper motion correction.
            - 'aberration' -> [tuple] Aberration correction parameters. Observer's velocity relative to Earth's center (vx, vy, vz) in km/s.
            - 'parallax' -> [None] If present, apply parallax correction.
            - 'deflection' -> [None] If present, apply light deflection correction.
    Outputs:
        df -> [pd.DataFrame] DataFrame containing stars within the rectangular search area.
        level -> [str] HEALPix level used for the search.
        nside -> [int] NSIDE parameter of the HEALPix grid used for the search.
        ids -> [list] List of HEALPix pixel IDs covered by the search box.
        pixel_size -> [float] Approximate size of each HEALPix pixel in degrees.
        fov_min -> [float] Field of view parameters in degrees.
    """
    ra_min, dec_min, ra_max, dec_max = radec_box
    dec_c,ra_c = np.mean([dec_max,dec_min]),np.mean([ra_max,ra_min])
    ra_range,dec_range = ra_max-ra_min,dec_max-dec_min
    if fov_min is None: fov_min = min(ra_range*np.cos(np.radians(dec_c)),dec_range)

    vertices = np.array([
        [ra_min, dec_min],  # Bottom left
        [ra_max, dec_min],  # Bottom right
        [ra_max, dec_max],  # Top right
        [ra_min, dec_max]   # Top left
    ])

    vertices_uec = hp.ang2vec(vertices[:,0],vertices[:,1],lonlat=True)

    # Find the HEALPix parameters for a given field of view (FOV) in degrees
    level, nside, npix, pixel_size = find_healpix_level(fov_min)

    # Query the HEALPix pixels covered by the search box
    ids = hp.query_polygon(nside, vertices_uec, inclusive=True,fact=64)

    # Query the database to retrieve specific data columns for pixel numbers at a given level from a specified star catalog table. 
    K5_indices_list = index_query_sql(catalog_indices_db, tb_name, level, ids)

    # Load and concatenate data from relevant tile files
    dfs = _load_files(K5_indices_list,dir_sc,sc_name,_mode,max_num_per_tile)
    
    df = pd.concat(dfs,ignore_index=True)
    if df.empty: 
        return df,level,nside,ids,pixel_size,fov_min

    # Ensure required columns are numeric
    required_columns = ['ra', 'dec', 'mag', 'epoch']
    pm_columns = ['pm_ra', 'pm_dec']
    dist_columns = ['dist']

    if set(pm_columns).issubset(df.columns):
        required_columns.extend(pm_columns)
    if set(dist_columns).issubset(df.columns):
        required_columns.extend(dist_columns)

    df[required_columns] = df[required_columns].apply(pd.to_numeric, errors='coerce')

    # Filter stars within the search box
    ra_flag = np.abs(df['ra'] - ra_c) < ra_range / 2
    dec_flag = np.abs(df['dec'] - dec_c) < dec_range / 2
    df = df[ra_flag & dec_flag]

    ra_rad, dec_rad = np.radians(df[['ra', 'dec']].values).T

    if astrometry_corrections:
        df = apply_astrometry_corrections(df, astrometry_corrections, ra_rad, dec_rad)

    df = df.sort_values(by='mag')
    df.reset_index(drop=True, inplace=True)

    return df,level,nside,ids,pixel_size,fov_min

def search_cone_simplified(center, radius, dir_sc, sc_name, _mode, tb_name, catalog_indices_db, fov_min, max_num_per_tile=None,astrometry_corrections={}):
    """
    This function performs a cone search of stars in specified simplified star catalogs within a given RA/Dec center and radius.
    It applies various astrometry corrections based on the input parameters.

    Inputs:
        center -> [list] Center of the cap in form of [Ra, Dec] in degrees.
        radius -> [float] Angular radius of the cap in degrees.
        dir_sc -> [str] Directory of the star catalog tile files.
        sc_name -> [str] Name of the star catalog.
        _mode -> [str] Mode that determines the number of header rows to skip; 'raw' mode skips two rows, other modes skip one row.
        tb_name -> [str] Name of the star catalog table.
        catalog_indices_db -> [str] Path to the database containing the star catalog indices.
        fov_min -> [float] Field of view parameters in degrees. It determines the hierarchical division of the sky region in HEALPix,
        ensuring that each tile's size is greater than one-quarter of the FOV and less than or equal to half of the FOV.
        max_num_per_tile -> [int, optional, default=None] Maximum number of stars to include in each tile, sorted by brightness. If None, includes all stars meeting the criteria.
        astrometry_corrections -> [dict, optional, default={}] Dictionary specifying the types of astrometry corrections to apply.
            - 't' -> [str] Observation time in UTC, such as '2019-02-26T20:11:14.347'.
            - 'proper-motion' -> [None] If present, apply proper motion correction.
            - 'aberration' -> [tuple] Aberration correction parameters. Observer's velocity relative to Earth's center (vx, vy, vz) in km/s.
            - 'parallax' -> [None] If present, apply parallax correction.
            - 'deflection' -> [None] If present, apply light deflection correction.
    Outputs:
        df -> [pd.DataFrame] DataFrame containing stars within the cone search area.
        level -> [str] HEALPix level used for the search.
        nside -> [int] NSIDE parameter of the HEALPix grid used for the search.
        ids -> [list] List of HEALPix pixel IDs covered by the search cone.
        pixel_size -> [float] Approximate size of each HEALPix pixel in degrees.
        fov_min -> [float] Field of view parameters in degrees.
    """
    # Extract center coordinates for the cone search
    ra_c, dec_c = center

    # Find the HEALPix parameters for a given field of view (FOV) in degrees
    if fov_min is None: fov_min = radius*2

    level, nside, npix, pixel_size = find_healpix_level(fov_min)

    # Compute the unit vector for the cone center
    uec = hp.ang2vec(ra_c, dec_c, lonlat=True)

    # Query the HEALPix pixels covering the cone
    ids = hp.query_disc(nside, uec, np.radians(radius), inclusive=True,fact=64)

    # Query the database to retrieve specific data columns for pixel numbers at a given level from a specified star catalog table. 
    K5_indices_list = index_query_sql(catalog_indices_db, tb_name, level, ids)

    # Load and concatenate data from relevant tile files
    dfs = _load_files(K5_indices_list,dir_sc,sc_name,_mode,max_num_per_tile)
    
    df = pd.concat(dfs,ignore_index=True)

    if df.empty: 
        return df,level,nside,ids,pixel_size,fov_min

    # Ensure required columns are numeric
    required_columns = ['ra', 'dec', 'mag', 'epoch']
    pm_columns = ['pm_ra', 'pm_dec']
    dist_columns = ['dist']

    if set(pm_columns).issubset(df.columns):
        required_columns.extend(pm_columns)
    if set(dist_columns).issubset(df.columns):
        required_columns.extend(dist_columns)

    df[required_columns] = df[required_columns].apply(pd.to_numeric, errors='coerce')

    # Calculate the angular distance between the cap center and each grid point.
    ra_c_rad,dec_c_rad = np.radians([ra_c,dec_c])
    ra_rad,dec_rad = np.radians(df[['ra', 'dec']].values).T
    angular_distance_cos = separation(ra_c_rad,dec_c_rad,ra_rad,dec_rad)
    # Determine whether each grid point is inside the cap.
    inside_cone = angular_distance_cos >= np.cos(np.deg2rad(radius))
    df = df[inside_cone]
    ra_rad, dec_rad = ra_rad[inside_cone], dec_rad[inside_cone]

    if astrometry_corrections:
        df = apply_astrometry_corrections(df, astrometry_corrections, ra_rad, dec_rad)

    df = df.sort_values(by='mag')
    df.reset_index(drop=True, inplace=True)

    return df,level,nside,ids,pixel_size,fov_min