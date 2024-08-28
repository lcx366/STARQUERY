import os
import bisect
import pandas as pd
import numpy as np
import healpy as hp
from glob import glob
from colorama import Fore
from natsort import natsorted
from sqlalchemy import create_engine, MetaData, Table, Column, Integer, Index
from sqlalchemy.engine.reflection import Inspector
import h5py

from .wcs import xy_catalog
from .invariantfeatures import unique_triangles,unique_quads

K_RANGE = range(1, 12)  # Levels for star catalog indexing, from 1 to 11. Higher levels mean finer sky partitions.
K_RANGE_INVERSE = K_RANGE[::-1]
N_STARS = 6 # Number of stars to extract for each partition.

NSIDE = np.power(2,K_RANGE_INVERSE)  # Calculate nside for each K
NPIX = hp.nside2npix(NSIDE)  # Calculate the total number of pixels for each nside
PIXEL_SIZE = np.degrees(np.sqrt(4 * np.pi / NPIX))  # Calculate the pixel size in degrees, using the area of a pixel

# Define reasonable limits for FOV
LOWER_FOV = PIXEL_SIZE[0] * 2  # lower limits of reasonable field of view in degrees
UPPER_FOV = PIXEL_SIZE[-1] * 2  # upper limits of reasonable field of view in degrees

def find_healpix_level(fov_min):
    """
    Find the HEALPix parameters for a given field of view (FOV) in degrees
    such that each pixel's size is greater than one-quarter of the FOV and less than or equal to half of the FOV.

    Usage:
        >>> level, nside, npix, pixel_size = find_healpix_level(5.0)
    Inputs:
        fov_min -> [float] The field of view parameters in degrees from which the HEALPix parameters need to be determined.
    Outputs:
        level -> [str] Partition level, such as 'K4'.
        nside -> [int] The number of divisions along the side of a HEALPix base pixel.
        npix -> [int] The total number of pixels.
        pixel_size -> [float] Approximate pixel size in degrees.
    """
    # Check if the input FOV is within the reasonable range
    if fov_min <= LOWER_FOV or fov_min > UPPER_FOV:
        raise ValueError(f"The FOV must be between {LOWER_FOV} and {UPPER_FOV} degrees. Given FOV: {fov_min}")

    # Use bisect to find the index where pixel size is just above the minimum diameter
    index = bisect.bisect(PIXEL_SIZE, fov_min/2) - 1
    level = f'K{K_RANGE_INVERSE[index]}'

    return level, NSIDE[index], NPIX[index],  PIXEL_SIZE[index]

def calculate_star_indices(tile_path, nsides):
    """
    Calculate the indices for each star in a single tile file across different levels.

    Usage:
        >>> tile_path = 'starcatalogs/simplified/at-hyg24/at-hyg24-100.csv'
        >>> nsides = [2**i for i in range(4, 12)]  # From 2^4 to 2^11
        >>> indices = calculate_star_indices(tile_path, nsides)
    Inputs:
        tile_path -> [str] Path to the tile file.
        nsides -> [list] List of nsides, defining the resolution levels for Healpix tiling.
    Outputs:
        indices -> [numpy.ndarray] 2D array, each row represents a star, each column represents a level.
    """
    df = pd.read_csv(tile_path)
    n = len(df)
    # Initialize the index array
    indices = np.empty((n, len(nsides)+1), dtype=int)

    if n > 0:
        ra, dec = df['ra'], df['dec']
        # Compute the pixel number of all stars for each level
        for i, nside in enumerate(nsides):
            pixels = hp.ang2pix(nside, ra, dec, lonlat=True)
            indices[:, i] = pixels

        # Fill the last column with sequential IDs
        indices[:, -1] = np.arange(n)  
    
    return indices

def build_catalog_indices(dir_from,sc_name,tb_name):
    """
    This function aggregates indices from all tile files for each star across different levels, 
    generating a comprehensive csv-formatted catalog index file.

    Usage:
        >>> dir_from = 'starcatalogs/simplified/at-hyg24/'
        >>> sc_name,tb_name = 'at-hyg24','at-hyg24_mag12.0_epoch2019.5'
        >>> indices_path = build_catalog_indices(dir_from,sc_name,tb_name)
    Inputs:
        dir_from -> [str] Directory containing the tile files.
        sc_name -> [str] Star catalog name.
        tb_name -> [str] Table name for the indices.
    Outputs:
        indices_path -> [str] Path to the CSV file containing the catalog indices.
    """
    nsides = [2**i for i in K_RANGE]  # From 2^1 to 2^11

    # Construct file pattern
    file_pattern = os.path.join(dir_from, f'{sc_name}-*.csv')
    files = natsorted(glob(file_pattern))
    
    # Initialize index list
    indices = []
    
    # Iterate over all files
    for i, file in enumerate(files):
        print(f'Building star indices in healpix {Fore.BLUE}{i+1}{Fore.RESET} of {len(files)}', end='\r')

        # Calculate the indices for each star in a single tile file across different nsides.
        tile_indices = calculate_star_indices(file, nsides)
        tile_indices[:,4] = i
        indices.append(tile_indices)
    
    # Concatenate all results into a large array
    indices = np.vstack(indices)

    # Save to CSV file
    column_names = [f'K{i}' for i in K_RANGE] + ['K5_SUB']
    df = pd.DataFrame(indices, columns=column_names)
    df = df.sort_values(by=list(df.columns)).reset_index(drop=True)

    dir_indices = dir_from.split('starcatalogs')[0] + 'starcatalogs/indices/'   
    os.makedirs(dir_indices, exist_ok=True) # Ensure directory exists  
    indices_path = os.path.join(dir_indices, f'{tb_name}.csv')   

    df.to_csv(indices_path, index=False)
    return indices_path

def generate_catalog_db(catalog_indices_db,catalog_indices_csv):
    """
    Build and manage an SQL database, including adding, indexing data tables of star catalogs.

    Usage:
        >>> catalog_indices_db = 'catalogs.db'
        >>> catalog_indices_csv = 'starcatalogs/indices/at-hyg24_mag12.0_epoch2019.5.csv'
        >>> generate_catalog_db(catalog_indices_db,catalog_indices_csv)
    Inputs:
        catalog_indices_db -> [str] Path to the database file.
        catalog_indices_csv -> [str] Path to the star catalog indices file (CSV format).
    Outputs:
        Message indicating the result of the database operation.
    """
    engine = create_engine(f'sqlite:///{catalog_indices_db}')
    metadata = MetaData()

    # Define the table name based on the CSV filename
    tb_name = catalog_indices_csv.split('/')[-1].split('.csv')[0]

    df = pd.read_csv(catalog_indices_csv)
    df.to_sql(tb_name, engine, if_exists='replace', index=False)
    print(f"Table '{tb_name}' has been added to the database.")

    # Define the table with metadata and columns
    table = Table(
        tb_name, metadata,
        Column('K5_SUB', Integer, primary_key=True),
        *(Column(f'K{i}', Integer) for i in K_RANGE)
    )

    metadata.create_all(engine)  # Creates table with columns defined
    for i in K_RANGE:
        index = Index(f'idx_{tb_name}_K{i}', table.columns[f'K{i}'])
        index.create(engine)  # Create the index on the engine
        print(f"Index on {Fore.BLUE}'K{i}'{Fore.RESET} has been created.", end='\r')

def delete_table(catalog_indices_db, tb_name):
    """
    Manage an SQL database by deleting data tables of star catalogs.

    Usage:
        >>> catalog_indices_db,tb_name = 'catalogs.db','at-hyg24_mag12.0_epoch2019.5'
        >>> delete_table(catalog_indices_db,tb_name)
    Inputs:
        catalog_indices_db -> [str] Path to the database file.
        tb_name -> [str] Name of the table to delete.
    Outputs:
        Message indicating the result of the database operation.
    """

    # Create a database connection using SQLAlchemy
    engine = create_engine(f'sqlite:///{catalog_indices_db}')
    metadata = MetaData()
    metadata.bind = engine

    # Create an Inspector object
    inspector = Inspector.from_engine(engine)

    # Check if the table exists and then delete
    if inspector.has_table(tb_name):
        table = Table(tb_name, metadata, autoload_with=engine)
        table.drop(engine)
        print(f"Table '{tb_name}' has been deleted from the database.")
    else:
        print(f"Table '{tb_name}' does not exist in the database.")

def index_query_sql(catalog_indices_db, tb_name, level, ids):
    """
    Queries the database to retrieve specific data columns for pixel numbers at a given level
    from a specified star catalog table. This function returns a list of lists, where each inner list
    contains tuples of 'K5' and a list of 'K5_SUB' values, grouped by 'K5'.

    Usage:
        >>> catalog_indices_db = 'starcatalogs/indices/catalogs.db'
        >>> tb_name = 'at-hyg24_mag12.0_epoch2019.5'
        >>> level = 'K6'
        >>> ids = [49150, 49149, 49148]
        >>> K5_indices_list = index_query_sql(catalog_indices_db, tb_name, level, ids)
    Inputs:
        catalog_indices_db -> [str] Path to the star catalog index database containing HEALPix tiling data at multiple resolutions.
        tb_name -> [str] Name of the star catalog table, e.g., 'hyg37', 'hyg37_mag9.0_epoch2022.0'.
        level -> [str] The resolution level from which to fetch the pixel numbers, such as 'K1', 'K2', ..., 'K11'
        ids -> [list of int] List of pixel IDs at the specified level, e.g., [49150, 49149].
    Outputs:
        K5_indices_list -> [list] Each inner list contains tuples, each tuple includes a 'K5' value and a list of 'K5_SUB' values.
    """
    # Create a database connection using SQLAlchemy
    engine = create_engine(f'sqlite:///{catalog_indices_db}')

    # Convert the list of IDs into a string suitable for SQL queries
    id_list_str = ', '.join(map(str, ids))

    # Check if level is 'K5' to avoid duplicate columns
    if level == 'K5':
        query = f"SELECT K5, K5_SUB FROM '{tb_name}' WHERE K5 IN ({id_list_str})"
    else:
        query = f"SELECT {level}, K5, K5_SUB FROM '{tb_name}' WHERE {level} IN ({id_list_str})"

    filtered_df = pd.read_sql_query(query, engine)

    if filtered_df.empty:
        return []

    # Check if level is 'K5'
    if level == 'K5':
        grouped = filtered_df.groupby('K5')['K5_SUB'].apply(list).reset_index()
        K5_indices_list = grouped.apply(lambda row: [(row['K5'], row['K5_SUB'])], axis=1).tolist()

    else:
        # Use pivot_table to restructure data, first group by level, then by K5, and finally convert K5_SUB to list
        grouped = filtered_df.groupby([level, 'K5'])['K5_SUB'].apply(list).reset_index()
        pivoted = grouped.pivot(index=level, columns='K5', values='K5_SUB')
        K5_indices_list = pivoted.apply(lambda x: list(zip(x.dropna().index, x.dropna())), axis=1).tolist()
        
    return K5_indices_list

def fetch_partitions(db_path, tb_name, level, n_stars, dir_sc, sc_name):
    """
    Fetch and process star catalog data partitions from a database, grouping them by HEALPix levels.

    Usage:
        >>> db_path = 'starcatalogs/indices/catalogs.db'
        >>> tb_name = 'at-hyg24_mag12.0_epoch2019.5'
        >>> level = 'K6'
        >>> n_stars = 8
        >>> dir_sc = 'starcatalogs/simplified/at-hyg24/'
        >>> sc_name = 'at-hyg24'
        >>> results_dict = fetch_partitions(db_path, tb_name, level, n_stars, dir_sc, sc_name)
    Inputs:
        db_path -> [str] Path to the database file.
        tb_name -> [str] Name of the star catalog table.
        level -> [str] HEALPix level for partitioning.
        n_stars -> [int] Number of stars to extract for each partition.
        dir_sc -> [str] Directory containing the star catalog files.
        sc_name -> [str] Star catalog name.
    Outputs:
        results_dict -> [dict] Dictionary with K1 partitions as keys and 2D numpy arrays of RA/Dec values as values.
    """
    engine = create_engine(f'sqlite:///{db_path}')

    # Build SQL query
    columns = '"K1", "K5", "K5_SUB"' + (f', "{level}"' if level not in ['K1', 'K5'] else '')
    query = f"SELECT {columns} FROM \"{tb_name}\""
    full_df = pd.read_sql_query(query, engine)

    # Preload necessary star catalog data
    needed_files = full_df['K5'].unique()
    star_data = {}
    for k5 in needed_files:
        file_path = os.path.join(dir_sc, f"{sc_name}-{int(k5)}.csv")
        if os.path.exists(file_path):
            star_data[k5] = pd.read_csv(file_path)

    results_dict = {}

    # Process levels K6 and above
    if int(level[1:]) > 5:
        for (k1, k5, lvl), group in full_df.groupby(['K1', 'K5', level]):
    
            # Check if all indices are within bounds
            if group['K5_SUB'].max() >= len(star_data[k5]):
                print(f"Index out of bounds for K5: {k5} with max index {group['K5_SUB'].max()}")
                continue

            relevant_data = star_data[k5].iloc[group['K5_SUB']].copy()
            sorted_data = relevant_data.sort_values(by='mag').head(n_stars)
            results_dict.setdefault(k1, []).append(sorted_data[['ra', 'dec']].to_numpy())
    else:
        # Process levels K5 and below
        for (k1, lvl), group in full_df.groupby(['K1', level]):
            k5_list = group['K5'].unique()
            combined_df = pd.concat([star_data[k5] for k5 in k5_list])
            sorted_data = combined_df.sort_values(by='mag').head(n_stars)
            results_dict.setdefault(k1, []).append(sorted_data[['ra', 'dec']].to_numpy())

    # Convert list in K1 to a single 2D numpy array
    for k1 in results_dict:
        results_dict[k1] = np.unique(np.vstack(results_dict[k1]), axis=0)

    return results_dict

def sort_data_by_dec(data):
    """
    Sorts data by absolute declination, grouping by declination and sorting by right ascension within each group.

    Usage:
        >>> data = np.array([[15.0, 30.0], [45.0, -30.0], [60.0, 15.0], [75.0, -15.0]])
        >>> sorted_data, original_indices = sort_data_by_dec(data)
    Inputs:
        data -> [numpy.ndarray] Array of shape (n, 2) where data[:, 0] is right ascension (RA) and data[:, 1] is declination (Dec).
    Outputs:
        sorted_data -> [numpy.ndarray] Sorted array where data is sorted by absolute declination groups, and by RA within those groups.
        original_indices -> [numpy.ndarray] Original indices of the sorted data.
    """
    
    # Convert the input array to a pandas DataFrame for easier handling
    df = pd.DataFrame(data, columns=['RA', 'Dec'])
    df['Original_Index'] = df.index  # Save the original index for tracking
    
    # Sort the DataFrame by the absolute values of declination
    sorted_df = df.sort_values(by='Dec', key=lambda x: np.abs(x)).reset_index()
    
    # Group the sorted DataFrame by the absolute values of declination
    grouped = sorted_df.groupby(np.abs(sorted_df['Dec']))
    
    # Prepare a list to store the sorted output
    output = []
    
    # Iterate through each group, sorting by RA and keeping north and south separate
    for name, group in grouped:
        # Sort the north hemisphere data by RA
        north = group[group['Dec'] >= 0].sort_values(by='RA')
        # Sort the south hemisphere data by RA
        south = group[group['Dec'] < 0].sort_values(by='RA')
        
        # Append sorted groups to the output list
        output.extend(north[['RA', 'Dec', 'Original_Index']].values.tolist())
        output.extend(south[['RA', 'Dec', 'Original_Index']].values.tolist())

    output_np = np.array(output)
    
    return output_np[:,:2],output_np[:,2].astype(int)

def h5_hashes(db_path, tb_name, dir_sc, sc_name, k_min,k_max,mode_invariants,pixel_width=0.001, theta=0):
    """
    For each region of the sky, this function calculates geometric invariants of star configurations, such as triangle
    edge length ratios or quadrilateral invariants (based on methods from astrometry.net). These invariants are useful
    for tasks like star pattern recognition and matching.

    Usage:
        >>> db_path = 'starcatalogs/indices/catalogs.db'
        >>> tb_name = 'at-hyg24_mag12.0_epoch2019.5'
        >>> dir_sc = 'starcatalogs/simplified/at-hyg24/'
        >>> sc_name = 'at-hyg24'
        >>> k_min,k_max = 1,6
        >>> mode_invariants = 'triangles'
        >>> h5_file = h5_hashes(db_path, tb_name, dir_sc, sc_name,k_min,k_max, mode_invariants)
    Inputs:
        db_path -> [str] Path to the database file.
        tb_name -> [str] Name of the star catalog table.
        dir_sc -> [str] Directory containing the star catalog files.
        sc_name -> [str] Star catalog name.
        k_min -> [int] Minimum HEALPix hierarchy level.
        k_max -> [int] Maximum HEALPix hierarchy level.
        mode_invariants -> [str] Type of invariants to calculate ('triangles' or 'quads').
        pixel_width -> [float, optional, default=0.001] Pixel width in degrees.
        theta -> [float, optional, default=0] Rotation angle in degrees.
    Outputs:
        outh5 -> [str] Path to the generated HDF5 file.
    """
    outh5 = db_path.replace('catalogs.db', f"{tb_name}_{mode_invariants}_K{k_min}_K{k_max}.h5")

    # Check if file already exists
    if not os.path.exists(outh5):
        # Calculate the center of healpix polygon
        fp_radecs = hp.pix2ang(2, range(48), lonlat=True)
        fp_radecs = np.stack(fp_radecs).T

        # Sort fp_radecs by declination
        sorted_fp_radecs,sorted_indices = sort_data_by_dec(fp_radecs)

        with h5py.File(outh5, 'w') as f:
            for i in range(k_min, k_max + 1):  # Iterate through the levels in the range
                lvl = f"K{i}"
                print(f'Generating hdf5 data in level {Fore.BLUE}{lvl}{Fore.RESET} of {"K1 -> K11"}', end='\r')

                data_dict = fetch_partitions(db_path, tb_name, lvl, N_STARS, dir_sc, sc_name)

                for j,idx in enumerate(sorted_indices):
                    radec = data_dict[idx]
                    x, y, wcs = xy_catalog(fp_radecs[idx], radec, pixel_width, theta)
                    xy = np.stack([x, y]).T

                    if mode_invariants == 'triangles':
                        invariants, asterisms = unique_triangles(xy)
                    elif mode_invariants == 'quads':
                        invariants, asterisms = unique_quads(xy)
                    else:
                        raise ValueError(f'Unrecognized mode invariants type: {mode_invariants}')

                    grp = f.create_group(f"{lvl}/{j}")
                    grp.create_dataset('xy', data=xy)
                    grp.create_dataset('invariants', data=invariants)
                    grp.create_dataset('asterisms', data=asterisms)

            f.create_dataset("fp_radecs", data=sorted_fp_radecs)

    return outh5

def read_h5_hashes(infile):
    """
    Parse the HDF5 file containing geometric invariants of star configurations for each level and partition.

    Usage:
        >>> infile = 'starcatalogs/indices/at-hyg24_mag13.0_epoch2019.5_triangles_K1_K6.h5'
        >>> data = read_h5_hashes(infile)
    Inputs:
        infile -> [str] Path to the HDF5 file.
    Outputs:
        data -> [dict] A dictionary where keys are levels ('K1' to 'K11') and 'fp_radecs'. Values for levels
                      are lists of tuples containing data for each partition ('xy', 'invariants', 'asterisms'),
                      and value for 'fp_radecs' is the corresponding central pointing.
    """
    data = {}
    with h5py.File(infile, 'r') as file:
        # Read the fixed point radec positions
        fp_radecs = data['fp_radecs'] = file['fp_radecs'][:]
        keys = np.arange(len(fp_radecs)).astype(str)

        # Iterate over each level in the file
        for lvl in list(file.keys())[:-1]:  # Assuming levels K1 to K11 are in the file
            level_data = []
            if lvl in file: # Iterate over each partition in the level
                for k1 in keys:
                    partition = file[lvl][k1]
                    xy = partition['xy'][:]
                    invariants = partition['invariants'][:]
                    asterisms = partition['asterisms'][:]
                    level_data.append((xy, invariants, asterisms))
            data[lvl] = level_data
    return data           
