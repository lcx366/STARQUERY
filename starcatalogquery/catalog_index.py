import os,math,bisect,h5py,duckdb
import pandas as pd
import numpy as np
import healpy as hp

from glob import glob
from tqdm import tqdm
from natsort import natsorted
from functools import lru_cache

from concurrent.futures import ProcessPoolExecutor, as_completed
from .wcs import xy_catalog
from .invariantfeatures import calculate_invariantfeatures

K_RANGE = range(1, 11)  # Levels for star catalog indexing, from 1 to 10. Higher levels mean finer sky partitions.
K_RANGE_INVERSE = K_RANGE[::-1]

NSIDE_K1 = 2  # K1 corresponds to nside = 2

NSIDE = np.power(2,K_RANGE_INVERSE)  # Calculate nside for each K
NPIX = hp.nside2npix(NSIDE)  # Calculate the total number of pixels for each nside
PIXEL_SIZE = np.degrees(np.sqrt(4 * np.pi / NPIX))  # Calculate the pixel size in degrees, using the area of a pixel

# Define reasonable limits for FOV
LOWER_FOV = PIXEL_SIZE[0] * 2  # lower limits of reasonable field of view in degrees
UPPER_FOV = PIXEL_SIZE[-1] * 2  # upper limits of reasonable field of view in degrees

def find_healpix_level(lvl=None, fov_min=None):
    """
    Find the HEALPix parameters based on either a given field of view (FOV) in degrees
    or an explicit HEALPix level.

    If FOV is given, it selects a level such that each pixel’s angular size is within:
      - (¼ × FOV, ½ × FOV]

    Usage:
        >>> find_healpix_level(fov_min=5.0)
        >>> find_healpix_level(lvl=6)

    Inputs:
        lvl -> [int, optional, default=None] HEALPix level K (e.g., 6 for K6).
        fov_min -> [float, optional, default=None] Field of view (degrees) for which to choose the appropriate HEALPix level.

    Returns:
        level -> [str] Partition level as string, e.g., 'K4'.
        nside -> [int] HEALPix NSIDE (2^level).
        npix -> [int] Total number of pixels at that level.
        pixel_size -> [float] Approximate angular size of each pixel in degrees.
    """
    if lvl is not None:
        # Validate level
        if lvl not in K_RANGE_INVERSE:
            raise ValueError(f"Invalid level: {lvl}. Must be in {K_RANGE_INVERSE}")
        index = list(K_RANGE_INVERSE).index(lvl)

    elif fov_min is not None:
        # Validate FOV
        if fov_min <= LOWER_FOV or fov_min > UPPER_FOV:
            raise ValueError(f"The FOV must be between {LOWER_FOV} and {UPPER_FOV} degrees. Given FOV: {fov_min}")

        # Choose level such that pixel size ∈ (¼×FOV, ½×FOV]
        index = bisect.bisect(PIXEL_SIZE, fov_min / 2) - 1

    else:
        raise ValueError("Either 'lvl' or 'fov_min' must be provided.")

    level = f'K{K_RANGE_INVERSE[index]}'
    return level, NSIDE[index], NPIX[index], PIXEL_SIZE[index]

def estimate_mag_limit(k: int, alpha: float = 0.35, beta: float = 2.0, threshold: int = 30) -> float:
    """
    Estimate the minimum magnitude limit required for a given HEALPix order (k)
    to ensure that each pixel has, on average, at least `threshold` stars.

    Inputs:
        k -> [int] The HEALPix order (integer). The total number of pixels is 12 * 4^k.
        alpha -> [float, optional, default=0.35] The coefficient in the empirical relation N_star(m) ~ 10^(alpha * m + beta).
        beta -> [float, optional, default=2.0] The constant term in the empirical relation N_star(m) ~ 10^(alpha * m + beta).
        threshold -> [float, optional, default=30] The desired minimum average number of stars per pixel.

    Returns:
        m_limit -> [float] The estimated minimal magnitude limit that satisfies the requirement of
            `threshold` stars (on average) per HEALPix pixel at order k.
    Notes:
        1. This function uses a simple log-based empirical model for the total number
           of stars N_star(m) as a function of magnitude limit m:
               N_star(m) ~ 10^(alpha * m + beta)
        2. The condition for ensuring an average of `threshold` stars per pixel is:
               N_star(m) >= threshold * (12 * 4^k)
        3. The result is an approximate average-based estimate. In reality, star
           counts vary with Galactic latitude and other factors. For precise
           requirements (e.g., "every pixel" must have at least X stars), additional
           region-specific analysis or data is needed.
    """
    # Calculate the total number of pixels in a HEALPix grid of order k
    n_pixels = 12 * (4 ** k)

    # Calculate how many total stars are needed to achieve the threshold average in each pixel
    needed_stars = threshold * n_pixels

    # Solve N_star(m) >= needed_stars, i.e. 10^(alpha * m + beta) >= needed_stars
    # => alpha * m + beta >= log10(needed_stars)
    # => m >= (log10(needed_stars) - beta) / alpha
    m_limit = (math.log10(needed_stars) - beta) / alpha
    m_limit = math.ceil(m_limit * 2) / 2

    return m_limit

def estimate_level(m_limit: float,alpha: float = 0.35,beta: float = 2.0,threshold: int = 30) -> int:
    """
    Compute the highest HEALPix order *k* such that a star catalogue
    complete to magnitude *m_limit* still provides at least *threshold*
    stars **on average** in every pixel.

    The empirical star-count model used is
        N_star(m) ≈ 10^(alpha · m + beta)

    Inputs:
        m_limit -> [float] Limiting magnitude of the catalogue.
        alpha -> [float, optional, default=0.35] Slope in the log-star-count relation.
        beta -> [float, optional, default=2.0] Intercept in the log-star-count relation.
        threshold -> [int, optional, default=30] Desired average number of stars per pixel.
    Returns:
        k_max -> [int] Maximum integer HEALPix order *k* that satisfies the density requirement.
    """
    # Total stars predicted down to m_limit
    n_stars = 10 ** (alpha * m_limit + beta)

    # Maximum number of pixels compatible with the threshold
    max_pixels = n_stars / threshold

    # Solve 12 · 4^k ≤ max_pixels  ⇒  k ≤ 0.5 · log₂(max_pixels / 12)
    k_float = 0.5 * math.log(max_pixels / 12, 2)

    # Return the largest integer order that does not exceed the limit
    return max(1, math.floor(k_float))

def calculate_star_indices(tile_path, nsides):
    """
    Calculate the indices for each star in a single tile file across different levels.

    Usage:
        >>> tile_path = 'starcatalogs/simplified/at-hyg32/at-hyg32-100.parquet'
        >>> nsides = [2**i for i in range(4, 12)]  # From 2^4 to 2^11
        >>> indices = calculate_star_indices(tile_path, nsides)
    Inputs:
        tile_path -> [str] Path to the tile file.
        nsides -> [list] List of nsides, defining the resolution levels for Healpix tiling.
    Returns:
        indices -> [numpy.ndarray] 2D array, each row represents a star, each column represents a level.
    """
    df = pd.read_parquet(tile_path,columns=['ra', 'dec'])
    n = len(df)
    # Initialize the index array
    indices = np.empty((n, len(nsides)+1), dtype=int)

    if n > 0:
        ra, dec = df['ra'], df['dec']
        # Compute the pixel number of all stars for each level
        for i, nside in enumerate(nsides):
            pixels = hp.ang2pix(nside, ra, dec, nest=True, lonlat=True)
            indices[:, i] = pixels

        # Fill the last column with sequential IDs
        indices[:, -1] = np.arange(n)  
    
    return indices

def build_catalog_indices(dir_sc, sc_name, indices_path, k_list):
    """
    Builds a multi-resolution spatial index for a HEALPix-tiled star catalog.

    This function aggregates precomputed HEALPix indices (e.g., K6, K7, ..., K10) from all tile files
    and generates a unified catalog index table sorted for efficient spatial search.

    Usage:
        >>> dir_sc = 'starcatalogs/simplified/at-hyg32/'
        >>> sc_name = 'at-hyg32'
        >>> indices_path = 'starcatalogs/indices/at-hyg32_mag12.0_epoch2025.5.parquet'
        >>> k_list = range(6,11)
        >>> indices_path = build_catalog_indices(dir_sc, sc_name, indices_path, k_list)

    Inputs:
        dir_sc -> [str] Directory containing HEALPix-tiled star catalog files (CSV).
        sc_name -> [str] Catalog name (e.g., 'at-hyg32').
        indices_path -> [str] Path to the generated catalog index file (in Parquet format).
        k_list -> [list] The resolution levels for Healpix.

    Outputs:
        - A sorted Parquet file located at:
              <project_root>/starcatalogs/indices/<tb_name>.parquet
        - The index file includes one row per star, with columns:
              K6, K7, K8, K9, K10, K5_SUB
    """
    nsides = [2 ** i for i in k_list]

    # Collect all tile files
    file_pattern = os.path.join(dir_sc, f'{sc_name}-*.parquet')
    files = natsorted(glob(file_pattern))

    indices = []

    # Build indices for all stars across all tiles
    for file in tqdm(files, desc="Building star indices", unit="file"):
        tile_indices = calculate_star_indices(file, nsides)
        indices.append(tile_indices)

    indices = np.vstack(indices)

    # Create DataFrame with columns: K6, ..., K9, K5_SUB
    column_list = [f'K{i}' for i in k_list] + ['K5_SUB']
    df = pd.DataFrame(indices, columns=column_list)

    # Use DuckDB for efficient sorting
    con = duckdb.connect(database=':memory:')
    con.register('idx', df)

    print("Sorting index table by hierarchical HEALPix levels...")
    column_str = ", ".join(column_list)
    sorted_df = con.execute(f"""
        SELECT * FROM idx ORDER BY {column_str}
    """).fetchdf().reset_index(drop=True)
    print("Sorting complete.")

    # Export sorted index as Parquet file
    con.register('sorted_idx', sorted_df)
    con.execute(f"""
        COPY sorted_idx TO '{indices_path}' (FORMAT 'parquet');
    """)
    con.close()

def get_parent_pixels(pix_A_list: list[int], nside_A: int, nside_B: int) -> list[int]:
    """
    Given a list of HEALPix pixel IDs at resolution A (nside_A),
    compute their parent pixel IDs at a coarser resolution B (nside_B), assuming nested indexing.

    Example:
        If nside_A = 128 (K7), nside_B = 32 (K5), then each parent covers 16 subpixels,
        and parent_id = pix_id >> 4

    Inputs:
        pix_A_list -> list[int] List of pixel IDs at finer resolution nside_A.
        nside_A -> [int] Finer HEALPix resolution (e.g., 128).
        nside_B -> [int] Coarser HEALPix resolution (e.g., 32).

    Returns:
        parent_ids -> list[int] Unique parent pixel IDs at nside_B.
    """
    # Determine difference in HEALPix levels
    level_diff = (nside_A.bit_length() - 1) - (nside_B.bit_length() - 1)

    # Use set to avoid duplicates
    parent_ids = set()
    for pix in pix_A_list:
        parent = pix >> (2 * level_diff)  # Equivalent to integer division by 4^level_diff
        parent_ids.add(parent)

    return list(parent_ids)

def get_subpixels(pix_A: int, nside_A: int, nside_B: int) -> list[int]:
    """
    Given a HEALPix pixel ID at resolution A, return the list of all child pixel IDs
    at a finer resolution B, assuming nested indexing.

    Inputs:
        pix_A -> [int] Pixel ID at coarser resolution (nside_A).
        nside_A -> [int] Coarser HEALPix resolution (must be power of 2).
        nside_B -> [int] Finer HEALPix resolution (must be power of 2, >= nside_A).

    Returns:
        subpixels -> [list of int] List of child pixel IDs at resolution nside_B.
    """
    if nside_A > nside_B:
        raise ValueError("nside_A must be less than or equal to nside_B for valid subdivision.")

    if not hp.isnsideok(nside_A) or not hp.isnsideok(nside_B):
        raise ValueError("nside must be a power of 2 and HEALPix-compatible.")

    npix_A = hp.nside2npix(nside_A)
    if pix_A < 0 or pix_A >= npix_A:
        raise ValueError(f"pix_A must be within [0, {npix_A - 1}].")

    # Number of subpixels per parent in nested indexing: (nside_B / nside_A)^2
    ratio = (nside_B // nside_A) ** 2
    start_pix_B = pix_A * ratio
    end_pix_B = start_pix_B + ratio

    return list(range(start_pix_B, end_pix_B))

def get_k5_indices(indices_path: str, query_level: str, pixel_ids: list[int]) -> dict:
    """
    Maps a list of pixel IDs at a given HEALPix query level (K1–K10)
    to their corresponding K5-level pixel IDs and sub-index ranges.

    This enables fast lookup of matching stars in a K5-partitioned Parquet catalog.

    Usage:
        >>> get_k5_indices("gaiadr3_mag16.0_epoch2024.5.parquet", "K7", [1234, 5678])

    Inputs:
        indices_path -> [str] Path to the Parquet index file (generated by build_catalog_indices).
        query_level -> [str] Query level, e.g., 'K3', 'K5', 'K7', ..., up to 'K10'.
        pixel_ids -> list[int] List of pixel IDs at the query level.

    Returns:
        dict: Mapping from query pixel to a list of (K5 pixel, sub_index_list) tuples:
              - For K1–K5: { pixel_id: [(k5_pixel, None),...] }
              - For K6–K10: { pixel_id: [(k5_pixel, [row_indices...]),...] }
    Note:
        - The sub_index_list indicates the row indices in the K5 parquet file corresponding to stars under the finer pixel.
        - If None, the full partition file should be read.
    """
    lvl_num = int(query_level.replace('K', '', 1))
    result = {}

    # A. Query level K1–K5: return direct pixel-to-k5 mappings
    if lvl_num <= 5:
        nside_query = 2 ** lvl_num
        nside_k5 = 2 ** 5
        for parent_pix in pixel_ids:
            if lvl_num < 5:
                k5_list = get_subpixels(parent_pix, nside_query, nside_k5)
            else:
                k5_list = [parent_pix]
            result[parent_pix] = [(k5_val, None) for k5_val in set(k5_list)]
        return result

    # B. Query level K6–K10: need to find sub-indices from index table
    else:
        pixel_ids_str = ",".join(map(str, pixel_ids))
        con = duckdb.connect(database=':memory:')
        query = f"""
        WITH filtered AS (
            SELECT 
                {query_level} AS qlevel,
                ({query_level} >> (2 * ({lvl_num} - 5))) AS computed_k5,
                K5_SUB AS sub
            FROM '{indices_path}'
            WHERE {query_level} IN ({pixel_ids_str})
        )
        SELECT qlevel, ANY_VALUE(computed_k5) AS computed_k5, list(sub ORDER BY sub) AS sub_list
        FROM filtered
        GROUP BY qlevel
        """
        df = con.execute(query).df()
        con.close()

        for row in df.itertuples(index=False):
            key = getattr(row, "qlevel")
            computed_k5 = getattr(row, "computed_k5")
            sub_list = getattr(row, "sub_list")
            result[key] = [(computed_k5, sub_list)]
        return result

@lru_cache(maxsize=1024)
def load_parquet(file_path: str) -> pd.DataFrame:
    """
    Load a Parquet file and cache the result in memory.

    Inputs:
        file_path -> [str] Path to the Parquet file.

    Returns:
        df -> [pandas.DataFrame] The DataFrame loaded from the Parquet file.
    """
    return pd.read_parquet(file_path)

def query_parquet_tile(file_path, sub_list=None, n_stars=None):
    """
    Load a star catalog Parquet tile with optional row filtering.

    Inputs:
        file_path -> [str] Path to the Parquet tile file.
        sub_list -> [list of int, optional, default=none] A list of row indices to extract. If None, all rows will be used.
        n_stars -> [int, optional, default=None] Maximum number of rows to return. If specified, only the first `n_stars` rows
            (after subsetting) will be returned.
    Returns:
        df -> [pandas.DataFrame] A DataFrame containing the requested subset of the catalog tile.
    """
    df = load_parquet(file_path)

    if sub_list is not None:
        df = df.iloc[sub_list]

    if n_stars is not None:
        df = df.head(n_stars)

    return df

def fetch_stars_from_k5_indices(k5_indices: dict, dir_sc: str, sc_name: str, n_stars: int = None) -> pd.DataFrame:
    """
    Fetches star records from pre-partitioned Parquet catalog files (at K5 resolution),
    based on a dictionary returned by get_k5_indices().

    For each query pixel:
        - Iterates over its (K5 pixel, sub_index) mappings
        - Loads only necessary partitions and rows
        - Merges and sorts stars by brightness (magnitude)

    Usage:
        >>> dir_sc = 'starcatalogs/simplified/at-hyg32/'
        >>> sc_name = 'at-hyg32'
        >>> n_stars = 5
        >>> df = fetch_stars_from_k5_indices(indices_dict, dir_sc, sc_name, n_stars)

    Inputs:
        k5_indices -> [dict] Output from get_k5_indices().
        dir_sc -> [str] Directory containing K5-parquet star data.
        sc_name -> [str] File name of the star catalog (e.g., 'at-hyg32').
        n_stars -> [int, optional, default=None] Max number of stars to return per query pixel. If None, return all.

    Returns:
        df -> pd.DataFrame Merged and sorted list of stars from all query pixels, sorted by 'mag'.

    Notes:
        - Assumes that catalog files are in Parquet format: <prefix>-<k5_id>.parquet
        - Performs per-query-pixel selection + global brightness sort
    """
    all_frames = []
    for parent_pix, tuple_list in k5_indices.items():
        frames_for_parent = []

        for (k5_val, sub_list) in tuple_list:
            file_path = os.path.join(dir_sc, f"{sc_name}-{k5_val}.parquet")
            if not os.path.exists(file_path): continue

            df_part = query_parquet_tile(file_path, sub_list=sub_list,n_stars=n_stars)

            if not df_part.empty:
                frames_for_parent.append(df_part)

        if not frames_for_parent: continue

        # Merge results for current parent pixel and sort
        df_parent = pd.concat(frames_for_parent, ignore_index=True)
        df_parent.sort_values(by='mag', inplace=True)

        if n_stars is not None:
            df_parent = df_parent.head(n_stars)

        all_frames.append(df_parent)

    if not all_frames:
        return pd.DataFrame()

    # Final global merge and sort
    df_all = pd.concat(all_frames, ignore_index=True)
    df_all.sort_values(by='mag', inplace=True, ignore_index=True)
    return df_all

def _process_parent_pixel(parent, indices_path, query_level, nside_query, dir_sc, sc_name, n_stars):
    """
    Processes a single K1-level HEALPix parent pixel:
      1. Derives all child pixels at the given query_level under this parent.
      2. Uses precomputed Parquet index to locate relevant K5 tiles and sub-rows.
      3. Loads star records from disk using fetch_stars_from_k5_indices.
      4. Returns a 2D numpy array of [RA, DEC].

    Inputs:
        parent -> [int] K1 pixel ID (0–47).
        indices_path -> [str] Path to the star catalog index file (Parquet).
        query_level -> [str] HEALPix level to query from (e.g., 'K3', 'K8').
        nside_query -> [int] Nside of the query level (e.g., 8 for K3, 256 for K8).
        dir_sc -> [str] Directory containing partitioned Parquet star files.
        sc_name -> [str] Name of star catalog (e.g., 'at-hyg32').
        n_stars -> [int] Max number of stars to return (brightest first).

    Returns:
        tuple[int, np.ndarray] or None:
            If stars are found, returns (parent_id, np.array([[ra, dec], ...])).
            If no stars found, returns None.
    """
    child_ids = get_subpixels(parent, NSIDE_K1, nside_query)

    # Retrieve matching K5 indices for all subpixels
    idx_dict = get_k5_indices(indices_path, query_level, child_ids)

    # Fetch matching star data
    df_stars = fetch_stars_from_k5_indices(idx_dict, dir_sc, sc_name, n_stars)

    if df_stars.empty: return None

    # Extract RA/DEC as numpy array
    arr = df_stars[['ra', 'dec']].to_numpy()
    return (parent, arr)

def fetch_partitions_parallel(indices_path: str, query_level: str, dir_sc: str, sc_name: str, n_stars: int) -> dict:
    """
    Parallel version of star region fetcher.

    For each of the 48 top-level HEALPix pixels (K1), this function:
      - Extracts child pixels at the specified query level (e.g., K3, K8),
      - Finds their K5 partition info using a star catalog index (Parquet),
      - Fetches corresponding stars from Parquet files,
      - Returns top N brightest stars as [RA, DEC] pairs.

    Inputs:
        indices_path -> [str] Path to the star catalog index file (Parquet format).
        query_level -> [str] Desired HEALPix level to query from (e.g., 'K3', 'K8').
        dir_sc -> [str] Directory containing star data partition files (.parquet).
        sc_name -> [str] Name of star catalog (e.g., 'at-hyg32').
        n_stars -> [int] Max number of stars to fetch per parent region (after sorting by mag).

    Returns:
        dict[int, np.ndarray]:
            A mapping from K1 pixel ID to a numpy array of [RA, DEC] values for stars.
            Only non-empty parent regions are included.
    """
    pix_K1 = list(range(48))       # K1-level pixels range from 0 to 47
    lvl_num = int(query_level.replace("K", "", 1))
    nside_query = 2 ** lvl_num     # For example, K8 -> 256

    results = {}
    max_workers = os.cpu_count() - 1

    # Parallel task execution using all CPU cores
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _process_parent_pixel,
                parent,
                indices_path,
                query_level,
                nside_query,
                dir_sc,
                sc_name,
                n_stars
            ): parent for parent in pix_K1
        }

        # Progress tracking with tqdm
        for future in tqdm(as_completed(futures), total=len(futures), desc=f"Processing {query_level} pixels"):
            try:
                res = future.result()
                if res is not None:
                    parent, arr = res
                    results[parent] = arr
            except Exception as e:
                print(f"Error processing pixel {futures[future]}: {e}")

    return results

def sort_data_by_dec(data):
    """
    Sorts an array of [RA, Dec] coordinates by increasing absolute declination,
    with the equator first, then northern hemisphere, then southern hemisphere,
    and within each group sorted by right ascension (RA).

    Inputs:
        data -> [numpy.ndarray]
            A (n, 2) array where each row represents a point,
            with the first column as Right Ascension (RA) in degrees,
            and the second column as Declination (Dec) in degrees.

    Returns:
        sorted_data -> [numpy.ndarray]
            A (n, 2) array sorted according to:
            - increasing absolute declination (|Dec|),
            - equator first, then north, then south at the same |Dec|,
            - ascending RA within each subgroup.

        original_indices -> [numpy.ndarray]
            A (n,) array containing the original indices of the sorted entries.
    """

    # Convert input array into a pandas DataFrame
    df = pd.DataFrame(data, columns=['RA', 'Dec'])

    # Calculate absolute declination, rounding to 6 decimals to avoid floating point issues
    df['absDec'] = np.round(np.abs(df['Dec']), 6)

    # Assign hemisphere order:
    # 0 = equator (Dec == 0), 1 = north (Dec > 0), 2 = south (Dec < 0)
    df['hemi_order'] = df['Dec'].apply(lambda x: 0 if x == 0 else (1 if x > 0 else 2))

    # Sort by absDec (ascending), then hemisphere order (equator → north → south),
    # then right ascension (RA) ascending
    df = df.sort_values(['absDec', 'hemi_order', 'RA']).reset_index()

    # Extract the sorted RA, Dec values and the original indices
    sorted_data = df[['RA', 'Dec']].to_numpy()
    original_indices = df['index'].to_numpy()

    return sorted_data, original_indices

def h5_hashes(indices_path, dir_sc, sc_name, tb_name, k_min,k_max,n_stars,num_nearest_neighbors,mode_invariants,pixel_width=0.001,theta=0):
    """
    For each region of the sky, this function calculates geometric invariants of star configurations, such as triangle
    edge length ratios or quadrilateral invariants (based on methods from astrometry.net). These invariants are useful
    for tasks like star pattern blind matching.

    Usage:
        >>> indices_path = 'starcatalogs/indices/at-hyg32_mag12.0_epoch2019.5.parquet'
        >>> dir_sc = 'starcatalogs/simplified/at-hyg32/'
        >>> sc_name = 'at-hyg32'
        >>> k_min,k_max = 1,6
        >>> mode_invariants = 'triangles'
        >>> h5_file = h5_hashes(indices_path, dir_sc, sc_name,k_min,k_max, mode_invariants)
    Inputs:
        indices_path -> [str] Path to the star catalog index file (Parquet).
        dir_sc -> [str] Directory containing the star catalog files.
        sc_name -> [str] Star catalog name.
        k_min -> [int] Minimum HEALPix hierarchy level.
        k_max -> [int] Maximum HEALPix hierarchy level.
        n_stars -> [int] Number of stars per tile for each level.
        num_nearest_neighbors -> [int] Number of nearest neighbors to consider for each point.
        mode_invariants -> [str] Type of invariants to calculate ('triangles' or 'quads').
        pixel_width -> [float, optional, default=0.001] Pixel width in degrees.
        theta -> [float, optional, default=0] Rotation angle in degrees.
    Output:
        HDF5 file
    Returns:
        outh5 -> [str] Path to the generated HDF5 file.
    """
    dir_hashes = os.path.dirname(indices_path).replace('indices','hashes')
    os.makedirs(dir_hashes, exist_ok=True)
    outh5 = os.path.join(dir_hashes, f"{tb_name}_{mode_invariants}_K{k_min}_K{k_max}.h5")

    # Check if file already exists
    if not os.path.exists(outh5):
        # Calculate the center of healpix polygon
        fp_radecs = hp.pix2ang(2, range(48), nest=True, lonlat=True)
        fp_radecs = np.stack(fp_radecs).T

        # Sort fp_radecs by declination
        sorted_fp_radecs,sorted_indices = sort_data_by_dec(fp_radecs)

        hashed_data = {'fp_radecs':sorted_fp_radecs}

        with h5py.File(outh5, 'w') as f:
            # Iterate through the levels in the range
            print("Generating hdf5 data ...")
            for i in range(k_min, k_max + 1):
                query_level = f"K{i}"
                data_dict = fetch_partitions_parallel(indices_path, query_level, dir_sc, sc_name, n_stars)

                lvl_list = []

                for j,idx in enumerate(sorted_indices):
                    radec = np.unique(data_dict[idx],axis=0)
                    x, y, wcs = xy_catalog(fp_radecs[idx], radec, pixel_width, theta)
                    xy = np.stack([x, y]).T
                    invariants,asterisms,_ = calculate_invariantfeatures(xy, num_nearest_neighbors, mode_invariants)
                    grp = f.create_group(f"{query_level}/{j}")
                    grp.create_dataset('xy', data=xy)
                    grp.create_dataset('invariants', data=invariants)
                    grp.create_dataset('asterisms', data=asterisms)
                    lvl_list.append((xy, invariants, asterisms))

                hashed_data[query_level] = lvl_list
            f.create_dataset("fp_radecs", data=sorted_fp_radecs)

    return outh5,hashed_data

def read_h5_hashes(infile):
    """
    Parse the HDF5 file containing geometric invariants of star configurations for each level and partition.

    Usage:
        >>> infile = 'starcatalogs/hashes/at-hyg32_mag13.0_epoch2019.5_triangles_K1_K6.h5'
        >>> data = read_h5_hashes(infile)
    Inputs:
        infile -> [str] Path to the HDF5 file.
    Returns:
        data -> [dict] A dictionary where keys are levels ('K1' to 'K10') and 'fp_radecs'. Values for levels
                      are lists of tuples containing data for each partition ('xy', 'invariants', 'asterisms'),
                      and value for 'fp_radecs' is the corresponding central pointing.
    """
    data = {}
    with h5py.File(infile, 'r') as file:
        # Read the fixed point radec positions
        fp_radecs = data['fp_radecs'] = file['fp_radecs'][:]
        keys = np.arange(len(fp_radecs)).astype(str)

        # Iterate over each level in the file
        for lvl in list(file.keys())[:-1]:  # Assuming levels K1 to K10 are in the file
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
