import os
from glob import glob
from pyarrow.parquet import ParquetFile

NUM_TILES = 12288  # Default number of HEALPix pixels for (K=5, NSIDE=32)

def tiles_statistic(dir_to,file_type='csv'):
    """
    Counts the total size and number of tile files in a star catalog directory and assesses the catalog's validity.

    Usage:
        >>> dir_to = 'starcatalogs/raw/at-hyg24/'
        >>> file_num, dir_size, validity = tiles_statistic(dir_to)
    Inputs:
        dir_to -> [str] Directory where the star catalog files are stored.
    Outputs:
        file_num -> [int] Total number of the tile files for the star catalog.
        dir_size -> [str] Total size of the star catalog, formatted as a string with appropriate units (KB, MB, GB).
        validity -> [bool] Validity of the star catalog. It indicates whether the star catalog is complete and safe to use.
    """
    # List all files in the catalog directory
    file_list = glob(os.path.join(dir_to, '*'))
    file_num = len(file_list)

    # Calculate the total size of all files in the directory
    dir_size = sum([os.path.getsize(file) / 1024 for file in file_list])  # Size in KB

    # Count total lines across all files
    total_lines = 0
    if file_type == 'csv':
        chunk_size = 1024 * 1024  # 1 MB
        for file in file_list:
            with open(file, 'rb') as f:
                while chunk := f.read(chunk_size):
                    total_lines += chunk.count(b'\n')
    elif file_type == 'parquet':
        for file in file_list:
            pf = ParquetFile(file)
            total_lines += pf.metadata.num_rows
    else:
        raise ValueError('Unrecognized file type.')

    # Check if the number of files matches the expected count based on tile size
    expected_file_count = NUM_TILES  # Expected file count for a full catalog
    validity = file_num == expected_file_count

    # Format the directory size into a human-readable format
    if dir_size < 1024:
        dir_size = '{:.1f} KB'.format(dir_size)
    elif dir_size < 1024 ** 2:
        dir_size = '{:.1f} MB'.format(dir_size / 1024)
    else:
        dir_size = '{:.1f} GB'.format(dir_size / 1024 ** 2)

    return file_num, dir_size, total_lines, validity

def starcatalog_info(sc_name):
    """
    Retrieves basic information about a given star catalog.

    Usage:
        >>> star_num, mag, desc = starcatalog_info('2mass')
    Inputs:
        sc_name -> [str] Name of the star catalog.
        Supported catalogs include 'hyg37', 'at-hyg24', 'gaiadr3', 'gsc30', 'ucac5', 'usnob', '2mass', etc.
    Outputs:
        star_num -> [str] Total number of stars in the catalog.
        mag -> [str] Magnitude range of the stars in the catalog.
        desc -> [str] Brief description of the catalog.
    """
    catalogs = {
        'hyg41': {
            'star_num': '~ 0.12 Million',
            'mag': '< 9',
            'desc': 'The database is a subset of the data in three major catalogs: the Hipparcos Catalog, the Yale Bright Star Catalog (5th Edition), and the Gliese Catalog of Nearby Stars (3rd Edition).'
        },
        'at-hyg32': {
            'star_num': '~ 2.5 Million',
            'mag': '< 12',
            'desc': 'The database is an augmented Tycho-2 catalog, which is essentially complete to V = 11.0 and has many fainter stars down to about V = 12.5, with Gaia DR3 distances and proper motions for nearly all of them.'
        },
        'gaiadr3': {
            'star_num': '~ 1811 Million',
            'mag': '3 ~ 21',
            'desc': 'The Gaia catalogues are star catalogues created using the results obtained by Gaia space telescope. The catalogues are released in stages that will contain increasing amounts of information.'
        },
        'gsc30': {
            'star_num': '~ 998 Million',
            'mag': '< 19',
            'desc': 'The GSC (Guide Star Catalog) III is an all-sky catalog based on 1" resolution scans of the photographic Sky Survey plates, at two epochs and three bandpasses, from the Palomar and UK Schmidt telescopes (DSS).'
        },
        'ucac5': {
            'star_num': '~ 107 Million',
            'mag': '7.5 ~ 15.5',
            'desc': 'The UCAC5 (U.S. Naval Observatory CCD Astrograph Catalog, Fifth Edition) is an extensive star catalog that builds upon the previous UCAC4 catalog. It integrates data from the Gaia DR1, enhancing the accuracy and precision of the proper motions of stars listed in the catalog.'
        },
        'usnob': {
            'star_num': '~ 1042 Million',
            'mag': '< 21',
            'desc': 'USNO-B is an all-sky catalog that presents positions, proper motions, magnitudes in various optical passbands, and star/galaxy estimators for 1,042,618,261 objects derived from 3,643,201,733 separate observations.'
        },
        '2mass': {
            'star_num': '~ 300 Million',
            'mag': '< 14',
            'desc': 'The Two Micron All-Sky Survey (2MASS), is an astronomical survey of the whole sky in infrared light. It is conducted in the short-wavelength infrared at three distinct frequency bands (J, H, and K) near 2 micrometres.'
        }
    }

    if sc_name in catalogs:
        star_num = catalogs[sc_name]['star_num']
        mag = catalogs[sc_name]['mag']
        desc = catalogs[sc_name]['desc']
    else:
        supported_catalogs = ', '.join(catalogs.keys())
        raise Exception(
            f'Star catalog {sc_name} is not supported temporarily. Supported catalogs are: {supported_catalogs}')

    return star_num, mag, desc
