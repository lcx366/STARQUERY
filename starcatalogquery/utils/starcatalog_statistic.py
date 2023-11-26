from glob import glob
import os

def tiles_statistic(dir_to,tile_size):
    """
    Counts the total size and number of tile files in a star catalog directory and assesses the catalog's validity.

    Usage:
        >>> file_num,dir_size,validity = tiles_statistic(dir_to,2)
    Inputs:
        dir_to -> [str] Path where the star catalog files are stored.
        tile_size -> [int] Size of each tile in degrees.
    Outputs:
        file_num -> [int] Total number of tile files in the star catalog.
        dir_size -> [str] Total size of the star catalog, formatted as a string with appropriate units (KB, MB, GB).
        validity -> [bool] Indicates whether the star catalog is complete and safe to use.
    """
    # List all files in the specified directory
    file_list = glob(dir_to + '*')
    file_num = len(file_list)

    # Calculate the total size of all files in the directory
    dir_size = sum([os.path.getsize(file) / 1024 for file in file_list])  # Size in KB

    # Check if the number of files matches the expected count based on tile size
    expected_file_count = (180 // tile_size) * (360 // tile_size)
    validity = file_num == expected_file_count 

    # Format the directory size into a human-readable format
    if dir_size < 1024:
        dir_size = '{:.1f} KB'.format(dir_size)
    elif dir_size < 1024**2:
        dir_size = '{:.1f} MB'.format(dir_size / 1024)
    else:
        dir_size = '{:.1f} GB'.format(dir_size / 1024**2)     

    return file_num,dir_size,validity

def starcatalog_info(sc_name):
    """
    Retrieves basic information about a given star catalog.

    Usage:
        >>> star_num,mag,desc = starcatalog_info('2mass')
    Inputs:
        sc_name -> [str] Name of the star catalog.
        Available options include 'hygv3.7', 'at-hygv2.4', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.
    Outputs:
        star_num -> [str] Total number of stars in the catalog.
        mag -> [str] Magnitude range of the stars in the catalog.
        desc -> [str] Brief description of the catalog.
    """
    if sc_name == 'hygv3.7':
        star_num = '~ 0.12 Million'
        mag = '< 9'    
        desc = 'The database is a subset of the data in three major catalogs: the Hipparcos Catalog,the Yale Bright Star Catalog (5th Edition), and the Gliese Catalog of Nearby Stars (3rd Edition).'
    elif sc_name == 'at-hygv2.4':
        star_num = '~ 2.5 Million'
        mag = '< 12'    
        desc = 'The database is an augmented Tycho-2 catalog, which is essentially complete to V = 11.0 and has many fainter stars down to about V = 12.5, with Gaia DR3 distances and proper motions for nearly all of them.'
    elif sc_name == 'gsc12':
        star_num = '~ 19 Million'
        mag = '6 ~ 15' 
        desc = 'The GSC(Guide Star Catalog) I is primarily based on an all-sky, single-epoch collection of Schmidt plates.'  
    elif sc_name == 'gsc242':
        star_num = '~ 998 Million'
        mag = '< 19' 
        desc = 'The GSC(Guide Star Catalog) II is an all-sky catalog based on 1" resolution scans of the photographic Sky Survey plates, at two epochs and three bandpasses, from the Palomar and UK Schmidt telescopes (DSS).'
    elif sc_name == 'gaiadr3':
        star_num = '~ 1811 Million'
        mag = '3 ~ 21' 
        desc = 'The Gaia catalogues are star catalogues created using the results obtained by Gaia space telescope. The catalogues are released in stages that will contain increasing amounts of information.'    
    elif sc_name == '2mass':
        star_num = '~ 300 Million'
        mag = '< 14' 
        desc = 'The Two Micron All-Sky Survey(2MASS), is an astronomical survey of the whole sky in infrared light. It is conducted in the short-wavelength infrared at three distinct frequency bands (J, H, and K) near 2 micrometres.'  
    elif sc_name == 'ucac5':
        star_num = '~ 107 Million'
        mag = '< 15' 
        desc = 'The USNO CCD Astrograph Catalog (UCAC)5 is an astrometric star catalog of the United States Naval Observatory.\
        UCAC5 is based on a re-reduction of the UCAC images using the TGAS objects as reference stars to determine UCAC5 object positions. Proper motion data is then gained by comparing the UCAC5 positions with those in the GAIA DR1 catalog.'   
    elif sc_name == 'usnob':
        star_num = '~ 1042 Million'
        mag = '< 21' 
        desc = 'USNO-B is an all-sky catalog that presents positions, proper motions, magnitudes in various optical passbands, and star/galaxy estimators for 1,042,618,261 objects derived from 3,643,201,733 separate observations.'  
    else:
        raise Exception('Star catalog {:s} is not supported temporarily.'.format(sc_name))
    return star_num,mag,desc    