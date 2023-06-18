from glob import glob
import os

def tiles_statistic(dir_to,tile_size):
    """
    Count the total size and number of tiles of the star catalog, and judge the validity of the star catalog accordingly.

    Usage:
        >>> dir_to = '/Volumes/TOSHIBA/starcatalog/reduced/ucac5/res2/''
        >>> file_num,dir_size,validity = tiles_statistic(dir_to,2)

    Inputs:
        dir_to -> [str] The storage path of the star catalog file
        tile_size -> [int] Size of tile in [deg]

    Outputs:
        file_num -> [int] The total number of tiles of the star catalog
        dir_size -> [float] The total size of the star catalog
        validity -> [bool] If True, the star catalog is safe to use, otherwise the integrity of the star catalog needs to be checked.
    """
    file_list = glob(dir_to+'*')
    file_num = len(file_list)

    dir_size = sum([os.path.getsize(file)/1024 for file in file_list])

    if file_num == 180//tile_size*360//tile_size:
        validity = True
    else:
        validity = False

    if dir_size < 1024:
        dir_size = '{:.1f} KB'.format(dir_size) 
    elif dir_size < 1024**2:
        dir_size = '{:.1f} MB'.format(dir_size/1024) 
    else:
        dir_size = '{:.1f} GB'.format(dir_size/1024**2)    

    return file_num,dir_size,validity

def starcatalog_info(sc_name):
    """
    Given a star catalog name, get basic information about the catalog.

    Usage:
        >>> star_num,mag,description = starcatalog_info('2mass')

    Inputs:
        sc_name -> [str] Name of the starcatalog. Available starcatalogs include 'hygv35', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.

    Outputs:
        star_num -> [str] The total number of stars contained in the catalog
        mag -> [str] The magnitude range of the star catalog
        description -> [str] A brief description of the catalog      
    """
    if sc_name == 'hygv35':
        star_num = '~ 0.12 Million'
        mag = '< 9'    
        description = 'The database is a subset of the data in three major catalogs: the Hipparcos Catalog,the Yale Bright Star Catalog (5th Edition), and the Gliese Catalog of Nearby Stars (3rd Edition).'
    elif sc_name == 'gsc12':
        star_num = '~ 19 Million'
        mag = '6 ~ 15' 
        description = 'The GSC(Guide Star Catalog) I is primarily based on an all-sky, single-epoch collection of Schmidt plates.'  
    elif sc_name == 'gsc242':
        star_num = '~ 998 Million'
        mag = '< 19' 
        description = 'The GSC(Guide Star Catalog) II is an all-sky catalog based on 1" resolution scans of the photographic Sky Survey plates, at two epochs and three bandpasses, from the Palomar and UK Schmidt telescopes (DSS).'
    elif sc_name == 'gaiadr3':
        star_num = '~ 1811 Million'
        mag = '3 ~ 21' 
        description = 'The Gaia catalogues are star catalogues created using the results obtained by Gaia space telescope. The catalogues are released in stages that will contain increasing amounts of information.'    
    elif sc_name == '2mass':
        star_num = '~ 300 Million'
        mag = '< 14' 
        description = 'The Two Micron All-Sky Survey(2MASS), is an astronomical survey of the whole sky in infrared light. It is conducted in the short-wavelength infrared at three distinct frequency bands (J, H, and K) near 2 micrometres.'  
    elif sc_name == 'ucac5':
        star_num = '~ 107 Million'
        mag = '< 15' 
        description = 'The USNO CCD Astrograph Catalog (UCAC)5 is an astrometric star catalog of the United States Naval Observatory.\
        UCAC5 is based on a re-reduction of the UCAC images using the TGAS objects as reference stars to determine UCAC5 object positions. Proper motion data is then gained by comparing the UCAC5 positions with those in the GAIA DR1 catalog.'   
    elif sc_name == 'usnob':
        star_num = '~ 1042 Million'
        mag = '< 21' 
        description = 'USNO-B is an all-sky catalog that presents positions, proper motions, magnitudes in various optical passbands, and star/galaxy estimators for 1,042,618,261 objects derived from 3,643,201,733 separate observations.'  
    else:
        raise Exception('Star catalog {:s} is not supported temporarily.'.format(sc_name))
    return star_num,mag,description    