import os
import numpy as np
import pandas as pd
from pathlib import Path
from gzip import GzipFile

from .utils.try_download import wget_download
from .utils.starcatalog_statistic import tiles_statistic

def catalog_download(scname,tile_size=2,dir_to=None):
    """
    Download star catalog tile files from https://outerspace.stsci.edu/display/GC.
    The celestial sphere is divided into multiple 'tiles', each representing a small sky region.
    This script downloads the tile files for a specified star catalog (e.g., GAIA DR3, UCAC5, 2MASS).

    Usage: 
        >>> dir_to,tile_size = catalog_download('gsc12')
    Inputs:
        scname -> [str] Name of the star catalog to download(e.g., 'gsc12', '2mass', 'gaiadr3', 'ucac5', 'usnob').
        tile_size -> [int, optional, default=2] Size of each tile in degrees.
        dir_to -> [str, optional, default=None] Download path of the star catalog.
    Outputs: 
        dir_to -> [str] Path where the star catalog tile files are stored. If None, defaults to 'starcatalogs/raw/<scname>/res<x>/', where
            - 'raw' indicates the original star catalog, relative to the reduced and simplified star catalog later.
            - '<scname>' is the catalog name.
            - '<x>' is the tile size in degrees.
        tile_size -> [int] Size of each tile in degrees.

    Note: The download path should not contain spaces to prevent issues in command line operations. 
    For example, a path like '/Volumes/TOSHIBA EXT/StarCatalog/' may cause problems due to the space in 'TOSHIBA EXT'.
    """  

    # Validate tile size for compatibility with celestial sphere division
    if tile_size not in [1, 2, 3, 4, 5, 6, 9, 10]: 
        raise Exception('Tile size must be in [1, 2, 3, 4, 5, 6, 9, 10] degrees')      

    # Handling different star catalogs and their respective tile sizes
    if scname == 'gsc12':
        # GSC-I: Contains approximately 20,000,000 stars with magnitudes of 6 to 15.
        if tile_size > 7.5: tile_size = 6    
    elif scname in ['gsc242','2mass']:
        # GSC-II: Current version for HST & JWST operations 2021- TBD, with 945,592,683 stars.
        # 2MASS: All-Sky IR survey Point Source Catalog.
        if tile_size > 5: tile_size = 5
    elif scname == 'gaiadr3':
        # GAIA DR3: Astrometric catalog of stellar positions and proper motions.
        tile_size = 1 
    elif scname == 'ucac5':
        # UCAC5: USNO CCD Astrograph Catalog v5.
        if tile_size > 7.5: tile_size = 6 
    elif scname == 'usnob':
        # USNO-B: Astrometric catalog v1.0.
        tile_size = 2 
    else:
        raise Exception("Star catalog '{:s}' is not supported.".format(scname))       

    # Default directory setup for saving star catalog files
    if dir_to is None:
        dir_to = 'starcatalogs/raw/{:s}/res{:d}/'.format(scname, tile_size)
        # 'raw' indicates original star catalog, '<scname>' is catalog name, '<resx>' is tile size.
    if not os.path.exists(dir_to):
        os.makedirs(dir_to)

    # Generating URLs for downloading the catalog files
    print('Generating the URL list for {:s}'.format(scname), end='...')  
    k = 0
    n_lat, n_lon = 180 // tile_size, 360 // tile_size  # Latitude and longitude divisions based on tile size

    url_file = 'url_{:s}.txt'.format(scname)  # Filename for storing URLs
    with open(url_file, 'w') as f:
        for i in range(n_lat):    
            for j in range(n_lon):
                # Calculate tile boundaries
                codec_min = tile_size * i
                codec_max = tile_size * (i + 1)
                dec_min = 90 - codec_max
                dec_max = 90 - codec_min
                ra_min = tile_size * j
                ra_max = tile_size * (j + 1)

                # Construct and write URL for each tile
                url = 'http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?bbox={:d},{:d},{:d},{:d}&format=csv&catalog={:s}'.format(ra_min,dec_min,ra_max,dec_max,scname)
                f.write(dir_to + '{:s}-{:d}.csv {:s}\n'.format(scname,k,url))
                k += 1
    print('URL list generated with {} entries.'.format(k))
               
    # Start downloading the catalog tile files using 16 threads
    print('Downloading the catalog tile files for {:s}'.format(scname),end = '...')  
    os.system('cat url_{:s}.txt | xargs -P 16 -n 2 wget -c -O'.format(scname))
    print('Finished')  
    # Clean up by removing the URL file after downloading is complete
    os.remove(url_file)

    return dir_to,tile_size  

def hyg_download(scname,tile_size=2,dir_to=None):
    """
    Downloads the HYG and AT-HYG databases from Astronexus(https://www.astronexus.com/hyg), and converts it to a tile-based format. 
    The tile-based format divides the celestial sphere into smaller regions for easier handling and analysis.

    Usage: 
        >>> dir_to,dir_size,file_num,validity,tile_size = hyg_download(5)
    Inputs:
        scname -> [str] Name of the star catalog to download. Available options include
            -  'hygv3.7': Current version, containing all stars in Hipparcos, Yale Bright Star, and Gliese catalogs (almost 120,000 stars, 14 MB)
            -  'at-hygv2.2': Current version, containing all valid stars in Tycho-2, with Gaia DR3 distances when available, as well as many fields from HYG for stars that have either HIP or Henry Draper IDs (2.55 million stars, 200 MB).
        tile_size -> [int,optinal,default=2] Size of each tile in degrees.
        dir_to -> [str,optional,default=None] Path where the star catalog is saved. Defaults to 'starcatalogs/raw/hygv37/res<x>/', where '<x>' is the tile size in degrees.
    Outputs: 
        dir_to -> [str] Path of the star catalog tile files. 
        dir_size -> [float] Size of the star catalog.
        file_num -> [int] Total number of tile files.
        validity -> [bool] Boolean indicating the validity of the star catalog.
        tile_size -> [int] Size of each tile in degrees.
    """
    # Validate tile size for compatibility with celestial sphere division.
    # N = 180//tile_size must be integer.
    if tile_size not in [1,2,3,4,5,6,9,10]: raise Exception('Tile size must be in [1,2,3,4,5,6,9,10] degrees')  

    if scname == 'hygv3.7':
        # URL of the HYG Database
        url = 'https://astronexus.com/downloads/catalogs/hygdata_v37.csv.gz'
        raw_file = 'hygdata_v37.csv'
        desc = "Downloading the HYG v3.7 database '{:s}' from The Astronomy Nexus.".format(raw_file)
    elif scname == 'at-hygv2.4':    
        # URL of the AT-HYG Database
        url = 'https://www.astronexus.com/downloads/catalogs/athyg_v24.csv.gz'
        raw_file = 'athyg_v24.csv'
        desc = "Downloading the AT-HYG v2.4 database '{:s}' from The Astronomy Nexus.".format(raw_file)
    else:
        raise Exception("For HYG/AT-HYG databases, only 'hygv3.7' and 'at-hygv2.4' are available.")
               
    # Define the directory and file name for the raw catalog
    raw_dir_to = str(Path.home()) + '/src/starcatalogs/data/'
    raw_dir_file = raw_dir_to + raw_file    
    raw_dir_file_gz = raw_dir_to + url.split('/')[-1]

    # Create the directory if it doesn't exist and download the catalog if not already present
    if not os.path.exists(raw_dir_to): os.makedirs(raw_dir_to)
    if not os.path.exists(raw_dir_file):
        wget_out = wget_download(url,raw_dir_file_gz,desc)
        g_file = GzipFile(wget_out)
        open(raw_dir_file, "wb").write(g_file.read())
        g_file.close()
        os.remove(wget_out) 
    else:
        print("Star catalog file '{:s}' is already in {:s}".format(raw_file,raw_dir_to)) 

    # Dividing the star catalog into tiles
    print('Dividing the star catalog into tiles',end='...')
    if dir_to is None:
        dir_to = 'starcatalogs/raw/{:s}/res{:d}/'.format(scname,tile_size)  
    Path(dir_to).mkdir(parents=True, exist_ok=True) 

    # Processing the catalog data
    df = pd.read_csv(raw_dir_file,skiprows=[1],dtype=str) # Skip the sun
    ra,dec = df['ra'].astype(float)*15,df['dec'].astype(float) # Convert RA to degrees

    # Latitude and longitude divisions based on tile size
    count_ra,count_dec = 360//tile_size,180//tile_size
    # Iterate over each sections and separate the catalog into tiles 
    for i in range(count_dec):
        for j in range(count_ra):
            ra_min,ra_max = tile_size*j,tile_size*(j+1)
            dec_min,dec_max = 90-tile_size*(i+1),90-i*tile_size
            
            # Selecting stars within the current tile
            ra_flag = np.abs(ra - (ra_min + ra_max)/2) < (ra_max - ra_min)/2
            dec_flag = np.abs(dec- (dec_min + dec_max)/2) < (dec_max - dec_min)/2
            
            flag = ra_flag & dec_flag 
            df_sec = df[flag]

            # Writing the tile data to a CSV file
            tile_file = dir_to + '{:s}-{:d}.csv'.format(scname,count_ra * i + j)
            with open(tile_file, 'w') as fn:
                fn.write('#Objects found: {:d}\n'.format(len(df_sec)))
                df_sec.to_csv(fn, index=False)
    print('Finished processing the catalog.')     

    # Calculating the total size and number of tile files    
    file_num,dir_size,validity = tiles_statistic(dir_to,tile_size)

    return dir_to,dir_size,file_num,validity,tile_size