import os
from glob import glob
import numpy as np
from colorama import Fore

def tile2seq(tile_name):
    """
    Converts a star catalog tile file name to its corresponding sequence file name.
    The sequence file is named based on the tile's coordinates in the star catalog.

    Usage:
        >>> tile_name = 'CatalogSearch.aspx?bbox=266,-90,268,-88&format=csv&catalog=ucac5'
        >>> seq_name = tile2seq(tile_name)
        >>> print(seq_name) # starcatalogs/raw/ucac5/res2/ucac5-16153.csv
    Inputs:
        tile_name -> [str] Tile file name, e.g., 'CatalogSearch.aspx?bbox=266,-90,268,-88&format=csv&catalog=ucac5'.
    Outputs:
        seq_name -> [str] Corresponding sequence file name, e.g., 'starcatalogs/raw/ucac5/ucac5-14.csv'.
    """
    # Extracting relevant information from the tile file name
    tile_name_sep = tile_name.split('=')
    box_range = tile_name_sep[1][:-7]  # Extracting the bounding box coordinates
    scname = tile_name_sep[3]  # Extracting the star catalog name

    # Calculating sequence number based the bounding box coordinates
    ra_min, dec_min, ra_max, dec_max = np.array(box_range.split(','), dtype=int)
    tile_size = ra_max - ra_min
    seq = int(((90 - dec_max) // tile_size) * 360 / tile_size + ra_min // tile_size)

    # Constructing the sequence file name
    seq_name = f'starcatalogs/raw/{scname}/res{tile_size}/{scname}-{seq}.csv'
    
    return seq_name

def seq2tile(seq_name, tile_size):
    """
    Converts a sequence file name to the corresponding star catalog tile file format.

    Usage:
        >>> seq_name = 'starcatalogs/raw/ucac5/ucac5-14.csv'
        >>> tile_size = 2 # in degrees
        >>> tile_name = seq2tile(seq_name,tile_size)
    Inputs:
        seq_name -> [str] Sequence file name, e.g., 'starcatalogs/raw/ucac5/ucac5-14.csv'.
        tile_size -> [int] Tile size in degrees.
    Outputs:
        tile_name -> [str] Star catalog tile file name, e.g., 'CatalogSearch.aspx?bbox=266,-90,268,-88&format=csv&catalog=ucac5'.
    """
    # Extracting the sequence number and star catalog name
    seq = int(seq_name.split('-')[-1].split('.')[0])
    scname = seq_name.split('/')[-2]

    # Calculating the bounding box coordinates from the sequence number
    count_lon = 360 // tile_size
    ra_min = (seq % count_lon) * tile_size
    dec_min = 90 - ((seq // count_lon + 1) * tile_size)
    ra_max, dec_max = ra_min + tile_size, dec_min + tile_size

    # Constructing the tile file name
    tile_name = f'CatalogSearch.aspx?bbox={ra_min},{dec_min},{ra_max},{dec_max}&format=csv&catalog={scname}'
    
    return tile_name

def catalog_check(scname,tile_size,dir_from=None):
    """
    Checks and validates the star catalog tile files in a directory. 
    Re-fetch invalid files from the remote server.

    Usage:
        >>> dir_from = '/Volumes/TOSHIBA/starcatalogs/raw/ucac5/res2/'
        >>> tile_size = 2 # in [deg]
        >>> dir_from,dir_size,file_num,validity = catalog_check(scname,tile_size,dir_from)
    Inputs:
        scname -> [str] Star catalog name. Available options include 'hygv3.7', 'at-hygv2.2', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.
        tile_size -> [int] Tile size in degrees.
        dir_from -> [str,optional,default=None] Path where the star catalog tile files are stored. If None, defaults to 'starcatalogs/raw/<scname>/res<x>/', where
            - 'raw' indicates the original star catalog, relative to the reduced and simplified star catalog later.
            - '<scname>' is the catalog name.
            - '<x>' is the tile size in degrees.
    Outputs:
        dir_from -> [str] Path where the star catalog tile files are stored. 
        dir_size -> [float] Size of the star catalog
        file_num -> [int] Total number of tile files
        validity -> [bool] Validity flag

    Note: The path of the star catalog tile files should not contain spaces to prevent issues in command line operations. 
    For example, a path like '/Volumes/TOSHIBA EXT/StarCatalog/' may cause problems due to the space in 'TOSHIBA EXT'.
    """
    # Prepare a file to record URLs of invalid files for re-download
    url_file = 'url_{:s}.txt'.format(scname)
    outputf = open(url_file,'w')

    if dir_from is None:
        dir_from = 'starcatalogs/raw/{:s}/res{:d}/'.format(scname,tile_size) 

    # List all files in the directory
    file_list = glob(dir_from + '*')
    file_num = len(file_list)

    # Initialize variables for directory size and issue flag
    dir_size = 0
    issue_flag = False

    # Calculate the expected number of tiles
    num_tiles = int(180 * 360 / tile_size**2)
    print('Checking and validating the star catalog tile files for {:s}'.format(scname))

    # Iterate over each expected tile
    for i in range(num_tiles):
        file = dir_from + '{0:s}-{1:d}.csv'.format(scname, i)

        # Progress display
        desc = 'Checking {:s}{:d}{:s} of {:d}'.format(Fore.BLUE,i,Fore.RESET,num_tiles)
        print(desc,end='\r')

        # Check if file exists and is valid
        if file in file_list:
            file_size = os.path.getsize(file)
            dir_size += file_size/1024 # Convert to KB
            
            # Open and read the first line of the file
            with open(file, 'r') as fn:
                firstline = fn.readline()
                firstline_parts = firstline.strip().split(' : ')
                star_count = sum(1 for _ in fn)       

            # Validate file based on count
            if firstline_parts[0] == '#Objects found' and int(firstline_parts[1]) == star_count:
                validflag = True
            else:
                validflag, issue_flag = False, True
                os.remove(file)
        else:   
            validflag,issue_flag = False,True

        # Record the URL for re-download if invalid
        if not validflag:
            url = 'http://gsss.stsci.edu/webservices/vo/' + seq2tile(file,tile_size)
            outputf.write(file + ' ' + url + '\n')

    outputf.close()
    print()

    # Re-download invalid files if necessary
    if issue_flag:
        cmd = 'cat {:s} | xargs -P 16 -n 2 wget -c -O'.format(url_file)
        os.system(cmd)
        print('Warnings: there may still be invalid files, which need to be confirmed manually.')
        os.remove(url_file)  
    else:
        print('All star catalog tile files are valid.')

    # Format the star catalog size
    if dir_size < 1024:
        dir_size = '{:.2f} KB'.format(dir_size) 
    elif dir_size < 1024**2:
         dir_size = '{:.2f} MB'.format(dir_size/1024) 
    else:
        dir_size = '{:.2f} GB'.format(dir_size/1024**2)
        
    validity = not issue_flag

    return dir_from,dir_size,file_num,validity 