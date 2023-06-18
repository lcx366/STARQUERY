import os
from glob import glob
import numpy as np
from colorama import Fore

def tile2seq(tile_name):
    """
    Parse the name of the star catalog tile file and return the name of the corresponding sequence file.

    Usage:
        >>> tile_name = 'CatalogSearch.aspx?bbox=266,-90,268,-88&format=csv&catalog=ucac5'
        >>> seq_name = tile2seq(tile_name)
        >>> print(seq_name) # 'starcatalogs/raw/ucac5/ucac5-14.csv'

    Inputs:
        tile_name -> [str] Name of the star catalog tile file to parse, such as 'CatalogSearch.aspx?bbox=266,-90,268,-88&format=csv&catalog=ucac5'.

    Outputs:
        seq_name -> [str] Name of the corresponding sequence file, such as 'starcatalogs/raw/ucac5/ucac5-14.csv'.
    """
    tile_name_sep = tile_name.split('=')
    box_range = tile_name_sep[1][:-7]
    scname = tile_name_sep[3]
    
    ra_min,dec_min,ra_max,dec_max = np.array(box_range.split(','),dtype=int)
    
    codec_min = 90 - dec_max
    codec_max = 90 - dec_min

    tile_size = ra_max - ra_min

    ra_seq =  ra_min//tile_size
    codec_seq = codec_min//tile_size
    seq = int(codec_seq*360/tile_size + ra_seq)
    seq_name = 'starcatalogs/raw/{0:s}/res{1:d}/{0:s}-{2:d}.csv'.format(scname,tile_size,seq)
    
    return seq_name

def seq2tile(seq_name,tile_size):
    """
    Parse the name of the sequence star catalog file and return the corresponding name of the star catalog tile file.

    Usage:
        >>> seq_name = 'starcatalogs/raw/ucac5/ucac5-14.csv'
        >>> tile_size = 2 # in [deg]
        >>> tile_name = seq2tile(seq_name,tile_size)
        >>> print(tile_name) # 'CatalogSearch.aspx?bbox=266,-90,268,-88&format=csv&catalog=ucac5'
    
    Inputs:
        seq_name -> [str] Name of the corresponding sequence file to parse, such as 'starcatalogs/raw/ucac5/ucac5-14.csv'.   
        tile_size -> [int] size of tile in [deg]

    Outputs:
        tile_name -> [str] Name of the star catalog tile file, such as 'CatalogSearch.aspx?bbox=266,-90,268,-88&format=csv&catalog=ucac5'.    
    """
    
    seq_name_sep = seq_name.split('-')
    scname = seq_name_sep[0].split('/')[-1]
    seq = int(seq_name_sep[1][:-4])  
    
    count_lon = 360//tile_size
    ra_min = seq%count_lon*tile_size
    ra_max = ra_min + tile_size
    codec_min = seq//count_lon*tile_size
    codec_max = codec_min + tile_size
    dec_min = 90 - codec_max
    dec_max = 90 - codec_min
    tile_name = 'CatalogSearch.aspx?bbox={:d},{:d},{:d},{:d}&format=csv&catalog={:s}'.format(ra_min,dec_min,ra_max,dec_max,scname)
    
    return tile_name

def catalog_check(scname,tile_size,dir_from=None):
    """
    Check the validity of the star catalog tile files and re-fetch invalid files from the remote server.

    Usage:
        >>> dir_from = '/Volumes/TOSHIBA/starcatalogs/raw/ucac5/res2/'
        >>> tile_size = 2 # in [deg]
        >>> dir_from,dir_size,file_num,validity = catalog_check(scname,tile_size,dir_from)

    Inputs:
        scname -> [str] Name of the starcatalog. Available starcatalogs include 'hygv35', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.
        tile_size -> [int] size of tile in [deg]
        dir_from -> [str,optional,default=None] Path of the starcatalog tile files. 

    Outputs:
        dir_from -> [str] Path of the starcatalog tile files. If None, the path is assigned to 'starcatalogs/raw/<scname>/<resx>/' by default.
        dir_size -> [float] The size of the star catalog.
        file_num -> [int] Total number of the tile files.
        validity -> [bool] The validity of the star catalog.

    Note: In order to avoid unnecessary troubles when downloading the starcatalog tile files from the command line, it is best not to contain spaces in <dir_from>, such as '/Volumes/TOSHIBA EXT/StarCatalog/'.
    """
    url_file = 'url_{:s}.txt'.format(scname)
    outputf = open(url_file,'w')

    if dir_from is None:
        dir_from = 'starcatalogs/raw/{:s}/res{:d}/'.format(scname,tile_size) 

    file_list = glob(dir_from+'*')
    file_num = len(file_list)

    dir_size = 0
    issue_flag = False

    num_tiles = int(180*360/tile_size**2)
    print('Checking the completeness and validity of the starcatalog tile files for {:s}'.format(scname)) 

    for i in range(num_tiles):
        file = dir_from+'{0:s}-{1:d}.csv'.format(scname,i)

        desc = 'Checking {:s}{:d}{:s} of {:d}'.format(Fore.BLUE,i,Fore.RESET,num_tiles)
        print(desc,end='\r')

        if file in file_list:
            file_size = os.path.getsize(file)
            dir_size += file_size/1024 # convert file_size in bytes to Kb
            
            fn = open(file,'r')
            firstline = fn.readline()
            firstline_parts = firstline.strip().split(' : ')

            for count, line in enumerate(fn): pass
            fn.close()

            if firstline_parts[0] == '#Objects found':
                if int(firstline_parts[1]) == count:
                    validflag = True
                else:
                    validflag,issue_flag = False,True 
                    os.remove(file) 
            else:
                validflag,issue_flag = False,True 
                os.remove(file)                  
        else:   
            validflag,issue_flag = False,True

        if not validflag:
            url = 'http://gsss.stsci.edu/webservices/vo/' + seq2tile(file,tile_size)
            outputf.write(file + ' ' + url + '\n')

    outputf.close()
    print()

    if issue_flag:
        cmd = 'cat {:s} | xargs -P 16 -n 2 wget -c -O'.format(url_file)
        os.system(cmd)
        print('Warnings: there may still be invalid files, which need to be confirmed manually.')
        os.remove(url_file)  
    else:
        print('All star catalog tile files are valid.')

    if dir_size < 1024:
        dir_size = '{:.2f} KB'.format(dir_size) 
    elif dir_size < 1024**2:
         dir_size = '{:.2f} MB'.format(dir_size/1024) 
    else:
        dir_size = '{:.2f} GB'.format(dir_size/1024**2)
    validity = not issue_flag

    return dir_from,dir_size,file_num,validity 