import os
import numpy as np
import pandas as pd
from pathlib import Path
from gzip import GzipFile

from .utils.try_download import wget_download
from .utils.starcatalog_statistic import tiles_statistic

def catalog_download(scname,tile_size=None,dir_to=None):
    """
    Download starcatalog tile files from https://outerspace.stsci.edu/display/GC

    Usage: 
        >>> dir_to,tile_size = catalog_download('gsc12')

    Inputs:
        scname -> [str] Name of the starcatalog to download. Available starcatalogs include 'hygv35', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.
        tile_size -> [int,optinal,default=None] size of tile in [deg]
        dir_to -> [str,optional,default=None] Download path of the starcatalog

    Outputs: 
        dir_to -> [str] Path of the starcatalog tile files. If None, the path is assigned to 'starcatalogs/raw/<scname>/<resx>/' by default.
        tile_size -> [int] Size of tile in [deg]. If None, the size of the tile is automatically assigned a feasible maximum according to the name of the star catalog.

    Note: In order to avoid unnecessary troubles when downloading the starcatalog tile files from the command line, it is best not to contain spaces in <dir_from>, such as '/Volumes/TOSHIBA EXT/StarCatalog/'.
    """

    if tile_size is None:
        tile_size = 2
    else:    
        # for pyshtools, N = 180//tile_size must be even.
        if tile_size not in [1,2,3,5,6,9,10]: raise Exception('tile size must be in [1,2,3,5,6,9,10] deg')   

    if scname == 'gsc12':
        # version for updated Astrometry 
        # GSC-I contained approximately 20,000,000 stars with apparent magnitudes of 6 to 15.
        if tile_size > 7.5: tile_size = 6    
    elif scname in ['gsc242','2mass']:
        # gsc242: Current version for HST & JWST operations 2021- TBD  
        # GSC-II contains 945,592,683 stars out to magnitude 21.
        # 2mass: 2MASS All-Sky IR survey Point Source Catalog   
        if tile_size > 5: tile_size = 5
    elif scname == 'gaiadr3':
        # gaiadr3: GAIA Astrometric catalog of stellar positions and proper motions DR3
        tile_size = 1 
    elif scname == 'ucac5':
        # ucac5: USNO CCD Astrograph Catalog v5
        if tile_size > 7.5: tile_size = 6 
    elif scname == 'usnob':
        # usnob: USNO-B Astrometric catalog v1.0
        tile_size = 2    

    print('Generating the url list for {:s}'.format(scname),end = '...')  

    if dir_to is None:
        dir_to = 'starcatalogs/raw/{:s}/res{:d}/'.format(scname,tile_size)
    if not os.path.exists(dir_to): os.makedirs(dir_to)
    
    k = 0
    n_lat,n_lon = 180//tile_size,360//tile_size

    url_file = 'url_{:s}.txt'.format(scname)
    f = open(url_file, 'w')

    for i in range(n_lat):    
        for j in range(n_lon):
            codec_min = tile_size*i
            codec_max = tile_size*(i+1)
            dec_min = 90 - codec_max
            dec_max = 90 - codec_min
            ra_min = tile_size*j
            ra_max = tile_size*(j+1)
            url = 'http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?bbox={:d},{:d},{:d},{:d}&format=csv&catalog={:s}'.format(ra_min,dec_min,ra_max,dec_max,scname)
            f.write(dir_to + '{:s}-{:d}.csv {:s}\n'.format(scname,k,url))
            k+=1    
    f.close()
    print('Finished')
               
    # Start the downloading by 16 threads
    print('Downloading the catalog tile files for {:s}'.format(scname),end = '...')  
    os.system('cat url_{:s}.txt | xargs -P 16 -n 2 wget -c -O'.format(scname))
    print('Finished')  
    os.remove(url_file)

    return dir_to,tile_size  

def hygv35_download(tile_size=None,dir_to=None):
    """
    Download the HYG v35 database and convert it to the tile mode.

    Usage: 
        >>> dir_to,tile_size = hygv35_download(5)

    Inputs:
        tile_size -> [int,optinal,default=None] size of tile in [deg]
        dir_to -> [str,optional,default=None] Download path of the starcatalog    

    Outputs: 
        dir_to -> [str] Path of the starcatalog tile files. If None, the path is assigned to 'starcatalogs/raw/hygv35/<resx>/' by default.
        dir_size -> [float] The size of the star catalog.
        file_num -> [int] Total number of the tile files.
        validity -> [bool] The validity of the star catalog. 
        tile_size -> [int] Size of tile in [deg]. If None, default = 2.
    """

    if tile_size is None:
        tile_size = 2
    else:    
        # for pyshtools, N = 180//tile_size must be even.
        if tile_size not in [1,2,3,5,6,9,10]: raise Exception('tile size must be in [1,2,3,5,6,9,10] deg')  

    raw_dir_to = str(Path.home()) + '/src/starcatalogs/data/'
    raw_file = 'hygdata_v35.csv'
    raw_dir_file = raw_dir_to + raw_file

    url = 'https://astronexus.com/downloads/catalogs/hygdata_v35.csv.gz'

    raw_dir_file_gz = raw_dir_to + url.split('/')[-1]

    if not os.path.exists(raw_dir_to): os.makedirs(raw_dir_to)
    if not os.path.exists(raw_dir_file):
        desc = "Downloading the HYG v35 database '{:s}' from The Astronomy Nexus.".format(raw_file)
        wget_out = wget_download(url,raw_dir_file_gz,desc)
        g_file = GzipFile(wget_out)
        open(raw_dir_file, "wb").write(g_file.read())
        g_file.close()
        os.remove(wget_out) 

    else:
        print("Star catalog file '{:s}' is already in {:s}".format(raw_file,raw_dir_to)) 

    # Now divide the entire star catalog file into multiple tile files
    print('Divide the entire star catalog file into multiple tile files',end='...')
    if dir_to is None:
        dir_to = 'starcatalogs/raw/hygv35/res{:d}/'.format(tile_size)  
    Path(dir_to).mkdir(parents=True, exist_ok=True) 

    df = pd.read_csv(raw_dir_file,skiprows=[1]) # remove the sun
    ra,dec = df['ra']*15,df['dec']  

    count_ra = 360//tile_size
    count_dec = 180//tile_size
    for i in range(count_dec):
        for j in range(count_ra):
            ra_min,ra_max = tile_size*j,tile_size*(j+1)
            dec_min,dec_max = 90-tile_size*(i+1),90-i*tile_size
            
            ra_flag = np.abs(ra - (ra_min + ra_max)/2) < (ra_max - ra_min)/2
            dec_flag = np.abs(dec- (dec_min + dec_max)/2) < (dec_max - dec_min)/2
            
            flag = ra_flag & dec_flag 
            df_sec = df[flag]

            fn = open(dir_to+'hygv35-{:d}.csv'.format(count_ra *i+j), 'w')
            fn.write('#Objects found : {:d}\n'.format(len(df_sec)))
            df_sec.to_csv(fn,index=False)   
            fn.close() 
    print('Finished')         

    # calculate total size and numbers of tile files    
    file_num,dir_size,validity = tiles_statistic(dir_to,tile_size)

    return dir_to,dir_size,file_num,validity,tile_size