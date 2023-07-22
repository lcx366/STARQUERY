import os,warnings,copy,h5py
from glob import glob
import pandas as pd
import numpy as np
from pathlib import Path
from astropy.time import Time
from astropy.coordinates import SkyCoord
from colorama import Fore
import healpy as hp
from scipy.spatial import KDTree

from .catalog_download import catalog_download,hygv35_download
from .catalog_check import catalog_check
from .utils.starcatalog_statistic import tiles_statistic,starcatalog_info
from .utils.df2info import df2info
from .catalog_query import search_box,search_cone,search_box_magpm,search_cone_magpm,box2seqs,cone2seqs,_load_files
from .tiles_draw import search_draw
from .wcs import xy_catalog
from .invariantfeatures import _generate_invariants

class StarCatalog(object):

    def get(sc_name,tile_size=None,dir_to=None):
        """
        Grab star catalog data files from remote servers.

        Usage:
            >>> from starcatalogquery import StarCatalog
            >>> hygv35_raw = StarCatalog.get('hygv35',5)
        Inputs:
            sc_name -> [str] Name of the star catalog. Available options include 'hygv35', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.
            tile_size -> [int,optinal,default=None] Geometric size of the tile in [deg]. If None, it is assigned an allowed maximum automatically based on the star catalog.
            dir_to -> [str,optional,default=None] Path of the star catalog files grabbed. If None, the it is assigned to a build-in directory by default.
        Outputs:
            Instance of class StarCatalogRaw
        """
        if sc_name == 'hygv35': # for HYG v35 star catalog
            dir_to,dir_size,file_num,validity,tile_size = hygv35_download(tile_size,dir_to) 

        else: # for other star catalogs

            if dir_to is None:
                dir_to = 'starcatalogs/raw/{:s}/res{:d}/'.format(sc_name,tile_size) 

            if os.path.exists(dir_to):
                dir_to,dir_size,file_num,validity = catalog_check(sc_name,tile_size,dir_to)
            else:
                dir_to,tile_size = catalog_download(sc_name,tile_size,dir_to)  
                dir_to,dir_size,file_num,validity = catalog_check(sc_name,tile_size,dir_to) 

            dir_to,dir_size,file_num,validity = catalog_check(sc_name,tile_size,dir_to)  

        stars_num,mag,description = starcatalog_info(sc_name)
        dict_values = dir_to,dir_size,file_num,validity,sc_name,'{:d} deg'.format(tile_size),'raw',stars_num,mag,description
        dict_keys = 'tiles_path','sc_size','tiles_num','validity','sc_name','tile_size','_mode','stars_num','mag','description'
        info = dict(zip(dict_keys, dict_values))

        return StarCatalogRaw(info)  

    def load(dir_from=None):
        """
        Load the star catalog files from the local database.

        Usage:
            >>> from starcatalogquery import StarCatalog
            >>> # load the raw star catalog HYGv3.5
            >>> dir_from_raw = '/Volumes/TOSHIBA/starcatalogs/raw/hygv35/res5/'
            >>> hygv35_raw = StarCatalog.load(dir_from_raw)
            >>>
            >>> # load the reduced star catalog HYGv3.5
            >>> dir_from_reduced = '/Volumes/TOSHIBA/starcatalogs/reduced/hygv35/res5/'
            >>> hygv35_reduced = StarCatalog.load(dir_from_reduced)
            >>>
            >>> # load the simplified star catalog HYGv3.5
            >>> dir_from_simplified = '/Volumes/TOSHIBA/starcatalogs/simplified/hygv35/res5/mag9.0/epoch2022.0/'
            >>> hygv35_simplified = StarCatalog.load(dir_from_simplified)
        Inputs:
            dir_from -> [str,optional,default=None] Directory of the star catalog files. If None, it is assigned to a build-in directory by default.
        Outputs:
            Instance of class StarCatalog
        """
        _mode,sc_name,tile_size = dir_from.split('starcatalogs/')[1].split('/')[:3]
        tile_size = int(tile_size[3:])

        # _mode -> [str] Type of the star catalog. Available options include 'raw', 'reduced', 'simplified', where
        # 'raw' represents the original star catalog, which covers all information about the stars,
        # 'reduced' represents the reduced star catalog, which contains the position, proper motion, apparent magnitude, epoch of the stars,
        # 'simplified' represents the minimalist star catalog, which only includes the position and apparent magnitude of the stars at a specific epoch.
        # sc_name -> [str] Name of the star catalog. Available options include 'hygv35', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.
        # tile_size -> [int] Geometric size of the tile in [deg]

        if _mode == 'raw':
            starcatalog = StarCatalogRaw.load(sc_name,tile_size,dir_from)
        elif _mode == 'reduced':
            starcatalog = StarCatalogReduced.load(sc_name,tile_size,dir_from)
        elif _mode == 'simplified':    
            mag_threshold,epoch = np.array(dir_from.split('mag')[1][:-1].split('/epoch'),dtype=float) 
            starcatalog = StarCatalogSimplified.load(sc_name,tile_size,mag_threshold,epoch,dir_from)

        return starcatalog 

    def read_h5_indices(infile):
        """
        Read and parse the h5-formatted star catalog indicea file, which records the center pointing, pixel coordinates of the stars, triangle invariants and asterism indices of each sky area.

        Usage:
            >>> from starcatalogquery import StarCatalog
            >>> infile_h5 = 'starcatalogs/indices/hygv35/k2_mag9.0_mcp30_2022.0.h5'
            >>> fp_radecs,stars_xy,stars_invariants,stars_asterisms = read_h5_indices(infile_h5)
        Inputs:
            infile_h5 -> [str] h5-formatted star catalog indices file  
        Outputs:
            fp_radecs -> [2d array of float] Center pointing of each sky area in form of [[RA0,DEC0],..[RAn,DECn]] in [deg]
            stars_xy -> [list of 2d array] Pixel coordinates of stars in each sky area
            stars_invariants -> [list of 2d array] Triangle invariants in each sky area
            stars_asterisms -> [list of 2d array] Asterism indices corresponding to the triangle invariants in each sky area       
        """
        fin = h5py.File(infile,'r')

        # read data
        fp_radecs = fin['fp_radecs'][:]
        stars_xy,stars_invariants, stars_asterisms = [],[],[]
        for j in range(len(fp_radecs)):
            stars_xy.append(fin['stars_xy/'+str(j)][:])
            stars_invariants.append(fin['stars_invariants/'+str(j)][:])
            stars_asterisms.append(fin['stars_asterisms/'+str(j)][:])
    
        # close the file
        fin.close()   
        return fp_radecs,stars_xy,stars_invariants,stars_asterisms  

class StarCatalogRaw(object):
    """
    Class StarCatalogRaw

    Attributes:
        - tiles_path: Path of the star catalog files. 
        - tiles_num: Total number of the tile files.
        - tile_size: Geometric size of the tile in [deg]
        - sc_size: File size of the star catalog.
        - validity: Validity of the star catalog.
        - sc_name: Name of the star catalog. Available options include 'hygv35', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.
        - _mode: Type of the star catalog. Available options include 'raw', 'reduced', 'simplified', where
            'raw' represents the original star catalog, which covers all information about the stars,
            'reduced' represents the reduced star catalog, which contains the position, proper motion, apparent magnitude, epoch of the stars,
            'simplified' represents the minimalist star catalog, which only includes the position and apparent magnitude of the stars at a specific epoch.
        - stars_num: Total number of stars in the catalog
        - mag: Range of apparent magnitudes for the catalog
        - description: Catalog Summary
    Methods:
        - load: Load the raw star catalog files from the local database.
        - reduce: Reduce the raw star catalog so that it only contains necessary information of stars such as the celestial position, proper motion, apparent magnitude, and epoch, etc.
        - search_box: Perform a rectangle search on the raw star catalog.
        - search_cone: Perform a cone search on the raw star catalog.
        - _search_draw: Visualize the scope of the search area and the coverage of the corresponding tiles.    
    """   
    def __init__(self,info):  

        self.info = info

        for key in info.keys():
            setattr(self, key, info[key])

    def __repr__(self):
    
        return 'Instance of class StarCatalogRaw'

    def load(sc_name,tile_size,dir_from=None):
        """
        Load the raw star catalog files from the local database.

        Usage:
            >>> from starcatalogquery import StarCatalogRaw
            >>> # load the raw star catalog HYGv3.5
            >>> dir_from_raw = '/Volumes/TOSHIBA/starcatalogs/raw/hygv35/res5/'
            >>> hygv35_raw = StarCatalogRaw.load('hygv35',5,dir_from_raw)
        Inputs:
            sc_name -> [str] Name of the star catalog. Available options include 'hygv35', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.
            tile_size -> [int] Geometric size of the tile in [deg]
            dir_from -> [str,optional,default=None] Directory of the star catalog files. If None, it is assigned to a build-in directory by default.
        Outputs:
            Instance of class StarCatalogRaw
        """
        if dir_from is None: dir_from = 'starcatalogs/raw/{:s}/res{:d}/'.format(sc_name,tile_size)    
        if not os.path.exists(dir_from): raise Error('The storage directory for the catalog {:s} does not exist.'.format(sc_name))  

        # calculate total size and numbers of tile files    
        file_num,dir_size,validity = tiles_statistic(dir_from,tile_size) 
        stars_num,mag,description = starcatalog_info(sc_name)

        dict_values = dir_from,dir_size,file_num,validity,sc_name,'{:d} deg'.format(tile_size),'raw',stars_num,mag,description
        dict_keys = 'tiles_path','sc_size','tiles_num','validity','sc_name','tile_size','_mode','stars_num','mag','description'
        info = dict(zip(dict_keys, dict_values))

        return StarCatalogRaw(info) 

    def reduce(self,dir_reduced=None):
        """
        Reduce the raw star catalog so that it only contains necessary information of stars such as the celestial position, proper motion, apparent magnitude, and epoch, etc.

        Usage:
            >>> from starcatalogquery import StarCatalogRaw
            >>> dir_from_raw = '/Volumes/TOSHIBA/starcatalogs/raw/hygv35/res5/'
            >>> hygv35_raw = StarCatalogRaw.load('hygv35',5,dir_from_raw)
            >>> hygv35_reduced = hygv35_raw.reduce()
        Inputs:
            dir_reduced -> [str,optional,default=None] Directory of the reduced star catalog files. If None, it is assigned to a build-in directory by default.
        Outputs:
            Instance of class StarCatalogReduced       
        """
        info = self.info.copy()
        tile_size = int(self.tile_size.split(' ')[0])
        tiles_path = self.tiles_path
        if dir_reduced is None:
            dir_reduced = tiles_path.replace('raw','reduced')
        Path(dir_reduced).mkdir(parents=True, exist_ok=True)    

        file_list = glob(tiles_path+'*') 
        sc_name = self.sc_name
        print('Reducing the star catalog {:s}, which may take a considerable amount of time'.format(sc_name))  

        if sc_name == 'hygv35':
            j = 1
            for tile_file in file_list:
                desc = 'Reducing {:s}{:d}{:s} of {:d}'.format(Fore.BLUE,j,Fore.RESET,self.tiles_num)
                print(desc,end='\r')
                df = pd.read_csv(tile_file,skiprows=1,dtype=str)
                # units: ra->hourangle, dec->deg, pmra->mas/a, pmdec->mas/a, epoch->2000.0 
                df_reduced = df.loc[:,['ra','dec','pmra','pmdec','mag']]
                df_reduced['epoch'] = '2000.0'
                df_reduced['ra'] = (df_reduced['ra'].astype(float)*15).round(6) # Convert hourangle to deg
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)
                j += 1
        elif sc_name == 'gsc12':
            j = 1
            for tile_file in file_list:
                desc = 'Reducing {:s}{:d}{:s} of {:d}'.format(Fore.BLUE,j,Fore.RESET,self.tiles_num)
                print(desc,end='\r')
                df = pd.read_csv(tile_file,skiprows=1,dtype=str)
                # units: ra->deg, dec->deg, RApm->arcs/0.1a, Decpm->arcs/0.1a, Epoch->1980.552
                df_reduced = df.loc[:,['ra','dec','RApm','Decpm','mag','Epoch']]
                columns_dict = {'RApm':'pmra', 'Decpm':'pmdec', 'Epoch':'epoch'}
                df_reduced.rename(columns=columns_dict, inplace=True)
                # Convert to standard proper motion in mas/a
                df_reduced['pmra'] = (df_reduced['pmra'].astype(float)*1e4).round(2) 
                df_reduced['pmdec'] = (df_reduced['pmdec'].astype(float)*1e4).round(2)
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)
                j += 1
        elif sc_name == 'gsc242':
            j = 1
            for tile_file in file_list:
                desc = 'Reducing {:s}{:d}{:s} of {:d}'.format(Fore.BLUE,j,Fore.RESET,self.tiles_num)
                print(desc,end='\r')
                df = pd.read_csv(tile_file,skiprows=1,dtype=str)
                # units: ra->deg, dec->deg, rapm->mas/a, decpm->mas/a, epoch->2012
                df_reduced = df.loc[:,['ra','dec','rapm','decpm','mag','epoch']]
                columns_dict = {'rapm':'pmra', 'decpm':'pmdec'}
                df_reduced.rename(columns=columns_dict, inplace=True)
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)  
                j += 1
        elif sc_name == 'gaiadr3':
            for tile_file in file_list:
                df = pd.read_csv(tile_file,skiprows=1,dtype=str)
                # units: ra->deg, dec->deg, pmra->mas/a, pmdec->mas/a, epoch->2016
                df_reduced = df.loc[:,['ra','dec','pmra','pmdec','mag','epoch']]
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)           
        elif sc_name == 'ucac5':
            j = 1
            for tile_file in file_list:
                desc = 'Reducing {:s}{:d}{:s} of {:d}'.format(Fore.BLUE,j,Fore.RESET,self.tiles_num)
                print(desc,end='\r')
                df = pd.read_csv(tile_file,skiprows=1,dtype=str)
                # units: ra->deg, dec->deg, pmur->mas/a, pmud->mas/a, epu->1998.754
                df_reduced = df.loc[:,['ra','dec','pmur','pmud','mag','epu']]
                columns_dict = {'pmur':'pmra', 'pmud':'pmdec','epu':'epoch'}
                df_reduced.rename(columns=columns_dict, inplace=True)
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)   
                j += 1 
        elif sc_name == 'usnob':
            j = 1
            for tile_file in file_list:
                desc = 'Reducing {:s}{:d}{:s} of {:d}'.format(Fore.BLUE,j,Fore.RESET,self.tiles_num)
                print(desc,end='\r')
                df = pd.read_csv(tile_file,skiprows=1,dtype=str)
                # units: ra->deg, dec->deg, pmRA->mas/a, pmDEC->mas/a, Epoch->1950,mag->unknown
                df_reduced = df.loc[:,['ra','dec','pmRA','pmDEC','mag','Epoch']]
                columns_dict = {'pmRA':'pmra', 'pmDEC':'pmdec','Epoch':'epoch'}
                df_reduced.rename(columns=columns_dict, inplace=True)
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)  
                j += 1
        elif sc_name == '2mass':
            j = 1
            for tile_file in file_list:
                desc = 'Reducing {:s}{:d}{:s} of {:d}'.format(Fore.BLUE,j,Fore.RESET,self.tiles_num)
                print(desc,end='\r')
                df = pd.read_csv(tile_file,skiprows=1,dtype=str)
                # units: ra->deg, dec->deg, jdate->2451063.6417
                df_reduced = df.loc[:,['ra','dec','mag']]
                df_reduced['epoch'] = Time(df['jdate'].astype(float), format='jd').jyear.round(2) 
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)  
                j += 1

        print('\nFinished')                                                   
                    
        file_num,dir_size,validity = tiles_statistic(dir_reduced,tile_size)  
        info['_mode'] = 'reduced'
        info['tiles_path'] = dir_reduced
        info['sc_size'] = dir_size

        return StarCatalogReduced(info)  
    
    def search_box(self,radec_box,mag_threshold,t_pm,max_num=None):
        """
        Perform a rectangle search of stars on raw star catalogs.

        Usage:
            >>> from starcatalogquery import StarCatalogRaw
            >>> dir_from_raw = '/Volumes/TOSHIBA/starcatalogs/raw/hygv35/res5/'
            >>> hygv35_raw = StarCatalogRaw.load('hygv35',5,dir_from_raw)
            >>> stars = hygv35_raw.search_box([20,30,30,40],9,2022.0)
        Inputs:
            radec_box -> [list] Rectangular search area in form of [ra_min,dec_min,ra_max,dec_max], where
                ra_min -> [float] Left border of RA in [deg].
                dec_min -> [float] Lower border of DEC in [deg].
                ra_max -> [float] Right border of RA in [deg].
                dec_max -> [float] Upper border of DEC in [deg].
            mag_threshold -> [float] Apparent magnitude limit  
            t_pm -> [float] Epoch to which the stars are unified
            max_num -> [int,optional,default=None] Maxinum number of the stars sorted by brightness for rectangle search. If None, all stars are counted in the search area by deault.
        Outputs:
            Instance of class Stars
        """
        ra_min,dec_min,ra_max,dec_max = radec_box
        tile_size = int(self.tile_size.split()[0])
        sc_path,sc_name = self.tiles_path,self.sc_name
        sc_indices = box2seqs(radec_box,tile_size) 

        df = pd.concat(_load_files(sc_indices,sc_path,sc_name,self._mode))

        if sc_name == 'hygv35':
            df['ra'] = (df['ra'].astype(float)*15).round(6) # Convert hourangle to deg
            df['epoch'] = 2000.0
        elif sc_name == 'gsc12':
            columns_dict = {'RApm':'pmra', 'Decpm':'pmdec', 'Epoch':'epoch'}
            df.rename(columns=columns_dict, inplace=True)
            # Convert to standard proper motion in mas/a
            df['pmra'] = (df['pmra'].astype(float)*1e4).round(2) 
            df['pmdec'] = (df['pmdec'].astype(float)*1e4).round(2)
        elif sc_name == 'gsc242':
            columns_dict = {'rapm':'pmra', 'decpm':'pmdec'}
            df.rename(columns=columns_dict, inplace=True)
        elif sc_name == 'gaiadr3':
            pass
        elif sc_name == 'ucac5':
            columns_dict = {'pmur':'pmra', 'pmud':'pmdec','epu':'epoch'}
            df.rename(columns=columns_dict, inplace=True) 
        elif sc_name == 'usnob':
            columns_dict = {'pmRA':'pmra', 'pmDEC':'pmdec','Epoch':'epoch'}
            df.rename(columns=columns_dict, inplace=True)
        elif sc_name == '2mass':
            df['epoch'] = Time(df['jdate'].astype(float), format='jd').jyear.round(2) 
                          
        df.drop_duplicates(subset=['ra','dec'],inplace=True)

        if {'pmra', 'pmdec'}.issubset(df.columns): 
            df[['ra','dec','pmra','pmdec','mag','epoch']] = df[['ra', 'dec','pmra','pmdec','mag','epoch']].apply(pd.to_numeric)
        else:    
            df[['ra','dec','mag','epoch']] = df[['ra','dec','mag','epoch']].apply(pd.to_numeric)

        mag_flag = (df['mag'] < mag_threshold)
        df = df[mag_flag].sort_values(by=['mag'])

        # calculate the proper motion
        dt = float(t_pm) - df['epoch']

        if {'pmra', 'pmdec'}.issubset(df.columns):    
            df['ra'] +=  df['pmra']/3.6e6 * dt   
            df['dec'] += df['pmdec']/3.6e6 * dt
        else:
            warnings.warn('Proper motion data for stars in catalog {:s} are not found.'.format(sc_name))

        ra,dec = df['ra'],df['dec']
        ra_flag = np.abs(ra - (ra_min + ra_max)/2) < (ra_max - ra_min)/2
        dec_flag = np.abs(dec- (dec_min + dec_max)/2) < (dec_max - dec_min)/2

        flag = ra_flag & dec_flag 
        df = df[flag]
        df['epoch'] = t_pm
        df.reset_index(drop=True,inplace=True)

        info = df2info(self.sc_name,center,df,max_num) 
        return Stars(info)

    def search_cone(self,center,radius,mag_threshold,t_pm,max_num=None):
        """
        Perform a cone search of stars on raw star catalogs.

        Usage:
            >>> from starcatalogquery import StarCatalogRaw
            >>> dir_from_raw = '/Volumes/TOSHIBA/starcatalogs/raw/hygv35/res5/'
            >>> hygv35_raw = StarCatalogRaw.load('hygv35',5,dir_from_raw)
            >>> stars = hygv35_raw.search_cone([20,30],10,9,2022.0)
        Inputs:
            center -> [list] Center of the cone in form of [ra_c,dec_c], where
                ra_c -> [float] RA, in [deg]
                dec_c -> [float] DEC, in [deg]
            radius -> [float] Angular radius of the cone, in [deg]
            mag_threshold -> [float] Apparent magnitude limit
            t_pm -> [float] Epoch to which the stars are unified
            max_num -> [int,optional,default=None] Maximum mumber of stars sorted by brightness for cone search. If None, all stars are counted in the search area by deault.
        Outputs:
            Instance of class Stars
        """
        ra_c,dec_c = center
        tile_size = int(self.tile_size.split()[0])
        sc_path,sc_name = self.tiles_path,self.sc_name
        sc_indices = cone2seqs(ra_c,dec_c,radius,tile_size) 

        df = pd.concat(_load_files(sc_indices,sc_path,sc_name,self._mode))

        if sc_name == 'hygv35':
            df['ra'] = (df['ra'].astype(float)*15).round(6) # Convert hourangle to deg
            df['epoch'] = 2000.0
        elif sc_name == 'gsc12':
            columns_dict = {'RApm':'pmra', 'Decpm':'pmdec', 'Epoch':'epoch'}
            df.rename(columns=columns_dict, inplace=True)
            # Convert to standard proper motion in mas/a
            df['pmra'] = (df['pmra'].astype(float)*1e4).round(2) 
            df['pmdec'] = (df['pmdec'].astype(float)*1e4).round(2)
        elif sc_name == 'gsc242':
            columns_dict = {'rapm':'pmra', 'decpm':'pmdec'}
            df.rename(columns=columns_dict, inplace=True)
        elif sc_name == 'gaiadr3':
            pass
        elif sc_name == 'ucac5':
            columns_dict = {'pmur':'pmra', 'pmud':'pmdec','epu':'epoch'}
            df.rename(columns=columns_dict, inplace=True) 
        elif sc_name == 'usnob':
            columns_dict = {'pmRA':'pmra', 'pmDEC':'pmdec','Epoch':'epoch'}
            df.rename(columns=columns_dict, inplace=True)
        elif sc_name == '2mass':
            df['epoch'] = Time(df['jdate'].astype(float), format='jd').jyear.round(2) 
                          
        df.drop_duplicates(subset=['ra','dec'],inplace=True)

        if {'pmra', 'pmdec'}.issubset(df.columns): 
            df[['ra','dec','pmra','pmdec','mag','epoch']] = df[['ra', 'dec','pmra','pmdec','mag','epoch']].apply(pd.to_numeric)
        else:    
            df[['ra','dec','mag','epoch']] = df[['ra','dec','mag','epoch']].apply(pd.to_numeric)

        mag_flag = (df['mag'] < mag_threshold)
        df = df[mag_flag].sort_values(by=['mag'])

        # calculate the proper motion
        dt = float(t_pm) - df['epoch']

        if {'pmra', 'pmdec'}.issubset(df.columns):    
            df['ra'] +=  df['pmra']/3.6e6 * dt   
            df['dec'] += df['pmdec']/3.6e6 * dt
        else:
            warnings.warn('Proper motion data for stars in catalog {:s} are not found.'.format(sc_name))

        ra,dec = df['ra'],df['dec']
        c1 = SkyCoord(ra,dec, unit='deg')
        c2 = SkyCoord(ra_c,dec_c, unit='deg')
        sep = c1.separation(c2).deg

        flag = sep < radius 
        df = df[flag]
        df['epoch'] = t_pm
        df.reset_index(drop=True,inplace=True)   

        info = df2info(self.sc_name,center,df,max_num) 
        return Stars(info)

    def _search_draw(self,search_area):
        """
        Visualize the scope of the search area and the coverage of the corresponding tiles.

        Usage:
            >>> from starcatalogquery import StarCatalogRaw
            >>> dir_from_raw = '/Volumes/TOSHIBA/starcatalogs/raw/hygv35/res5/'
            >>> hygv35_raw = StarCatalogRaw.load('hygv35',5,dir_from_raw)
            >>> cone_area = {'cone':[20,30,10]}
            >>> stars = hygv35_raw._search_draw(cone_area)
        Inputs:
            search_area -> [dict] Scope of the search area, such as {'cone':[20,30,10]} or {'box':[20,30,30,40]}, where
            for {'box':[ra_min,dec_min,ra_max,dec_max]}:
                ra_min -> [float] Left border of RA in [deg].
                dec_min -> [float] Lower border of DEC in [deg].
                ra_max -> [float] Right border of RA in [deg].
                dec_max -> [float] Upper border of DEC in [deg].
            for {'cone':[ra_c,dec_c,radius]}:
                ra_c -> [float] RA, in [deg].
                dec_c -> [float] DEC, in [deg].
                radius -> [float] Angular radius of the cap, in [deg].  
        Outputs:
            An image sketching the scope of the search area and the coverage of the corresponding tiles.          
        """
        tile_size = int(self.tile_size.split()[0]) 
        search_draw(tile_size,search_area)       

class StarCatalogReduced(object):
    """
    Class StarCatalogReduced

    Attributes:
        - tiles_path: Path of the star catalog files. 
        - tiles_num: Total number of the tile files.
        - tile_size: Geometric size of the tile in [deg]
        - sc_size: File size of the star catalog.
        - validity: Validity of the star catalog.
        - sc_name: Name of the star catalog. Available options include 'hygv35', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.
        - _mode: Type of the star catalog. Available options include 'raw', 'reduced', 'simplified', where
            'raw' represents the original star catalog, which covers all information about the stars,
            'reduced' represents the reduced star catalog, which contains the position, proper motion, apparent magnitude, epoch of the stars,
            'simplified' represents the minimalist star catalog, which only includes the position and apparent magnitude of the stars at a specific epoch.
        - stars_num: Total number of stars in the catalog
        - mag: Range of apparent magnitudes for the catalog
        - description: Catalog Summary
    Methods:
        - load: Load the reduced star catalog files from the local database.
        - simplify: Simplify the reduced star catalog so that it only contains key information of the stars, such as the celetial position and apparent magnitude, etc.
        - search_box: Perform a rectangle search on the reduced star catalog.
        - search_cone: Perform a cone search on the reduced star catalog.
        - _search_draw: Visualize the scope of the search area and the coverage of the corresponding tiles.
    """  
    def __init__(self,info):  

        self.info = info

        for key in info.keys():
            setattr(self, key, info[key])

    def __repr__(self):
    
        return 'Instance of class StarCatalogReduced'

    def load(sc_name,tile_size,dir_from=None):
        """
        Load the reduced star catalog files from the local database.

        Usage:
            >>> from starcatalogquery import StarCatalogReduced
            >>> # load the reduced star catalog HYGv3.5
            >>> dir_from_reduced = '/Volumes/TOSHIBA/starcatalogs/reduced/hygv35/res5/'
            >>> hygv35_reduced = StarCatalogReduced.load('hygv35',5,dir_from_reduced)
        Inputs:
            sc_name -> [str] Name of the star catalog. Available options include 'hygv35', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.
            tile_size -> [int] Geometric size of the tile in [deg]
            dir_from -> [str,optional,default=None] Directory of the star catalog files. If None, it is assigned to a build-in directory by default.
        Outputs:
            Instance of class StarCatalogReduced
        """
        if dir_from is None: dir_from = 'starcatalogs/reduced/{:s}/res{:d}/'.format(sc_name,tile_size)    
        if not os.path.exists(dir_from): raise Exception('Path of the star catalog {:s} does not exist.'.format(sc_name))  

        # calculate the total file size and numbers of tile files    
        file_num,dir_size,validity = tiles_statistic(dir_from,tile_size) 
        stars_num,mag,description = starcatalog_info(sc_name)

        dict_values = dir_from,dir_size,file_num,validity,sc_name,'{:d} deg'.format(tile_size),'reduced',stars_num,mag,description
        dict_keys = 'tiles_path','sc_size','tiles_num','validity','sc_name','tile_size','_mode','stars_num','mag','description'
        info = dict(zip(dict_keys, dict_values))

        return StarCatalogReduced(info)  

    def simplify(self,mag_threshold,t_pm,dir_simplified=None):
        """
        Simplify the reduced star catalog so that it only contains key information of the stars, such as the celetial position and apparent magnitude of stars.

        Usage:
            >>> from starcatalogquery import StarCatalogReduced
            >>> # load the reduced star catalog HYGv3.5
            >>> dir_from_reduced = '/Volumes/TOSHIBA/starcatalogs/reduced/hygv35/res5/'
            >>> hygv35_reduced = StarCatalogReduced.load('hygv35',5,dir_from_reduced)
            >>> hygv35_simplified = hygv35_reduced.simplify(9.0,2022.0)
        Inputs:
            mag_threshold -> [float] Apparent magnitude limit
            t_pm -> [float] Epoch to which the simplification is unified
            dir_simplified -> [str,optional,default=None] Directory of the simplified star catalog files. If None, it is assigned to a build-in directory by default.
        Outputs:
            Instance of class StarCatalogSimplified      
        """
        info = self.info.copy()
        tile_size = int(self.tile_size.split(' ')[0])
        tiles_path = self.tiles_path
        if dir_simplified is None:
            dir_simplified = tiles_path.replace('reduced','simplified')+'mag{:.1f}/epoch{:.1f}/'.format(mag_threshold,t_pm)
        Path(dir_simplified).mkdir(parents=True, exist_ok=True)    

        file_list = glob(tiles_path+'*') 
        sc_name = self.sc_name
        print('Simplifying the star catalog {:s}, which may take a considerable amount of time'.format(sc_name))  

        j = 1
        for tile_file in file_list:

            desc = 'Simplifying {:s}{:d}{:s} of {:d}'.format(Fore.BLUE,j,Fore.RESET,self.tiles_num)
            print(desc,end='\r')

            df_simplified = pd.read_csv(tile_file)
            mag = df_simplified['mag']
            mag_flag = mag < mag_threshold
            df_simplified = df_simplified[mag_flag].sort_values(by=['mag'])    

            # calculate proper motion    
            dt = float(t_pm) - df_simplified['epoch']  

            if {'pmra', 'pmdec'}.issubset(df_simplified.columns):    
                df_simplified['ra'] +=  df_simplified['pmra']/3.6e6 * dt   
                df_simplified['dec'] += df_simplified['pmdec']/3.6e6 * dt
                df_simplified = df_simplified.drop(['pmra','pmdec'],axis=1)
            else:
                warnings.warn('Proper motion data for stars in catalog {:s} are not found.'.format(sc_name))
                    
            # df_simplified['epoch'] = t_pm
            df_simplified.drop(columns=['epoch'],inplace=True)
            df_simplified['ra'] = df_simplified['ra'].round(6)
            df_simplified['dec'] = df_simplified['dec'].round(6)
            df_simplified['mag'] = df_simplified['mag'].round(1)

            df_simplified.to_csv(dir_simplified + tile_file.split('/')[-1],index=False)  

            j += 1

        print('\nFinished')

        file_num,dir_size,validity = tiles_statistic(dir_simplified,tile_size)  
        info['_mode'] = 'simplified'
        info['tiles_path'] = dir_simplified
        info['sc_size'] = dir_size
        info['mag_threshold'] = mag_threshold
        info['epoch'] = t_pm

        return StarCatalogSimplified(info)    

    def search_box(self,radec_box,mag_threshold,t_pm,max_num=None):
        """
        Perform a rectangle search of stars on the reduced star catalog.

        Usage:
            >>> from starcatalogquery import StarCatalogReduced
            >>> # load the reduced star catalog HYGv3.5
            >>> dir_from_reduced = '/Volumes/TOSHIBA/starcatalogs/reduced/hygv35/res5/'
            >>> hygv35_reduced = StarCatalogReduced.load('hygv35',5,dir_from_reduced)
            >>> stars = hygv35_reduced.search_box([20,30,30,40],9,2022.0)
        Inputs:
            radec_box -> [list] Rectangular search area in form of [ra_min,dec_min,ra_max,dec_max], where
                ra_min -> [float] Left border of RA in [deg]
                dec_min -> [float] Lower border of DEC in [deg]
                ra_max -> [float] Right border of RA in [deg]
                dec_max -> [float] Upper border of DEC in [deg]
            mag_threshold -> [float] Apparent magnitude limit
            t_pm -> [float] Epoch to which the stars are unified
            max_num -> [int,optional,default=None] Maxinum number of the stars sorted by brightness for rectangle search. If None, all stars are counted in the search area by deault.
        Outputs:
            Instance of class Stars
        """
        width = int(self.tile_size.split()[0])
        df = search_box_magpm(radec_box,self.tiles_path,self.sc_name,width,self._mode,mag_threshold,t_pm)
        info = df2info(self.sc_name,center,df,max_num) 
        return Stars(info)

    def search_cone(self,center,radius,mag_threshold,t_pm,max_num=None):   
        """
        Perform a cone search of stars on the reduced star catalog.

        Usage:
            >>> from starcatalogquery import StarCatalogReduced
            >>> # load the reduced star catalog HYGv3.5
            >>> dir_from_reduced = '/Volumes/TOSHIBA/starcatalogs/reduced/hygv35/res5/'
            >>> hygv35_reduced = StarCatalogReduced.load('hygv35',5,dir_from_reduced)
            >>> stars = hygv35_reduced.search_cone([20,30],10,9,2022.0)
        Inputs:
            center -> [list] Center of the cone in form of [ra_c,dec_c], where
                ra_c -> [float] RA, in [deg].
                dec_c -> [float] DEC, in [deg].
            radius -> [float] Angular radius of the cone, in [deg].
            mag_threshold -> [float] Apparent magnitude limit  
            t_pm -> [float] Epoch to which the stars are unified
            max_num -> [int,optional,default=None] Maxinum number of the stars sorted by brightness for cone search. If None, all stars are counted in the search area by deault.
        Outputs:
            Instance of class Stars
        """
        width = int(self.tile_size.split()[0])
        df = search_cone_magpm(center,radius,self.tiles_path,self.sc_name,width,self._mode,mag_threshold,t_pm)
        info = df2info(self.sc_name,center,df,max_num) 
        return Stars(info)

    def _search_draw(self,search_area):
        """
        Visualize the scope of the search area and the coverage of the corresponding tiles.

        Usage:
            >>> from starcatalogquery import StarCatalogReduced
            >>> # load the reduced star catalog HYGv3.5
            >>> dir_from_reduced = '/Volumes/TOSHIBA/starcatalogs/reduced/hygv35/res5/'
            >>> hygv35_reduced = StarCatalogReduced.load('hygv35',5,dir_from_reduced)
            >>> cone_area = {'cone':[20,30,10]}
            >>> stars = hygv35_reduced._search_draw(cone_area)
        Inputs:
            search_area -> [dict] Scope of the search area, such as {'cone':[20,30,10]} or {'box':[20,30,30,40]}, where
            for {'box':[ra_min,dec_min,ra_max,dec_max]}:
                ra_min -> [float] Left border of RA in [deg]
                dec_min -> [float] Lower border of DEC in [deg]
                ra_max -> [float] Right border of RA in [deg]
                dec_max -> [float] Upper border of DEC in [deg]
            for {'cone':[ra_c,dec_c,radius]}:
                ra_c -> [float] RA, in [deg]
                dec_c -> [float] DEC, in [deg]
                radius -> [float] Angular radius of the cone, in [deg] 
        Outputs:
            An image sketching the scope of the search area and the coverage of the corresponding tiles.          
        """
        width = int(self.tile_size.split()[0]) 
        search_draw(width,search_area)   

class StarCatalogSimplified(object):
    """
    Class StarCatalogSimplified

    Attributes:
        - tiles_path: Path of the star catalog files. 
        - tiles_num: Total number of the tile files.
        - tile_size: Geometric size of the tile in [deg]
        - sc_size: File size of the star catalog.
        - validity: Validity of the star catalog.
        - sc_name: Name of the star catalog. Available options include 'hygv35', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.
        - _mode: Type of the star catalog. Available options include 'raw', 'reduced', 'simplified', where
            'raw' represents the original star catalog, which covers all information about the stars,
            'reduced' represents the reduced star catalog, which contains the position, proper motion, apparent magnitude, epoch of the stars,
            'simplified' represents the minimalist star catalog, which only includes the position and apparent magnitude of the stars at a specific epoch.
        - stars_num: Total number of stars in the catalog
        - mag: Range of apparent magnitudes for the catalog
        - description: Catalog Summary
    Methods:
        - load: Load the simplified star catalog files from the local database.
        - search_box: Perform a rectangle search on the simplified star catalog.
        - search_cone: Perform a cone search on the simplified star catalog.
        - _search_draw: Visualize the scope of the search area and the coverage of the corresponding tiles.
        - h5_incices: Generate a h5-formatted star catalog indices file, which records the center pointing, pixel coordinates of the stars, triangle invariants and asterism indices of each sky area.
    """    
    def __init__(self,info):  

        self.info = info

        for key in info.keys():
            setattr(self, key, info[key])

    def __repr__(self):
    
        return 'Instance of class StarCatalogSimplified'        

    def load(sc_name,tile_size,mag_threshold,epoch,dir_from=None):
        """
        Load the simplified star catalog files from the local database.

        Usage:
            >>> from starcatalogquery import StarCatalogSimplified
            >>> # load the simplified star catalog HYGv3.5
            >>> dir_from_simplified = '/Volumes/TOSHIBA/starcatalogs/simplified/hygv35/res5/mag9.0/epoch2022.0/'
            >>> hygv35_simplified = StarCatalogSimplified.load('hygv35',5,9,2022,dir_from_simplified)
        Inputs:
            sc_name -> [str] Name of the star catalog. Available options include 'hygv35', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.
            tile_size -> [int] Geometric size of the tile in [deg]
            mag_threshold -> [float] Apparent magnitude limit
            epoch -> [float] Epoch of the star catalog
            dir_from -> [str,optional,default=None] Diectory of the star catalog files. If None, it is assigned to a build-in directory by default.
        Outputs:
            Instance of class StarCatalogSimplified
        """ 
        if dir_from is None: dir_from = 'starcatalogs/simplified/{:s}/res{:d}/mag{:.1f}/epoch{:.1f}/'.format(sc_name,tile_size,mag_threshold,epoch)
        if not os.path.exists(dir_from): raise Exception('Path of the star catalog {:s} does not exist.'.format(sc_name))  

        # calculate total size and numbers of tile files    
        file_num,dir_size,validity = tiles_statistic(dir_from,tile_size) 
        stars_num,mag,description = starcatalog_info(sc_name)

        dict_values = dir_from,dir_size,file_num,validity,sc_name,'{:d} deg'.format(tile_size),'simplified',stars_num,mag,description,mag_threshold,epoch
        dict_keys = 'tiles_path','sc_size','tiles_num','validity','sc_name','tile_size','_mode','stars_num','mag','description','mag_threshold','epoch'
        info = dict(zip(dict_keys, dict_values))

        return StarCatalogSimplified(info)   

    def search_box(self,radec_box,max_num=None):
        """
        Perform a rectangle search of stars on the simplified star catalog.

        Usage:
            >>> from starcatalogquery import StarCatalogSimplified
            >>> # load the simplified star catalog HYGv3.5
            >>> dir_from_simplified = '/Volumes/TOSHIBA/starcatalogs/simplified/hygv35/res5/mag9.0/epoch2022.0/'
            >>> hygv35_simplified = StarCatalogSimplified.load('hygv35',5,9,2022,dir_from_simplified)
            >>> stars = hygv35_simplified.search_box([20,30,30,40])
        Inputs:
            radec_box -> [list] Rectangular search area in form of [ra_min,dec_min,ra_max,dec_max], where
                ra_min -> [float] Left border of RA in [deg]
                dec_min -> [float] Lower border of DEC in [deg]
                ra_max -> [float] Right border of RA in [deg]
                dec_max -> [float] Upper border of DEC in [deg]
            max_num -> [int,optional,default=None] Maxinum number of the stars sorted by brightness for rectangle search. If None, all stars are counted in the search area by deault.
        Outputs:
            Instance of class Stars
        """
        width = int(self.tile_size.split()[0])
        df = search_box(radec_box,self.tiles_path,self.sc_name,width,self._mode) 
        info = df2info(self.sc_name,center,df,max_num) 
        return Stars(info)

    def search_cone(self,center,radius,max_num=None):   
        """
        Perform a cone search of stars on the simplified star catalog.

        Usage:
            >>> from starcatalogquery import StarCatalogSimplified
            >>> # load the simplified star catalog HYGv3.5
            >>> dir_from_simplified = '/Volumes/TOSHIBA/starcatalogs/simplified/hygv35/res5/mag9.0/epoch2022.0/'
            >>> hygv35_simplified = StarCatalogSimplified.load('hygv35',5,9,2022,dir_from_simplified)
            >>> stars = hygv35_simplified.search_cone([20,30],10)
        Inputs:
            center -> [int,array_like] Center of the cone in form of [ra_c,dec_c], where
                ra_c -> [float] RA, in [deg]
                dec_c -> [float] DEC, in [deg]
            radius -> [float] Angular radius of the cone, in [deg].
            max_num -> [int,optional,default=None] Maxinum number of the stars sorted by brightness for rectangle search. If None, all stars are counted in the search area by deault.
        Outputs:
            Instance of class Stars
        """
        width = int(self.tile_size.split()[0])
        df = search_cone(center,radius,self.tiles_path,self.sc_name,width,self._mode)
        info = df2info(self.sc_name,center,df,max_num) 
        return Stars(info)

    def _search_draw(self,search_area):
        """
        Visualize the scope of the search area and the coverage of the corresponding tiles.

        Usage:
            >>> from starcatalogquery import StarCatalogSimplified
            >>> # load the simplified star catalog HYGv3.5
            >>> dir_from_simplified = '/Volumes/TOSHIBA/starcatalogs/simplified/hygv35/res5/mag9.0/epoch2022.0/'
            >>> hygv35_simplified = StarCatalogReduced.load('hygv35',5,9,2022,dir_from_simplified)
            >>> cone_area = {'cone':[20,30,10]}
            >>> stars = hygv35_simplified._search_draw(cone_area)
        Inputs:
            search_area -> [dict] Scope of the search area, such as {'cone':[20,30,10]} or {'box':[20,30,30,40]}, where
            for {'box':[ra_min,dec_min,ra_max,dec_max]}:
                ra_min -> [float] Left border of RA in [deg]
                dec_min -> [float] Lower border of DEC in [deg]
                ra_max -> [float] Right border of RA in [deg]
                dec_max -> [float] Upper border of DEC in [deg]
            for {'cone':[ra_c,dec_c,radius]}:
                ra_c -> [float] RA, in [deg]
                dec_c -> [float] DEC, in [deg]
                radius -> [float] Angular radius of the cone, in [deg].  
        Outputs:
            An image sketching the scope of the search area and the coverage of the corresponding tiles.          
        """
        width = int(self.tile_size.split()[0]) 
        search_draw(width,search_area)   

    def h5_incices(self,fov,pixel_width,max_num=30):
        """
        Generate a h5-formatted star catalog incices file, which records the center pointing, pixel coordinates of the stars, triangle invariants and asterism indices of each sky area.

        Usage:
            >>> from starcatalogquery import StarCatalogSimplified
            >>> # load the simplified star catalog HYGv3.5
            >>> dir_from_simplified = '/Volumes/TOSHIBA/starcatalogs/simplified/hygv35/res5/mag9.0/epoch2022.0/'
            >>> hygv35_simplified = StarCatalogSimplified.load('hygv35',5,9,2022,dir_from_simplified)
            >>> outh5 = hygv35_simplified.h5_incices(8,0.01)
        Inputs:
            fov -> [float] FOV of camera
            pixel_width -> [float] Pixel width in [deg]
            max_num -> [int,optional,default=30] Maxinum number of the stars sorted by brightness used for a sky area.
        Outputs:
            outh5 -> h5-formatted star catalog indices file    
        """ 
        dir_h5 = 'starcatalogs/indices/{:s}/'.format(self.sc_name)
        Path(dir_h5).mkdir(parents=True, exist_ok=True)  

        if fov >= 28.4 and fov < 58.6: 
            k,nside = 0,1
            search_radius = 37.2
        elif fov >= 13.4 and fov < 28.4:
            k,nside = 1,2
            search_radius = 17.0
        elif fov >= 6.6 and fov < 13.4:
            k,nside = 2,4
            search_radius = 9.3
        elif fov >= 3.3 and fov < 6.6:
            k,nside = 3,8
            search_radius = 4.1  
        else:
            raise Exception('FOV should be between 3.3 and 58.6 deg')     

        outh5 = dir_h5 + 'k{:d}_mag{:.1f}_mcp{:d}_{:.1f}.h5'.format(k,self.mag_threshold,max_num,self.epoch)

        if os.path.exists(outh5): return outh5 

        # write to h5 file
        fout = h5py.File(outh5,'w')
        # create group 
        stars_xy_grp = fout.create_group("stars_xy")
        stars_invariants_grp = fout.create_group("stars_invariants")
        stars_asterisms_grp = fout.create_group("stars_asterisms")

        npix = hp.nside2npix(nside)
        fp_radecs = []
        for seq in range(npix):

            desc = 'Generating starcatalog sky area index {:s}{:d}{:s} of {:d} for level k{:d}'.format(Fore.BLUE,seq+1,Fore.RESET,npix,k)
            print(desc,end='\r')

            fp_radec = hp.pix2ang(nside,seq,lonlat=True)
            fp_radecs.append(fp_radec)

            stars = self.search_cone(fp_radec,search_radius,max_num)
            stars.pixel_xy(pixel_width)
            stars.invariantfeatures()

            stars_xy_grp.create_dataset(str(seq), data=stars.xy)
            stars_invariants_grp.create_dataset(str(seq), data=stars.invariants)
            stars_asterisms_grp.create_dataset(str(seq), data=stars.asterisms)  

        fout.create_dataset("fp_radecs", data=np.array(fp_radecs))
        fout.close() # close file
        return outh5    
            
class Stars(object):
    """
    Class Stars

    Attributes:
        - sc_name: Name of the star catalog
        - center: Center pointing in form of [ra_c,dec_c] in [deg]
        - df: Pandas dataframe of the stars
        - max_num: Number of stars in the dataframe
        - xy: Pixel coordinates of stars
        - radec: Celestial coordinates of stars
    Methods:
        - pixel_xy: Calculate the pixel coordinates of stars in a sky area.
    """    

    def __init__(self,info):  

        self.info = info

        for key in info.keys():
            setattr(self, key, info[key])

    def __repr__(self):
    
        return 'Instance of class Stars'

    def pixel_xy(self,pixel_width):
        """
        Calculate the pixel coordinates of stars in a sky area.

        Usage:
            >>> from starcatalogquery import StarCatalogSimplified
            >>> # load the simplified star catalog HYGv3.5
            >>> dir_from_simplified = '/Volumes/TOSHIBA/starcatalogs/simplified/hygv35/res5/mag9.0/epoch2022.0/'
            >>> hygv35_simplified = StarCatalogReduced.load(hygv35',5,9,2022,dir_from_simplified)
            >>> hygv35_simplified_stars = hygv35_simplified.search_cone([20,30],10)
            >>> stars = hygv35_simplified_stars.pixel_xy(0.01)
        Inputs:
            pixel_width -> [float] Pixel width in [deg]
        Outputs:
            Instance of class Stars       
        """
        # Define WCS transformation and convert the celestial coordinates to pixel coordinates
        x,y,wcs = xy_catalog(self.center,self.radec,pixel_width)

        df = self.info['df']
        df['pixelx'],df['pixely'] = x,y
        xy = np.stack([x,y]).T
        self.xy = self.info['xy'] = xy
        self.wcs = self.info['wcs'] = wcs
        return self

    def invariantfeatures(self):
        """
        1. Calculate the unique invariants (L2/L1,L1/L0), where L2 >= L1 >= L0 are the three sides of the triangle composed of stars.
        2. Construct the 2D Tree from the the unique invariants.
        3. Record an array of the indices of stars that correspond to each invariant.
        """
        if not hasattr(self,'xy'): raise Exception("The pixel coordinates of stars should be caiculated first by `.pixel_xy(pixel_width`)")
        inv_uniq, triang_vrtx_uniq = _generate_invariants(self.xy)
        inv_uniq_tree = KDTree(inv_uniq)
        self.info.update({'invariants':inv_uniq,'asterisms':triang_vrtx_uniq,'invariants_2dtree':inv_uniq_tree})
        self.invariants,self.asterisms,self.kdtree = inv_uniq,triang_vrtx_uniq,inv_uniq_tree
        return self  