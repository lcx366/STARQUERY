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

from .catalog_download import catalog_download,hygv3_download
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
        Grab star catalog data files from a remote server.

        Usage:
            >>> from starcatalogquery import StarCatalog
            >>> gaiadr3_raw = StarCatalog.get('gaiadr3',2)

        Inputs:
            sc_name -> [str] Name of the starcatalog to download. Available starcatalogs include 'hygv3', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.
            tile_size -> [int,optinal,default=None] size of the tile in [deg]. If None, the size of the tile is automatically assigned a feasible maximum according to the name of the star catalog.
            dir_to -> [str,optional,default=None] The download path of the star catalog files. If None, the path is automatically assigned to a suitable directory by default.

        Outputs:
            Instance of class StarCatalogRaw
        """
        if sc_name == 'hygv3': # for HYG v3 star catalog
            dir_to,dir_size,file_num,validity,tile_size = hygv3_download(tile_size,dir_to) 

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
            >>> # load the raw star catalog GAIADR3
            >>> dir_from_raw = '/Volumes/TOSHIBA/starcatalog/raw/gaiadr3/res2/'
            >>> gaiadr3_raw = StarCatalog.load(dir_from_raw)
            >>>
            >>> # load the reduced star catalog GAIADR3
            >>> dir_from_reduced = '/Volumes/TOSHIBA/starcatalog/reduced/gaiadr3/res2/'
            >>> gaiadr3_reduced = StarCatalog.load(dir_from_reduced)
            >>>
            >>> # load the simplified star catalog GAIADR3
            >>> dir_from_simplified = '/Volumes/TOSHIBA/starcatalog/simplified/gaiadr3/res2/mag8.0/epoch2023.0/'
            >>> gaiadr3_simplified = StarCatalog.load(dir_from_simplified)
        Inputs:
            dir_from -> [str,optional,default=None] The loading path of the star catalog files. If None, the path is automatically assigned to a suitable directory by default.
        Outputs:
            Instance of class StarCatalog
        """
        _mode,sc_name,tile_size = dir_from.split('starcatalogs/')[1].split('/')[:3]
        tile_size = int(tile_size[3:])

        # _mode -> [str] Types of star catalogs, including 'raw', 'reduced', 'simplified', where
        # 'raw' represents the original star catalog, which contains all information about the star
        # 'reduced' represents the reduced star catalog, which contains the position, proper motion, apparent magnitude, epoch of the star
        # 'simplified' represents the minimalist star catalog, which only includes the position and apparent magnitude of stars at a specific epoch
        # sc_name -> [str] Name of the starcatalog. Available starcatalogs include 'hygv3', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.
        # tile_size -> [int] Size of the tile in [deg]

        if _mode == 'raw':
            starcatalog = StarCatalogRaw.load(sc_name,tile_size,dir_from)
        elif _mode == 'reduced':
            starcatalog = StarCatalogReduced.load(sc_name,tile_size,dir_from)
        elif _mode == 'simplified':    
            mag_cutoff,epoch = np.array(dir_from.split('mag')[1][:-1].split('/epoch'),dtype=float) 
            starcatalog = StarCatalogSimplified.load(sc_name,tile_size,mag_cutoff,epoch,dir_from)

        return starcatalog 

    def read_h5_indices(infile):
        """
        Read in h5-formatted star catalog file, which records the center pointing of each sky area, the pixel coordinates of the stars, the triangle invariants and the asterism indices.

        Usage:
            >>> from starcatalogquery import StarCatalog
            >>> infile_h5 = 'starcatalogs/indices/hygv3/fov20_mag8_mcp40_2023.0.h5'
            >>> fp_radecs,stars_xy,stars_invariants,stars_asterisms = read_h5_indices(infile_h5)
        
        Inputs:
            infile_h5 -> [str] h5-formatted star catalog file  

        Outputs:
            fp_radecs -> [2d array of float] The center pointing of each sky area in format of [[RA0,DEC0],..[RAn,DECn]] in [deg]
            stars_xy -> [list of 2d array of float] The pixel coordinates of the stars in each sky area
            stars_invariants -> [list of 2d array of float] The triangle invariants in each sky area
            stars_asterisms -> [list of 2d array of int] The asterism indices corresponding to the triangle invariants in each sky area       
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
        - tiles_path: Path of the starcatalog tile files. 
        - sc_size: The size of the star catalog.
        - tiles_num: Total number of the tile files.
        - validity: The validity of the star catalog.
        - sc_name: Name of the starcatalog. Available starcatalogs include 'hygv3', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.
        - tile_size: Size of the tile in [deg]
        - _mode: Types of star catalogs, including 'raw', 'reduced', 'simplified', where
            'raw' represents the original star catalog, which contains all information about the star
            'reduced' represents the reduced star catalog, which contains the position, proper motion, apparent magnitude, epoch of the star
            'simplified' represents the minimalist star catalog, which only includes the position and apparent magnitude of stars
        - stars_num: Total number of stars in the catalog
        - mag: Apparent magnitude of stars in the catalog
        - description: Catalog Summary

    Methods:
        - load: Load the raw star catalog files from the local database.
        - reduce: Reduce the original star catalog so that the reduced star catalog only contains necessary information such as the position, proper motion, apparent magnitude, epoch, etc.
        - search_box: Perform a rectangle search of stars on the raw star catalog and return an instance of class Stars.
        - search_cone: Perform a cone search of stars on the raw star catalog and return an instance of class Stars.   
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
            >>> # load the raw star catalog GAIADR3
            >>> dir_from_raw = '/Volumes/TOSHIBA/starcatalog/raw/gaiadr3/res2/'
            >>> gaiadr3_raw = StarCatalogRaw.load('raw','gaiadr3',2,dir_from_raw)

        Inputs:
            sc_name -> [str] Name of the starcatalog. Available starcatalogs include 'hygv3', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.
            tile_size -> [int] Size of the tile in [deg]
            dir_from -> [str,optional,default=None] The loading path of the star catalog files. If None, the path is automatically assigned to a suitable directory by default.

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
        Reduce the original star catalog so that the reduced star catalog only contains necessary information such as the position, proper motion, apparent magnitude, epoch, etc.

        Usage:
            >>> from starcatalogquery import StarCatalogRaw
            >>> dir_from_raw = '/Volumes/TOSHIBA/starcatalog/raw/gaiadr3/res2/'
            >>> gaiadr3_raw = StarCatalogRaw.load('raw','gaiadr3',2,dir_from_raw)
            >>> gaiadr3_reduced = gaiadr3_raw.reduce()

        Inputs:
            dir_reduced -> [str,optional,default=None] The path of the reduced star catalog files to store. If None, the path is automatically assigned to a suitable directory by default.

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

        if sc_name == 'hygv3':
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
    
    def search_box(self,radec_box,mag_threshold,t_pm,max_control_points=None):
        """
        Perform a rectangle search of stars on the raw star catalog and return an instance of class Stars.

        Usage:
            >>> from starcatalogquery import StarCatalogRaw
            >>> dir_from_raw = '/Volumes/TOSHIBA/starcatalog/raw/gaiadr3/res2/'
            >>> gaiadr3_raw = StarCatalogRaw.load('raw','gaiadr3',2,dir_from_raw)
            >>> stars = gaiadr3_raw.search_box([20,30,30,40],8,2023.0)

        Inputs:
            radec_box -> [int,array_like] Rectangular search area in format of [ra_min,dec_min,ra_max,dec_max], where
                ra_min -> [float] Left border of RA in [deg].
                dec_min -> [float] Lower border of DEC in [deg].
                ra_max -> [float] Right border of RA in [deg].
                dec_max -> [float] Upper border of DEC in [deg].
            mag_threshold -> [float] Apparent magnitude limit of the detector  
            t_pm -> [float] The epoch when the search was performed
            max_control_points -> [int,optional,default=None] Number of brightest stars in the search area. If None, it is the number of all stars in the search area by deault.

        Outputs:
            Instance of class Stars
        """
        ra_min,dec_min,ra_max,dec_max = radec_box
        tile_size = int(self.tile_size.split()[0])
        sc_path,sc_name = self.tiles_path,self.sc_name
        sc_indices = box2seqs(radec_box,tile_size) 

        df = pd.concat(_load_files(sc_indices,sc_path,sc_name,self._mode))

        if sc_name == 'hygv3':
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

        # calculate proper motion
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

        info = df2info(self.sc_name,center,df,max_control_points) 
        return Stars(info)

    def search_cone(self,center,radius,mag_threshold,t_pm,max_control_points=None):
        """
        Perform a cone search of stars on the raw star catalog and return an instance of class Stars.

        Usage:
            >>> from starcatalogquery import StarCatalogRaw
            >>> dir_from_raw = '/Volumes/TOSHIBA/starcatalog/raw/gaiadr3/res2/'
            >>> gaiadr3_raw = StarCatalogRaw.load('raw','gaiadr3',2,dir_from_raw)
            >>> stars = gaiadr3_raw.search_cone([20,30],10,8,2023.0)

        Inputs:
            center -> [int,array_like] Center of the cap in format of [ra_c,dec_c], where
                ra_c -> [float] RA, in [deg].
                dec_c -> [float] DEC, in [deg].
            radius -> [float] Angular radius of the cap, in [deg].
            mag_threshold -> [float] Apparent magnitude limit of the detector  
            t_pm -> [float] The epoch when the search was performed
            max_control_points -> [int,optional,default=None] Number of brightest stars in the search area. If None, it is the number of all stars in the search area by deault.

        Outputs:
            Instance of class Stars
        """
        ra_c,dec_c = center
        tile_size = int(self.tile_size.split()[0])
        sc_path,sc_name = self.tiles_path,self.sc_name
        sc_indices = cone2seqs(ra_c,dec_c,radius,tile_size) 

        df = pd.concat(_load_files(sc_indices,sc_path,sc_name,self._mode))

        if sc_name == 'hygv3':
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

        # calculate proper motion
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

        info = df2info(self.sc_name,center,df,max_control_points) 
        return Stars(info)

    def _search_draw(self,search_area):
        """
        Visualize the scope of the search area and the coverage of the corresponding tiles.

        Usage:
            >>> from starcatalogquery import StarCatalogRaw
            >>> dir_from_raw = '/Volumes/TOSHIBA/starcatalog/raw/gaiadr3/res2/'
            >>> gaiadr3_raw = StarCatalogRaw.load('raw','gaiadr3',2,dir_from_raw)
            >>> cone_area = {'cone':[20,30,10]}
            >>> stars = gaiadr3_raw._search_draw(cone_area)

        Inputs:
            search_area -> [dict] The scope of the search area, such as {'cone':[20,30,10]} or {'box':[20,30,30,40]}, where
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
            An image that shows the scope of the search area and the coverage of the corresponding tiles.          
        """
        tile_size = int(self.tile_size.split()[0]) 
        search_draw(tile_size,search_area)       

class StarCatalogReduced(object):
    """
    Class StarCatalogReduced

    Attributes:
        - tiles_path: Path of the starcatalog tile files. 
        - sc_size: The size of the star catalog.
        - tiles_num: Total number of the tile files.
        - validity: The validity of the star catalog.
        - sc_name: Name of the starcatalog. Available starcatalogs include 'hygv3', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.
        - tile_size: Size of the tile in [deg]
        - _mode: Types of star catalogs, including 'raw', 'reduced', 'simplified', where
            'raw' represents the original star catalog, which contains all information about the star
            'reduced' represents the reduced star catalog, which contains the position, proper motion, apparent magnitude, epoch of the star
            'simplified' represents the minimalist star catalog, which only includes the position and apparent magnitude of stars
        - stars_num: Total number of stars in the catalog
        - mag: Apparent magnitude of stars in the catalog
        - description: Catalog Summary

    Methods:
        - load: Load the reduced star catalog files from the local database.
        - simplify: Simplify the reduced star catalog so that the simplified star catalog only contains key information: the position and apparent magnitude of stars.
        - search_box: Perform a rectangle search of stars on the reduced star catalog and return an instance of class Stars.
        - search_cone: Perform a cone search of stars on the reduced star catalog and return an instance of class Stars.   
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
            >>> # load the reduced star catalog GAIADR3
            >>> dir_from_reduced = '/Volumes/TOSHIBA/starcatalog/reduced/gaiadr3/res2/'
            >>> gaiadr3_reduced = StarCatalogReduced.load('reduced','gaiadr3',2,dir_from_reduced)

        Inputs:
            sc_name -> [str] Name of the starcatalog. Available starcatalogs include 'hygv3', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.
            tile_size -> [int] Size of the tile in [deg]
            dir_from -> [str,optional,default=None] The loading path of the star catalog files. If None, the path is automatically assigned to a suitable directory by default.

        Outputs:
            Instance of class StarCatalogReduced
        """
        if dir_from is None: dir_from = 'starcatalogs/reduced/{:s}/res{:d}/'.format(sc_name,tile_size)    
        if not os.path.exists(dir_from): raise Exception('Path of the star catalog {:s} does not exist.'.format(sc_name))  

        # calculate total size and numbers of tile files    
        file_num,dir_size,validity = tiles_statistic(dir_from,tile_size) 
        stars_num,mag,description = starcatalog_info(sc_name)

        dict_values = dir_from,dir_size,file_num,validity,sc_name,'{:d} deg'.format(tile_size),'reduced',stars_num,mag,description
        dict_keys = 'tiles_path','sc_size','tiles_num','validity','sc_name','tile_size','_mode','stars_num','mag','description'
        info = dict(zip(dict_keys, dict_values))

        return StarCatalogReduced(info)  

    def simplify(self,mag_threshold,t_pm,dir_simplified=None):
        """
        Simplify the reduced star catalog so that the simplified star catalog only contains key information: the position and apparent magnitude of stars.

        Usage:
            >>> from starcatalogquery import StarCatalogReduced
            >>> # load the reduced star catalog GAIADR3
            >>> dir_from_reduced = '/Volumes/TOSHIBA/starcatalog/reduced/gaiadr3/res2/'
            >>> gaiadr3_reduced = StarCatalogReduced.load('reduced','gaiadr3',2,dir_from_reduced)
            >>> gaiadr3_simplified = gaiadr3_reduced.simplify()

        Inputs:
            mag_threshold -> [float] Apparent magnitude limit of the detector  
            t_pm -> [float] The epoch when the search was performed
            dir_simplified -> [str,optional,default=None] The path of the simplified star catalog files to store. If None, the path is automatically assigned to a suitable directory by default.

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
                    
            #df_simplified['epoch'] = t_pm
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
        info['mag_cutoff'] = mag_threshold
        info['epoch'] = t_pm

        return StarCatalogSimplified(info)    

    def search_box(self,radec_box,mag_threshold,t_pm,max_control_points=None):
        """
        Perform a rectangle search of stars on the reduced star catalog and return an instance of class Stars.

        Usage:
            >>> from starcatalogquery import StarCatalogReduced
            >>> # load the reduced star catalog GAIADR3
            >>> dir_from_reduced = '/Volumes/TOSHIBA/starcatalog/reduced/gaiadr3/res2/'
            >>> gaiadr3_reduced = StarCatalogReduced.load('reduced','gaiadr3',2,dir_from_reduced)
            >>> stars = gaiadr3_reduced.search_box([20,30,30,40],8,2023.0)

        Inputs:
            radec_box -> [int,array_like] Rectangular search area in format of [ra_min,dec_min,ra_max,dec_max], where
                ra_min -> [float] Left border of RA in [deg].
                dec_min -> [float] Lower border of DEC in [deg].
                ra_max -> [float] Right border of RA in [deg].
                dec_max -> [float] Upper border of DEC in [deg].
            mag_threshold -> [float] Apparent magnitude limit of the detector  
            t_pm -> [float] The epoch when the search was performed    
            max_control_points -> [int,optional,default=None] Number of brightest stars in the search area. If None, it is the number of all stars in the search area by deault.

        Outputs:
            Instance of class Stars
        """
        width = int(self.tile_size.split()[0])
        df = search_box_magpm(radec_box,self.tiles_path,self.sc_name,width,self._mode,mag_threshold,t_pm)
        info = df2info(self.sc_name,center,df,max_control_points) 
        return Stars(info)

    def search_cone(self,center,radius,mag_threshold,t_pm,max_control_points=None):   
        """
        Perform a cone search of stars on the reduced star catalog and return an instance of class Stars.

        Usage:
            >>> from starcatalogquery import StarCatalogReduced
            >>> # load the reduced star catalog GAIADR3
            >>> dir_from_reduced = '/Volumes/TOSHIBA/starcatalog/reduced/gaiadr3/res2/'
            >>> gaiadr3_reduced = StarCatalogReduced.load('reduced','gaiadr3',2,dir_from_reduced)
            >>> stars = gaiadr3_reduced.search_cone([20,30],10,8,2023.0)

        Inputs:
            center -> [int,array_like] Center of the cap in format of [ra_c,dec_c], where
                ra_c -> [float] RA, in [deg].
                dec_c -> [float] DEC, in [deg].
            radius -> [float] Angular radius of the cap, in [deg].
            mag_threshold -> [float] Apparent magnitude limit of the detector  
            t_pm -> [float] The epoch when the search was performed
            max_control_points -> [int,optional,default=None] Number of brightest stars in the search area. If None, it is the number of all stars in the search area by deault.

        Outputs:
            Instance of class Stars
        """
        width = int(self.tile_size.split()[0])
        df = search_cone_magpm(center,radius,self.tiles_path,self.sc_name,width,self._mode,mag_threshold,t_pm)
        info = df2info(self.sc_name,center,df,max_control_points) 
        return Stars(info)

    def _search_draw(self,search_area):
        """
        Visualize the scope of the search area and the coverage of the corresponding tiles.

        Usage:
            >>> from starcatalogquery import StarCatalogReduced
            >>> # load the reduced star catalog GAIADR3
            >>> dir_from_reduced = '/Volumes/TOSHIBA/starcatalog/reduced/gaiadr3/res2/'
            >>> gaiadr3_reduced = StarCatalogReduced.load('reduced','gaiadr3',2,dir_from_reduced)
            >>> cone_area = {'cone':[20,30,10]}
            >>> stars = gaiadr3_reduced._search_draw(cone_area)

        Inputs:
            search_area -> [dict] The scope of the search area, such as {'cone':[20,30,10]} or {'box':[20,30,30,40]}, where
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
            An image that shows the scope of the search area and the coverage of the corresponding tiles.          
        """
        width = int(self.tile_size.split()[0]) 
        search_draw(width,search_area)   

class StarCatalogSimplified(object):
    """
    Class StarCatalogSimplified

    Attributes:
        - tiles_path: Path of the starcatalog tile files. 
        - sc_size: The size of the star catalog.
        - tiles_num: Total number of the tile files.
        - validity: The validity of the star catalog.
        - sc_name: Name of the starcatalog. Available starcatalogs include 'hygv3', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.
        - tile_size: Size of the tile in [deg]
        - _mode: Types of star catalogs, including 'raw', 'reduced', 'simplified', where
            'raw' represents the original star catalog, which contains all information about the star
            'reduced' represents the reduced star catalog, which contains the position, proper motion, apparent magnitude, epoch of the star
            'simplified' represents the minimalist star catalog, which only includes the position and apparent magnitude of stars
        - stars_num: Total number of stars in the catalog
        - mag: Apparent magnitude of stars in the catalog
        - description: Catalog Summary

    Methods:
        - load: Load the simplified star catalog files from the local database.
        - search_box: Perform a rectangle search of stars on the simplified star catalog and return an instance of class Stars.
        - search_cone: Perform a cone search of stars on the simplified star catalog and return an instance of class Stars.   
        - _search_draw: Visualize the scope of the search area and the coverage of the corresponding tiles.  
        - h5_incices: Generate a h5-formatted star catalog file, which records the center pointing of each sky area, the pixel coordinates of the stars, the triangle invariants and the asterism indices.
    """    
    def __init__(self,info):  

        self.info = info

        for key in info.keys():
            setattr(self, key, info[key])

    def __repr__(self):
    
        return 'Instance of class StarCatalogSimplified'        

    def load(sc_name,tile_size,mag_cutoff,epoch,dir_from=None):
        """
        Load the simplified star catalog files from the local database.

        Usage:
            >>> from starcatalogquery import StarCatalogSimplified
            >>> # load the simplified star catalog GAIADR3
            >>> dir_from_simplified = '/Volumes/TOSHIBA/starcatalog/simplified/gaiadr3/res2/'
            >>> gaiadr3_simplified = StarCatalogSimplified.load('simplified','gaiadr3',2,dir_from_simplified)

        Inputs:
            sc_name -> [str] Name of the starcatalog. Available starcatalogs include 'hygv3', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.
            tile_size -> [int] Size of the tile in [deg]
            mag_cutoff -> [float] The truncated magnitude of the simplified star catalog
            epoch -> [float] The epoch of the simplified star catalog
            dir_from -> [str,optional,default=None] The loading path of the star catalog files. If None, the path is automatically assigned to a suitable directory by default.

        Outputs:
            Instance of class StarCatalogSimplified
        """ 
        if dir_from is None: dir_from = 'starcatalogs/simplified/{:s}/res{:d}/mag{:.1f}/epoch{:.1f}/'.format(sc_name,tile_size,mag_cutoff,epoch)
        if not os.path.exists(dir_from): raise Exception('Path of the star catalog {:s} does not exist.'.format(sc_name))  

        # calculate total size and numbers of tile files    
        file_num,dir_size,validity = tiles_statistic(dir_from,tile_size) 
        stars_num,mag,description = starcatalog_info(sc_name)

        dict_values = dir_from,dir_size,file_num,validity,sc_name,'{:d} deg'.format(tile_size),'simplified',stars_num,mag,description,mag_cutoff,epoch
        dict_keys = 'tiles_path','sc_size','tiles_num','validity','sc_name','tile_size','_mode','stars_num','mag','description','mag_cutoff','epoch'
        info = dict(zip(dict_keys, dict_values))

        return StarCatalogSimplified(info)   

    def search_box(self,radec_box,max_control_points=None):
        """
        Perform a rectangle search of stars on the simplified star catalog and return an instance of class Stars.

        Usage:
            >>> from starcatalogquery import StarCatalogSimplified
            >>> # load the simplified star catalog GAIADR3
            >>> dir_from_simplified = '/Volumes/TOSHIBA/starcatalog/simplified/gaiadr3/res2/'
            >>> gaiadr3_simplified = StarCatalogSimplified.load('simplified','gaiadr3',2,dir_from_simplified)
            >>> stars = gaiadr3_simplified.search_box([20,30,30,40])

        Inputs:
            radec_box -> [int,array_like] Rectangular search area in format of [ra_min,dec_min,ra_max,dec_max], where
                ra_min -> [float] Left border of RA in [deg].
                dec_min -> [float] Lower border of DEC in [deg].
                ra_max -> [float] Right border of RA in [deg].
                dec_max -> [float] Upper border of DEC in [deg].
            max_control_points -> [int,optional,default=None] Number of brightest stars in the search area. If None, it is the number of all stars in the search area by deault.

        Outputs:
            Instance of class Stars
        """
        width = int(self.tile_size.split()[0])
        df = search_box(radec_box,self.tiles_path,self.sc_name,width,self._mode) 
        info = df2info(self.sc_name,center,df,max_control_points) 
        return Stars(info)

    def search_cone(self,center,radius,max_control_points=None):   
        """
        Perform a cone search of stars on the simplified star catalog and return an instance of class Stars.

        Usage:
            >>> from starcatalogquery import StarCatalogSimplified
            >>> # load the simplified star catalog GAIADR3
            >>> dir_from_simplified = '/Volumes/TOSHIBA/starcatalog/simplified/gaiadr3/res2/'
            >>> gaiadr3_simplified = StarCatalogSimplified.load('simplified','gaiadr3',2,dir_from_simplified)
            >>> stars = gaiadr3_simplified.search_cone([20,30],10)

        Inputs:
            center -> [int,array_like] Center of the cap in format of [ra_c,dec_c], where
                ra_c -> [float] RA, in [deg].
                dec_c -> [float] DEC, in [deg].
            radius -> [float] Angular radius of the cap, in [deg].
            max_control_points -> [int,optional,default=None] Number of brightest stars in the search area. If None, it is the number of all stars in the search area by deault.

        Outputs:
            Instance of class Stars
        """
        width = int(self.tile_size.split()[0])
        df = search_cone(center,radius,self.tiles_path,self.sc_name,width,self._mode)
        info = df2info(self.sc_name,center,df,max_control_points) 
        return Stars(info)

    def _search_draw(self,search_area):
        """
        Visualize the scope of the search area and the coverage of the corresponding tiles.

        Usage:
            >>> from starcatalogquery import StarCatalogSimplified
            >>> # load the simplified star catalog GAIADR3
            >>> dir_from_simplified = '/Volumes/TOSHIBA/starcatalog/simplified/gaiadr3/res2/'
            >>> gaiadr3_simplified = StarCatalogReduced.load('simplified','gaiadr3',2,dir_from_simplified)
            >>> cone_area = {'cone':[20,30,10]}
            >>> stars = gaiadr3_simplified._search_draw(cone_area)

        Inputs:
            search_area -> [dict] The scope of the search area, such as {'cone':[20,30,10]} or {'box':[20,30,30,40]}, where
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
            An image that shows the scope of the search area and the coverage of the corresponding tiles.          
        """
        width = int(self.tile_size.split()[0]) 
        search_draw(width,search_area)   

    def h5_incices(self,fov,pixel_width,max_control_points=60):
        """
        Generate a h5-formatted star catalog file, which records the center pointing of each sky area, the pixel coordinates of the stars, the triangle invariants and the asterism indices.

        Usage:
            >>> from starcatalogquery import StarCatalogSimplified
            >>> # load the simplified star catalog GAIADR3
            >>> dir_from_simplified = '/Volumes/TOSHIBA/starcatalog/simplified/gaiadr3/res2/'
            >>> gaiadr3_simplified = StarCatalogSimplified.load('simplified','gaiadr3',2,dir_from_simplified)
            >>> outh5 = gaiadr3_simplified.h5_incices(20,0.01)

        Inputs:
            fov -> [float] FOV of camera
            pixel_width -> [float] Pixel width in [deg]
            max_control_points -> [int,optional,default=60] Number of brightest stars in the search area.

        Outputs:
            h5-formatted star catalog file    
        """ 
        dir_h5 = 'starcatalogs/indices/{:s}/'.format(self.sc_name)
        Path(dir_h5).mkdir(parents=True, exist_ok=True)  

        outh5 = dir_h5 + 'fov{:d}_mag{:.1f}_mcp{:d}_{:.1f}.h5'.format(fov,self.mag_cutoff,max_control_points,self.epoch)

        if os.path.exists(outh5): return outh5 

        # write to h5 file
        fout = h5py.File(outh5,'w')
        # create group 
        stars_xy_grp = fout.create_group("stars_xy")
        stars_invariants_grp = fout.create_group("stars_invariants")
        stars_asterisms_grp = fout.create_group("stars_asterisms")

        if fov > 43: 
            nside = 1
            search_radius = 34
        elif fov > 22:
            nside = 2
            search_radius = 17
        else:
            nside = 4
            search_radius = 9

        npix = hp.nside2npix(nside)
        fp_radecs = []
        for seq in range(npix):

            desc = 'Generating starcatalog sky area index {:s}{:d}{:s} of {:d}'.format(Fore.BLUE,seq+1,Fore.RESET,npix)
            print(desc,end='\r')

            fp_radec = hp.pix2ang(nside,seq,lonlat=True)
            fp_radecs.append(fp_radec)

            stars = self.search_cone(fp_radec,search_radius,max_control_points)
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
        - sc_name: Name of the starcatalog.
        - center: Center pointing in format of [ra_c,dec_c] in [deg]
        - df: Dataframe of the stars
        - max_control_points: Number of stars in the dataframe
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
            >>> # load the simplified star catalog GAIADR3
            >>> dir_from_simplified = '/Volumes/TOSHIBA/starcatalog/simplified/gaiadr3/res2/'
            >>> gaiadr3_simplified = StarCatalogReduced.load('simplified','gaiadr3',2,dir_from_simplified)
            >>> gaiadr3_simplified_stars = gaiadr3_simplified.search_cone([20,30],10)
            >>> stars = gaiadr3_simplified_stars.pixel_xy(0.01)

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