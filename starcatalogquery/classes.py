# Load external packages
import os,warnings
from glob import glob
import numpy as np
import pandas as pd
from astropy.time import Time
from colorama import Fore
from natsort import natsorted
import sqlite3

# Load internal function library
from .catalog_index import build_catalog_indices,delete_table,generate_catalog_db,read_h5_hashes,h5_hashes
from .catalog_download import stsci_download,hyg_download
from .catalog_check import stsci_check
from .utils.starcatalog_statistic import tiles_statistic,starcatalog_info
from .utils.df2info import df2info
from .catalog_query import search_box_raw,search_cone_raw,search_box_reduced,search_cone_reduced,search_box_simplified,search_cone_simplified
from .tiles_draw import search_draw
from .wcs import xy_catalog
from .invariantfeatures import calculate_invariantfeatures
from .astrometry_corrections import parallax2dist

TILE_SIZE = 3.66  # Pixel size (in degrees) for HEALPix hierarchy level 'K4'

def format_columns(df, format_dict = {'ra': 8, 'dec': 8, 'pm_ra': 3, 'pm_dec': 3, 'dist': 8, 'mag': 3, 'epoch': 3}):
    """
    Formats the specified columns of a DataFrame according to the provided format dictionary.

    Inputs:
        df -> [pandas.DataFrame] The DataFrame to be formatted.
        format_dict -> [dict] A dictionary where the keys are column names and the values are the number of decimal places to round to.
    Outputs:
        df -> [pandas.DataFrame] The formatted DataFrame.
    """
    for col, decimals in format_dict.items():
        if col in df.columns:
            df[col] = df[col].round(decimals)
    return df

class StarCatalog(object):
    """
    StarCatalog class for handling star catalogs with varying levels of details, including 'raw', 'reduced', and 'simplified' versions.

    Levels of Catalogs:
        - 'raw': The original star catalog, encompassing comprehensive information about the stars,
        including their celestial coordinates, magnitudes, proper motions, parallax, and other astrophysical parameters, such as spectral type.
        - 'reduced': A more compact version of the star catalog, retaining key data such as position, proper motion, apparent magnitude, epoch, and parallax.
        - 'simplified': A basic catalog created from the reduced catalog by applying magnitude truncation and proper motion correction, suitable for quick lookups.

    Methods:
        - get: Downloads 'raw' star catalog files from remote servers or loads them from a local directory.
        - load: Load star catalogs at various levels of hierarchy from local storage.
    """

    def get(sc_name,dir_to=None):
        """
        Downloads 'raw' star catalog files from remote servers or loads them from a local directory.

        Usage:
            >>> from starcatalogquery import StarCatalog
            >>> sc_raw = StarCatalog.get('at-hyg24')
        Inputs:
            sc_name -> [str] Name of the star catalog. 
            Available star catalog include 'hyg37','at-hyg24','gaiadr3','gsc30','ucac5','usnob','2mass', etc.
            For more star catalogs, please refer to the Space Telescope Science Institute:
            https://outerspace.stsci.edu/display/GC/WebServices+for+Catalog+Access
            dir_to -> [str, optional, default=None] Directory for storing the downloaded catalog.
            If None, a built-in directory is used by default.
        Outputs:
            sc_raw -> [StarCatalogRaw object] An instance of the StarCatalogRaw class.
        """
        # List of supported star catalogs
        valid_catalogs = ['hyg37','at-hyg24','gaiadr3','gsc30','ucac5','usnob','2mass']
        if sc_name not in valid_catalogs:
            raise ValueError(f'Star catalog {sc_name} is not supported. Valid catalogs are: {valid_catalogs}')

        dir_starcatalogs = f'starcatalogs/raw/{sc_name}/'
        if dir_to is None:
            dir_to = dir_starcatalogs # Set default download directory if not specified
        else:
            dir_to = os.path.join(dir_to,dir_starcatalogs)   

        # Handle specific catalogs that require different download methods
        if sc_name in valid_catalogs[:2]:
            # For HYG and AT-HYG star catalogs
            dir_size, file_num, validity = hyg_download(sc_name,dir_to)
        else:
            # For GAIA DER3, GSC 30, UCAC5, USNOB, 2MASS star catalogs
            dir_url = dir_to.split('starcatalogs')[0] + 'starcatalogs/url/'
            url_file = os.path.join(dir_url, f'{sc_name}.txt')
            os.makedirs(dir_url, exist_ok=True)

            if os.path.exists(dir_to):
                # Check existing star catalog for validity
                dir_size,file_num,validity = stsci_check(sc_name,dir_to,url_file)
            else:
                # Download and check star catalog if not already present
                stsci_download(sc_name, dir_to, url_file)
                dir_size,file_num,validity = stsci_check(sc_name,dir_to,url_file)

        # Gather star catalog information
        stars_num,mag,description = starcatalog_info(sc_name)

        indices_dir = os.path.join(dir_to.split('starcatalogs')[0], 'starcatalogs', 'indices')
        tb_name = sc_name
        indices_path = os.path.join(indices_dir, f'{tb_name}.csv')
        db_path = os.path.join(indices_dir, 'catalogs.db')

        info = {
            'tiles_dir': dir_to,
            'sc_size': dir_size,
            'tiles_num': file_num,
            'validity': validity,
            'sc_name': sc_name,
            'tile_size': f'{TILE_SIZE} deg',
            '_mode': 'raw',
            'stars_num': stars_num,
            'mag': mag,
            'description': description,
            '_indices_dir': indices_dir,
            '_tb_name': tb_name,
            '_indices_path': indices_path,
            '_db_path': db_path
        }

        return StarCatalogRaw(info)  

    def load(dir_from):
        """
        Load star catalogs at various levels of hierarchy from local storage.

        Usage:
            >>> from starcatalogquery import StarCatalog
            >>> # load the raw star catalog AT-HYG v2.4
            >>> dir_from_raw = 'starcatalogs/raw/at-hyg24/'
            >>> at_hyg24_raw = StarCatalog.load(dir_from_raw)
            >>>
            >>> # load the reduced star catalog AT-HYG v2.4
            >>> dir_from_reduced = 'starcatalogs/reduced/at-hyg24/'
            >>> at_hyg24_reduced = StarCatalog.load(dir_from_reduced)
            >>>
            >>> # load the simplified star catalog AT-HYG v2.4
            >>> dir_from_simplified = 'starcatalogs/simplified/at-hyg24/mag12.0/epoch2019.5/'
            >>> at_hyg24_simplified = StarCatalog.load(dir_from_simplified)
        Inputs:
            dir_from -> [str] Directory for the star catalog files.
        Outputs:
            StarCatalogRaw | StarCatalogReduced | StarCatalogSimplified: 
            An instance of the appropriate star catalog class, depending on the catalog level identified in the directory.
        """
        # Check if the directory exists, raise an exception if it does not
        if not os.path.exists(dir_from): 
            raise Exception('Directory of the star catalog to load does not exist.')

        # Parse the directory to get components of the path
        parse_dir = dir_from.strip('/').split('/')

        if 'epoch' in parse_dir[-1]:
            _mode, sc_name, mag_threshold, epoch = parse_dir[-4:]
            tb_name = f'{sc_name}_{mag_threshold}_{epoch}'
        else:
            _mode, sc_name = parse_dir[-2:]
            tb_name = sc_name

        # load the star catalog
        if _mode == 'raw':
            starcatalog = StarCatalogRaw._load(sc_name,dir_from)
        elif _mode == 'reduced':
            starcatalog = StarCatalogReduced._load(sc_name,dir_from)
        elif _mode == 'simplified':
            mag_threshold = float(mag_threshold[3:])
            epoch = float(epoch[5:])
            starcatalog = StarCatalogSimplified._load(sc_name,mag_threshold,epoch,dir_from)
        else:
            raise Exception(f"Invalid star catalog level: {_mode}.")

        indices_dir = os.path.join(dir_from.split('starcatalogs')[0], 'starcatalogs', 'indices')
        indices_path = os.path.join(indices_dir, f'{tb_name}.csv')
        db_path = os.path.join(indices_dir, 'catalogs.db')

        starcatalog._tb_name = tb_name
        starcatalog._indices_path = indices_path
        starcatalog._db_path = db_path
        starcatalog._indices_dir = indices_dir

        return starcatalog

class StarCatalogRaw(object):
    """
    This version of the catalog is the original star catalog, encompassing comprehensive information about the stars,
    including their celestial coordinates, magnitudes, proper motions, parallax, and other astrophysical parameters.

    Attributes:
        - tiles_dir (str): Directory where the star catalog files are stored.
        - tiles_num (int): Total number of the tile files in the catalog.
        - tile_size (str): Geometric size of each tile in the catalog, in degrees.
        - sc_size (str): Overall size of the star catalog (e.g., in MB or GB).
        - validity (bool): Indicates whether the star catalog is complete and valid for use.
        - sc_name (str): Name of the star catalog, such as 'at-hyg24'.
        - _mode (str): Level of the star catalog, set to 'raw' for this class.
        - stars_num (str): Total number of the stars included in the catalog.
        - mag (str): Range of stars magnitudes in the catalog.
        - description (str): A brief description or summary of the catalog.

    Methods:
        - _load: Loads the raw star catalog from a specified directory.
        - reduce: Reduces the raw star catalog to contain only essential information.
        - search_box: Searches for stars within a specified rectangular area.
        - search_cone: Searches for stars within a specified conical area.
    """   
    def __init__(self,info): 
        """
        Initializes a new instance of the StarCatalogRaw class.

        Inputs:
            info -> [dict] A dictionary containing key-value pairs of catalog attributes.
        """
        self.__dict__.update(info)

    def __repr__(self):
        """
        Returns a string representation of the StarCatalogRaw instance.

        Returns:
            A string describing the StarCatalogRaw instance.
        """
    
        return "<StarCatalogRaw object: CATALOG_NAME = '{:s}' CATALOG_SIZE = '{:s}' TILES_NUM = {:d} TILE_SIZE = '{:s}' STARS_NUM = '{:s}' MAG = '{:s}'>".format(
            self.sc_name,self.sc_size,self.tiles_num,self.tile_size,self.stars_num,self.mag)

    @classmethod
    def _load(cls, sc_name, dir_from):
        """
        Load the raw star catalog files from a specified directory.

        Usage:
            >>> from starcatalogquery import StarCatalogRaw
            >>> # load the raw star catalog AT-HYG v2.4
            >>> dir_from_raw = 'starcatalogs/raw/at-hyg24/'
            >>> at_hyg24_raw = StarCatalogRaw.load('at-hyg24',dir_from_raw)
        Inputs:
            sc_name -> [str] Name of the star catalog (e.g., 'at-hyg24').
            dir_from -> [str] Directory for the star catalog files.
        Outputs:
            at_hyg24_raw -> [StarCatalogRaw object] Instance of the class StarCatalogRaw.
        """

        # Calculate total file size and total number of the tile files
        file_num, dir_size, validity = tiles_statistic(dir_from)
        stars_num, mag, description = starcatalog_info(sc_name)

        info = {
            'tiles_dir': dir_from,
            'sc_size': dir_size,
            'tiles_num': file_num,
            'validity': validity,
            'sc_name': sc_name,
            'tile_size': f'{TILE_SIZE} deg',
            '_mode': 'raw',
            'stars_num': stars_num,
            'mag': mag,
            'description': description
        }

        return cls(info)    

    def reduce(self):
        """
        Reduce the raw star catalog, retaining only essential star information.

        Usage:
            >>> at_hyg24_reduced = at_hyg24_raw.reduce()
        Outputs:
            at_hyg24_reduced -> [StarCatalogReduced object] The reduced star catalog.
        """
        # Copying existing information
        info = self.__dict__.copy()

        # Setting up the directory for the reduced catalog
        tiles_dir = self.tiles_dir
        dir_reduced = tiles_dir.replace('raw','reduced')
        os.makedirs(dir_reduced, exist_ok=True)

        file_list = glob(tiles_dir+'*') 
        sc_name = self.sc_name 
        print(f'Reducing the star catalog {sc_name}, which may take a considerable amount of time...')

        if sc_name in ['hyg37','at-hyg24']:
            # Processing each file in the raw star catalog directory
            for j, tile_file in enumerate(file_list, start=1):
                desc = f'Reducing {Fore.BLUE}{j}{Fore.RESET} of {self.tiles_num}'
                print(desc,end='\r')
                df = pd.read_csv(tile_file,skiprows=1,dtype=str,na_values=[' ', ''])
                # Applying specific transformations based on the catalog name
                # units: ra->hourangle, dec->deg, pmra->mas/a, pmdec->mas/a, epoch->2000.0
                columns_dict = {'pmra':'pm_ra', 'pmdec':'pm_dec'}
                df.rename(columns=columns_dict, inplace=True)
                df_reduced = df.loc[:,['ra','dec','pm_ra','pm_dec','dist','mag']] # distance is in parsecs
                df_reduced['epoch'] = '2000.0'
                df_reduced = df_reduced.astype(float)
                df_reduced['ra'] = df_reduced['ra']*15 # Convert hourangle to deg
                df_reduced['dist'] = df_reduced['dist'] / 1e3  # Convert pc to kpc
                df_reduced = format_columns(df_reduced)
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)
        elif sc_name == 'gaiadr3':
            for j, tile_file in enumerate(file_list, start=1):
                desc = f'Reducing {Fore.BLUE}{j}{Fore.RESET} of {self.tiles_num}'
                print(desc,end='\r')
                df = pd.read_csv(tile_file,skiprows=1,dtype=str,na_values=[' ', ''])
                # units: ra->deg, dec->deg, pmra->mas/a, pmdec->mas/a, epoch->2016
                df_reduced = df.loc[:,['ra','dec','pmra','pmdec','parallax','mag','epoch']]
                df_reduced = df_reduced.astype(float)
                # Calculate distance from parallax
                df_reduced['dist'] = parallax2dist(df_reduced['parallax'])
                parallax_index = df_reduced.columns.get_loc('parallax')
                # Remove the 'parallax' column from the DataFrame
                df_reduced.drop(columns=['parallax'], inplace=True)
                # Insert the 'dist' column in the position of the original 'parallax' column
                df_reduced.insert(parallax_index, 'dist', df_reduced.pop('dist'))
                # Rename columns
                columns_dict = {'pmra': 'pm_ra', 'pmdec': 'pm_dec'}
                df_reduced.rename(columns=columns_dict, inplace=True)
                df_reduced = format_columns(df_reduced)
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)
        elif sc_name == 'gsc30':
            for j, tile_file in enumerate(file_list, start=1):
                desc = f'Reducing {Fore.BLUE}{j}{Fore.RESET} of {self.tiles_num}'
                print(desc,end='\r')
                df = pd.read_csv(tile_file,skiprows=1,dtype=str,na_values=[' ', ''])
                # units: ra->deg, dec->deg, rapm->mas/a, decpm->mas/a, epoch->2012
                df_reduced = df.loc[:,['ra','dec','rapm','decpm','parallax','mag','epoch']]
                df_reduced = df_reduced.astype(float)
                # Calculate distance from parallax
                df_reduced['dist'] = parallax2dist(df_reduced['parallax'])
                # Get the position of the 'parallax' column
                parallax_index = df_reduced.columns.get_loc('parallax')
                # Remove the 'parallax' column from the DataFrame
                df_reduced.drop(columns=['parallax'], inplace=True)
                # Insert the 'dist' column in the position of the original 'parallax' column
                df_reduced.insert(parallax_index, 'dist', df_reduced.pop('dist'))
                # Rename columns
                columns_dict = {'rapm':'pm_ra', 'decpm':'pm_dec'}
                df_reduced.rename(columns=columns_dict, inplace=True)
                df_reduced = format_columns(df_reduced)
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)
        elif sc_name == 'ucac5':
            for j, tile_file in enumerate(file_list, start=1):        
                desc = f'Reducing {Fore.BLUE}{j}{Fore.RESET} of {self.tiles_num}'
                print(desc,end='\r')
                df = pd.read_csv(tile_file,skiprows=1,dtype=str,na_values=[' ', ''])
                # units: ra->deg, dec->deg, pmur->mas/a, pmud->mas/a, epu->1998.754
                df_reduced = df.loc[:,['ra','dec','pmur','pmud','mag','epu']]
                df_reduced = df_reduced.astype(float)
                columns_dict = {'pmur':'pm_ra', 'pmud':'pm_dec','epu':'epoch'}
                df_reduced.rename(columns=columns_dict, inplace=True)
                df_reduced = format_columns(df_reduced)
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)
        elif sc_name == 'usnob':
            for j, tile_file in enumerate(file_list, start=1):        
                desc = f'Reducing {Fore.BLUE}{j}{Fore.RESET} of {self.tiles_num}'
                print(desc,end='\r')
                df = pd.read_csv(tile_file,skiprows=1,dtype=str,na_values=[' ', ''])
                # units: ra->deg, dec->deg, pmRA->mas/a, pmDEC->mas/a, Epoch->1950,mag->unknown
                df_reduced = df.loc[:,['ra','dec','pmRA','pmDEC','mag','Epoch']]
                df_reduced = df_reduced.astype(float)
                columns_dict = {'pmRA':'pm_ra', 'pmDEC':'pm_dec','Epoch':'epoch'}
                df_reduced.rename(columns=columns_dict, inplace=True)
                df_reduced = format_columns(df_reduced)
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)  
        elif sc_name == '2mass':
            for j, tile_file in enumerate(file_list, start=1):   
                desc = f'Reducing {Fore.BLUE}{j}{Fore.RESET} of {self.tiles_num}'
                print(desc,end='\r')
                df = pd.read_csv(tile_file,skiprows=1,dtype=str,na_values=[' ', ''])
                # units: ra->deg, dec->deg, jdate->2451063.6417
                df_reduced = df.loc[:,['ra','dec','mag','jdate']]
                df_reduced = df_reduced.astype(float)
                # Convert Julian date to epoch year
                df_reduced['epoch'] = Time(df_reduced['jdate'], format='jd').jyear
                df_reduced.pop('jdate')
                df_reduced = format_columns(df_reduced)
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)

        print('\nFinished')
                    
        file_num,dir_size,validity = tiles_statistic(dir_reduced)  
        info['_mode'] = 'reduced'
        info['tiles_dir'] = dir_reduced
        info['sc_size'] = dir_size

        return StarCatalogReduced(info)  
    
    def search_box(self,radec_box,mag_threshold,t_pm,fov_min=None,max_num=None,max_num_per_tile=None):
        """
        Perform a rectangle search of stars on raw star catalogs.

        Usage:
            >>> stars = at_hyg24_raw.search_box([20,30,30,40],9,2022.0)
        Inputs:
            radec_box -> [list] Rectangular search area in form of [ra_min, dec_min, ra_max, dec_max] in degrees.
            mag_threshold -> [float] Apparent magnitude limit.
            t_pm -> [float] Epoch to which the stars are unified.
            fov_min -> [float,optional,default=None] Field of view parameters in degrees.
            It determines the hierarchical division of the sky region in HEALPix,
            ensuring that each tile's size is greater than one-quarter of the FOV and less than or equal to half of the FOV.
            If None, it takes the short side of the rectangle as the value.
            max_num -> [int,optional,default=None] Maximum number of stars to include in the search, sorted by brightness.
            If None, includes all stars meeting the magnitude criteria.
            max_num_per_tile [int,optional,default=None] Maximum number of stars per tile to include in the search.
        Outputs:
            stars -> [Stars object] An instance of Stars with the search results.
        """
        # Calculating the center of the search box
        ra_min,dec_min,ra_max,dec_max = radec_box
        center = [(ra_min+ra_max)/2,(dec_min+dec_max)/2]

        search_area = {'box':radec_box}

        dir_sc = self.tiles_dir
        sc_name,_mode = self.sc_name,self._mode
        tb_name = self._tb_name
        catalog_indices_db = self._db_path

        # Performs a rectangular search on the raw star catalog considering magnitude and proper motion.
        df,level,nside,ids,pixel_size,fov_min = search_box_raw(radec_box, dir_sc, sc_name, _mode, tb_name, catalog_indices_db,mag_threshold,t_pm,fov_min,max_num_per_tile)
        info = df2info(sc_name,center,df,max_num,level,nside,ids,pixel_size,search_area,fov_min)
        return Stars(info)

    def search_cone(self,center,radius,mag_threshold,t_pm,fov_min=None,max_num=None,max_num_per_tile=None):
        """
        Perform a conical search of stars on raw star catalogs.

        Usage:
            >>> stars = at_hyg24_raw.search_cone([30, 40], 5, 9, 2022.0)
        Inputs:
            center -> [list] Center of the cone search area in the form of [ra, dec] in degrees.
            radius -> [float] Radius of the cone search area in degrees.
            mag_threshold -> [float] Apparent magnitude limit.
            t_pm -> [float] Epoch to which the stars are unified.
            fov_min -> [float,optional,default=None] Field of view parameters in degrees.
            It determines the hierarchical division of the sky region in HEALPix,
            ensuring that each tile's size is greater than one-quarter of the FOV and less than or equal to half of the FOV.
            If None, it takes the search diameter as the value.
            max_num -> [int,optional,default=None] Maximum number of stars to include in the search, sorted by brightness.
            If None, includes all stars meeting the magnitude criteria.
            max_num_per_tile [int,optional,default=None] Maximum number of stars per tile to include in the search.
        Outputs:
            Stars: An instance of Stars with the search results.
        """

        search_area = {'cone':(center,radius)}

        dir_sc = self.tiles_dir
        sc_name,_mode = self.sc_name,self._mode
        tb_name = self._tb_name

        catalog_indices_db = self._db_path

        # Performs a conical search of stars on the raw star catalog considering magnitude and proper motion.
        df,level,nside,ids,pixel_size,fov_min = search_cone_raw(center,radius,dir_sc,sc_name,_mode,tb_name,catalog_indices_db,mag_threshold,t_pm,fov_min,max_num_per_tile)
        info = df2info(sc_name,center,df,max_num,level,nside,ids,pixel_size,search_area,fov_min)
        return Stars(info)      

class StarCatalogReduced(object):
    """
    This version of the catalog represents the reduced star catalog, which retains key data such as
    the position, proper motion, apparent magnitude, parallax, and epoch of the stars, making it a more concise version of the raw catalog.

    Attributes:
        - tiles_dir (str): Directory where the star catalog files are stored.
        - tiles_num (int): Total number of the tile files in the catalog.
        - tile_size (str): Geometric size of each tile in the catalog, in degrees.
        - sc_size (str): Overall size of the star catalog (e.g., in MB or GB).
        - validity (bool): Indicates whether the star catalog is complete and valid for use.
        - sc_name (str): Name of the star catalog, such as 'at-hyg24'.
        - _mode (str): Level of the star catalog, set to 'reduced' for this class.
        - stars_num (str): Total number of the stars included in the catalog.
        - mag (str): Range of star magnitudes in the catalog.
        - description (str): A brief description or summary of the catalog.

    Methods:
        - _load: Loads the reduced star catalog from a specified directory.
        - build_indices: Aggregates indices from all tile files and generates a comprehensive catalog index file.
        - simplify: Simplifies the reduced star catalog to a more basic version.
        - search_box: Searches for stars within a specified rectangular area.
        - search_cone: Searches for stars within a specified conical area.
    """ 
    def __init__(self,info):
        """
        Initializes a new instance of the StarCatalogReduced class.

        Inputs:
            info -> [dict] A dictionary containing key-value pairs of catalog attributes.
        """
        self.__dict__.update(info)

    def __repr__(self):
        """
        Returns a string representation of the StarCatalogReduced instance.

        Returns:
            A string describing the StarCatalogReduced instance.
        """
    
        return "<StarCatalogReduced object: CATALOG_NAME = '{:s}' CATALOG_SIZE = '{:s}' TILES_NUM = {:d} TILE_SIZE = '{:s}' STARS_NUM = '{:s}' MAG = '{:s}'>".format(
            self.sc_name,self.sc_size,self.tiles_num,self.tile_size,self.stars_num,self.mag)

    @classmethod
    def _load(cls,sc_name,dir_from):
        """
        Load the reduced star catalog files from a specified directory.

        Usage:
            >>> from starcatalogquery import StarCatalogReduced
            >>> # load the reduced star catalog HYG v3.7
            >>> dir_from_reduced = 'starcatalogs/reduced/at-hyg24/'
            >>> at_hyg24_reduced = StarCatalogReduced.load('at-hyg24',dir_from_reduced)
        Inputs:
            sc_name -> [str] Name of the star catalog (e.g., 'at-hyg24').
            dir_from -> [str] Directory for the star catalog files.
        Outputs:
            at_hyg24_reduced -> [StarCatalogReduced object] Instance of the class StarCatalogReduced.
        """ 

        # Calculate total file size and total number of the tile files
        file_num,dir_size,validity = tiles_statistic(dir_from) 
        stars_num,mag,description = starcatalog_info(sc_name)

        info = {
            'tiles_dir': dir_from,
            'sc_size': dir_size,
            'tiles_num': file_num,
            'validity': validity,
            'sc_name': sc_name,
            'tile_size': f'{TILE_SIZE} deg',
            '_mode': 'reduced',
            'stars_num': stars_num,
            'mag': mag,
            'description': description
        }

        return cls(info)  

    def build_indices(self):
        """
        Aggregate indices from all tile files for each star across different HEALPix hierarchy levels,
        generating a comprehensive csv-formatted catalog index file.

        Usage:
            >>> at_hyg24_reduced.build_indices()
        """
        indices_path = build_catalog_indices(self.tiles_dir,self.sc_name,self._tb_name)

    def simplify(self,mag_threshold,t_pm):
        """
        Simplify the reduced star catalog, applying magnitude truncation and proper motion correction，
        making it a basic version suitable for quick lookups.
        
        Usage:
            >>> at_hyg24_simplified = at_hyg24_reduced.simplify(9.0,2022.0)
        Inputs:
            mag_threshold -> [float] Apparent magnitude limit.
            t_pm -> [float] Epoch to which the stars are unified.
        Outputs:
            Instance of class StarCatalogSimplified      
        """
        # Copying existing information
        info = self.__dict__.copy()

        # Setting up the directory for the reduced catalog
        tiles_dir = self.tiles_dir
        dir_simplified = tiles_dir.replace('reduced','simplified')+'mag{:.1f}/epoch{:.1f}/'.format(mag_threshold,t_pm)
        os.makedirs(dir_simplified, exist_ok=True)

        # Construct file pattern
        sc_name = self.sc_name
        file_pattern = os.path.join(tiles_dir, f'{sc_name}-*.csv')
        file_list = natsorted(glob(file_pattern))

        print('Simplifying the star catalog {:s}, which may take a considerable amount of time'.format(sc_name))  

        # Processing each file in the raw star catalog directory
        for j, tile_file in enumerate(file_list, start=1):
            desc = 'Simplifying {:s}{:d}{:s} of {:d}'.format(Fore.BLUE,j,Fore.RESET,self.tiles_num)
            print(desc,end='\r')

            df_simplified = pd.read_csv(tile_file).dropna()
            mag = df_simplified['mag']
            mag_flag = np.abs(mag) < mag_threshold
            df_simplified = df_simplified[mag_flag].sort_values(by=['mag'])

            # Correct proper motion    
            dt = float(t_pm) - df_simplified['epoch']  

            if {'pm_ra', 'pm_dec'}.issubset(df_simplified.columns):
                df_simplified['ra'] +=  df_simplified['pm_ra']/3.6e6 * dt
                df_simplified['dec'] += df_simplified['pm_dec']/3.6e6 * dt
                df_simplified['pm_ra'] = df_simplified['pm_ra']
                df_simplified['pm_dec'] = df_simplified['pm_dec']
                df_simplified['epoch'] = t_pm
            else:
                warnings.warn('Proper motion data for stars in catalog {:s} are not found.'.format(sc_name))

            if 'dist' in df_simplified.columns:
                dist_flag = df_simplified['dist'] > 0
                df_simplified = df_simplified[dist_flag]

            df_simplified = format_columns(df_simplified)
            df_simplified.to_csv(dir_simplified + tile_file.split('/')[-1],index=False)

        print('\nFinished')

        tb_name = '{:s}_mag{:.1f}_epoch{:.1f}'.format(sc_name, mag_threshold, t_pm)
        indices_path = os.path.join(self._indices_dir, f'{tb_name}.csv')

        file_num,dir_size,validity = tiles_statistic(dir_simplified)  
        info['_mode'] = 'simplified'
        info['tiles_dir'] = dir_simplified
        info['sc_size'] = dir_size
        info['mag_threshold'] = mag_threshold
        info['epoch'] = t_pm
        info['_tb_name'] = tb_name
        info['_indices_path'] = indices_path

        return StarCatalogSimplified(info)    

    def search_box(self,radec_box,mag_threshold,t_pm,fov_min=None,max_num=None,max_num_per_tile=None):
        """
        Perform a rectangle search of stars on reduced star catalogs.

        Usage:
            >>> stars = at_hyg24_reduced.search_box([20,30,30,40],9,2022.0)
        Inputs:
            radec_box -> [list] Rectangular search area in form of [ra_min, dec_min, ra_max, dec_max] in degrees.
            mag_threshold -> [float] Apparent magnitude limit.
            t_pm -> [float] Epoch to which the stars are unified.
            fov_min -> [float,optional,default=None] Field of view parameters in degrees.
            It determines the hierarchical division of the sky region in HEALPix,
            ensuring that each tile's size is greater than one-quarter of the FOV and less than or equal to half of the FOV.
            If None, it takes the search diameter as the value.
            max_num -> [int,optional,default=None] Maximum number of stars to include in the search, sorted by brightness.
            If None, includes all stars meeting the magnitude criteria.
            max_num_per_tile [int,optional,default=None] Maximum number of stars per tile to include in the search.
        Outputs:
            stars -> [Stars object] An instance of Stars with the search results.
        """
        # Calculating the center of the search box
        ra_min,dec_min,ra_max,dec_max = radec_box
        center = [(ra_min+ra_max)/2,(dec_min+dec_max)/2]

        search_area = {'box':radec_box}

        dir_sc = self.tiles_dir
        sc_name,_mode = self.sc_name,self._mode
        tb_name = self._tb_name
        catalog_indices_db = self._db_path

        # Performs a rectangular search on the reduced star catalog considering magnitude and proper motion.
        df,level,nside,ids,pixel_size,fov_min = search_box_reduced(radec_box, dir_sc, sc_name, _mode, tb_name, catalog_indices_db,mag_threshold,t_pm,fov_min,max_num_per_tile)
        info = df2info(sc_name,center,df,max_num,level,nside,ids,pixel_size,search_area,fov_min)
        return Stars(info)

    def search_cone(self,center,radius,mag_threshold,t_pm,fov_min=None,max_num=None,max_num_per_tile=None):
        """
        Perform a conical search of stars on reduced star catalogs.

        Usage:
            >>> stars = at_hyg24_reduced.search_cone([20,30],10,9,2022.0)
        Inputs:
            center -> [list] Center of the cone search area in form of [ra_c, dec_c] in degrees.
            radius -> [float] Angular radius of the cone.
            mag_threshold -> [float] Apparent magnitude limit.
            t_pm -> [float] Epoch to which the stars are unified.
            fov_min -> [float,optional,default=None] Field of view parameters in degrees.
            It determines the hierarchical division of the sky region in HEALPix,
            ensuring that each tile's size is greater than one-quarter of the FOV and less than or equal to half of the FOV.
            If None, it takes the search diameter as the value.
            max_num -> [int,optional,default=None] Maximum number of stars to include in the search, sorted by brightness.
            If None, includes all stars meeting the magnitude criteria.
            max_num_per_tile [int,optional,default=None] Maximum number of stars per tile to include in the search.
        Outputs:
            stars -> [Stars object] An instance of Stars with the search results.
        """
        search_area = {'cone':(center,radius)}

        dir_sc = self.tiles_dir
        sc_name,_mode = self.sc_name,self._mode
        tb_name = self._tb_name
        catalog_indices_db = self._db_path

        # Performs a conical search of stars on the reduced star catalog considering magnitude and proper motion.
        df,level,nside,ids,pixel_size,fov_min = search_cone_reduced(center,radius,dir_sc,sc_name,_mode,tb_name,catalog_indices_db,mag_threshold,t_pm,fov_min,max_num_per_tile)
        info = df2info(sc_name,center,df,max_num,level,nside,ids,pixel_size,search_area,fov_min)
        return Stars(info)  

class StarCatalogSimplified(object):
    """
    This version of star catalog is streamlined to apply magnitude truncation and proper motion correction at a specific epoch.

    Attributes:
        - tiles_dir (str): Directory where the star catalog files are stored.
        - tiles_num (int): Total number of the tile files in the catalog.
        - tile_size (str): Size of each tile in the catalog, in degrees.
        - sc_size (str): Overall size of the star catalog (e.g., in MB or GB).
        - validity (bool): Indicates whether the star catalog is complete and valid for use.
        - sc_name (str): Name of the star catalog, such as 'hygv3.7'.
        - _mode (str): Level of the catalog, set to 'simplified' for this class.
        - stars_num (str): Total number of the stars included in the catalog.
        - mag (str): Range of star magnitudes in the catalog.
        - description (str): A brief description or summary of the catalog.
        - _mode_invariants -> [str] Mode of geometric invariants after executing h5_hashes, e.g., 'triangles', 'quads'.
        - hashed_h5 -> [str] A h5-formatted hashed file for star catalog geometric invariants
        - hashed_data -> [dict] A dictionary containing the geometric invariants data from the h5 file.

    Methods:
        - _load: Loads the simplified star catalog from a specified directory.
        - search_box: Searches for stars within a specified rectangular area.
        - search_cone: Searches for stars within a specified conical area.
        - build_indices: Aggregates indices from all tile files and generates a comprehensive catalog index file.
        - h5_hashes: Generates a h5-formatted star catalog geometric invariants hashed file.
        - read_h5_hashes: Reads the h5-formatted hashed file, essential for blind star pattern matching.
    """   
    def __init__(self,info):
        """
        Initializes a new instance of the StarCatalogSimplified class.

        Inputs:
            info -> [dict] A dictionary containing key-value pairs of catalog attributes.
        """

        self.__dict__.update(info)

    def __repr__(self):
        """
        Returns a string representation of the StarCatalogSimplified instance.

        Returns:
            A string describing the StarCatalogSimplified instance.
        """

        return "<StarCatalogSimplified object: CATALOG_NAME = '{:s}' CATALOG_SIZE = '{:s}' TILES_NUM = {:d} TILE_SIZE = '{:s}' STARS_NUM = '{:s}' MAG = '{:s}' MAG_CUTOFF = {:.1f} EPOCH = {:.1f}>".format(
            self.sc_name,self.sc_size,self.tiles_num,self.tile_size,self.stars_num,self.mag,self.mag_threshold,self.epoch)  

    @classmethod
    def _load(cls,sc_name,mag_threshold,epoch,dir_from):
        """
        Load the simplified star catalog files from a specified directory.

        Usage:
            >>> from starcatalogquery import StarCatalogSimplified
            >>> # load the simplified star catalog AT-HYG v2.4
            >>> dir_from_simplified = 'starcatalogs/simplified/at_hyg24/mag12.0/epoch2019.5/'
            >>> at_hyg24_simplified = StarCatalogSimplified.load('at-hyg24',5,9,2022,dir_from_simplified)
        Inputs:
            sc_name -> [str] Name of the star catalog (e.g., 'at-hyg24').
            dir_from -> [str] Directory for the star catalog files.
        Outputs:
            at_hyg24_simplified -> [StarCatalogSimplified object] Instance of the class StarCatalogSimplified with loaded data.
        """  

        # Calculate total size and number of tile files   
        file_num,dir_size,validity = tiles_statistic(dir_from) 
        stars_num,mag,description = starcatalog_info(sc_name)

        info = {
            'tiles_dir': dir_from,
            'sc_size': dir_size,
            'tiles_num': file_num,
            'validity': validity,
            'sc_name': sc_name,
            'tile_size': f'{TILE_SIZE} deg',
            '_mode': 'simplified',
            'stars_num': stars_num,
            'mag': mag,
            'description': description,
            'mag_threshold': mag_threshold,
            'epoch': epoch
        }

        return cls(info) 

    def build_indices(self):
        """
        Aggregate indices from all tile files for each star across different HEALPix hierarchy levels,
        generating a comprehensive csv-formatted catalog index file.

        Usage:
            >>> at_hyg24_simplified.build_indices()
        """
        indices_path = build_catalog_indices(self.tiles_dir,self.sc_name,self._tb_name)

    def search_box(self,radec_box,fov_min=None,max_num=None,max_num_per_tile=None,astrometry_corrections={}):
        """
        Perform a rectangle search of stars on simplified star catalogs.

        Usage:
            >>> stars = at_hyg24_simplified.search_box([20,30,30,40])
        Inputs:
            radec_box -> [list] Rectangular search area in form of [ra_min, dec_min, ra_max, dec_max] in degrees.
            fov_min -> [float,optional,default=None] Field of view parameters in degrees.
            It determines the hierarchical division of the sky region in HEALPix,
            ensuring that each tile's size is greater than one-quarter of the FOV and less than or equal to half of the FOV.
            If None, it takes the short side of the rectangle as the value.
            max_num -> [int,optional,default=None] Maximum number of stars to include in the search, sorted by brightness.
            If None, includes all stars meeting the magnitude criteria.
            max_num_per_tile -> [int,optional,default=None] Maximum number of stars per tile to include in the search.
            astrometry_corrections -> [dict, optional, default={}] Dictionary specifying the types of astrometry corrections to apply.
                - 't' -> [str] Observation time in UTC, such as '2019-02-26T20:11:14.347'.
                - 'proper-motion' -> [None] If present, apply proper motion correction.
                - 'aberration' -> [tuple] Aberration correction parameters. Observer's velocity relative to Earth's center (vx, vy, vz) in km/s.
                - 'parallax' -> [None] If present, apply parallax correction.
                - 'deflection' -> [None] If present, apply light deflection correction.
        Outputs:
            stars -> [Stars object] An instance of Stars with the search results.
        """
        # Calculating the center of the search box
        ra_min,dec_min,ra_max,dec_max = radec_box
        center = [(ra_min+ra_max)/2,(dec_min+dec_max)/2]

        search_area = {'box':radec_box}

        dir_sc = self.tiles_dir
        sc_name,_mode = self.sc_name,self._mode
        tb_name = self._tb_name
        catalog_indices_db = self._db_path

        # Performs a rectangular search on the simplified star catalog
        df,level,nside,ids,pixel_size,fov_min = search_box_simplified(radec_box, dir_sc, sc_name, _mode, tb_name, catalog_indices_db,fov_min, max_num_per_tile,astrometry_corrections)
        info = df2info(sc_name,center,df,max_num,level,nside,ids,pixel_size,search_area,fov_min)
        return Stars(info)

    def search_cone(self,center,radius,fov_min=None,max_num=None,max_num_per_tile=None,astrometry_corrections={}):
        """
        Perform a conical search of stars on simplified star catalogs.

        Usage:
            >>> stars = at_hyg24_simplified.search_cone([20,30],10)
        Inputs:
            center -> [list] Center of the cone in form of [ra_c, dec_c] in degrees.
            radius -> [float] Angular radius of the cone.
            fov_min -> [float,optional,default=None] Field of view parameters in degrees.
            max_num -> [int,optional,default=None] Maximum number of stars to include in the search, sorted by brightness.
            If None, includes all stars meeting the magnitude criteria.
            max_num_per_tile [int,optional,default=None] Maximum number of stars per tile to include in the search.
            astrometry_corrections -> [dict, optional, default={}] Dictionary specifying the types of astrometry corrections to apply.
                - 't' -> [str] Observation time in UTC, such as '2019-02-26T20:11:14.347'.
                - 'proper-motion' -> [None] If present, apply proper motion correction.
                - 'aberration' -> [tuple] Aberration correction parameters. Observer's velocity relative to Earth's center (vx, vy, vz) in km/s.
                - 'parallax' -> [None] If present, apply parallax correction.
                - 'deflection' -> [None] If present, apply light deflection correction.
        Outputs:
            stars -> [Stars object] An instance of class Stars with the search results.
        """
        search_area = {'cone':(center,radius)}

        dir_sc = self.tiles_dir
        sc_name,_mode = self.sc_name,self._mode
        tb_name = self._tb_name
        catalog_indices_db = self._db_path

        # Performs a cone search on the simplified star catalog
        df,level,nside,ids,pixel_size,fov_min = search_cone_simplified(center,radius,dir_sc,sc_name,_mode,tb_name,catalog_indices_db,fov_min,max_num_per_tile,astrometry_corrections)
        info = df2info(sc_name,center,df,max_num,level,nside,ids,pixel_size,search_area,fov_min)
        return Stars(info)  

    def h5_hashes(self,k_min,k_max,mode_invariants='triangles'):
        """
        Generate a h5-formatted star catalog geometric invariants hashed file.

        Usage:
            >>> at_hyg24_simplified.h5_hashes(1, 6)
        Inputs:
            k_min -> [int] Minimum HEALPix hierarchy level.
            k_max -> [int] Maximum HEALPix hierarchy level.
            mode_invariants -> [str, optional, default='triangles'] Mode of geometric invariants to use, e.g., 'triangles', 'quads'.
        Outputs:
            hashed_h5 -> Path to the h5-formatted hashed file.
        """
        hashed_h5 = h5_hashes(self._db_path, self._tb_name, self.tiles_dir, self.sc_name, k_min,k_max, mode_invariants)
        self.hashed_h5 = hashed_h5
        self._mode_invariants = mode_invariants

        return self

    def read_h5_hashes(self,infile=None):
        """
        Reads an h5-formatted star catalog hashed file that contains geometric invariants of star configurations.
        These invariants include triangle edge length ratios or quadrilateral invariants, which are used for
        tasks like blind star pattern recognition and matching.

        Inputs:
            infile -> [str,optional,default=None] Path to the h5 hashed file.
        Outputs:
            data -> [dict] A dictionary containing the geometric invariants data from the h5 file.
        """

        if infile is None:
            if not hasattr(self,'hashed_data'):
                if not hasattr(self,'hashed_h5'):
                    raise Exception('No h5 hashed file provided.')
                else:
                    infile = self.hashed_h5
                    self.hashed_data = read_h5_hashes(infile)
        else:
            self.hashed_data = read_h5_hashes(infile)

class CatalogDB(object):
    """
    SQLite Database of the star catalog indices.

    Attributes:
        - db_path (str): Path to the database file.
    Methods:
        - add_table: Adds a new table to the database.
        - del_table: Deletes an existing table from the database.
        - table_list: Lists all tables in the database.
    """   

    def __init__(self,db_path):
        """
        Initializes a new instance of the CatalogDB class.

        Inputs:
            db_path -> [str] Path to the database file.
        """ 

        self.db_path = db_path

    def __repr__(self):

        """
        String representation of the CatalogDB instance.

        Returns:
            A formatted string that provides a summary of the CatalogDB instance.
        """

        tables = self.table_list(return_output=True)
        return f"<CatalogDB object: tables={tables}>"

    def add_table(self,table_path):
        """
        Adds a new table to the database from the specified table path.

        Inputs:
            table_path -> [str] Path to the CSV file containing the table data to be added.
        """
        generate_catalog_db(self.db_path,table_path)

    def del_table(self,tb_name):
        """
        Deletes an existing table from the database.

        Inputs:
            tb_name -> [str] Name of the table to be deleted from the database.
        """
        delete_table(self.db_path,tb_name)

    def table_list(self, return_output=False):
        """
        Lists all tables in the database.

        Inputs:
            return_output -> [bool, optional, default=False] If True, returns the list of table names instead of printing them.
        Outputs:
            If return_output is True, returns a list of table names. Otherwise, prints the names of all tables in the database. If no tables exist, prints a message indicating so.
        """
        cursor = sqlite3.connect(self.db_path).cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = cursor.fetchall()
        tb_names = [table[0] for table in tables]

        if return_output:
            return tb_names

        if tables:
            print("Tables in the database:")
            for table in tb_names:
                print(table)
        else:
            print("No tables in the database.")

class Stars(object):
    """
    A class to represent a collection of stars.

    Attributes:
        - sc_name -> [str] Name of the star catalog from which the data is derived.
        - center -> [tuple of float] Center of the region of interest in the sky, expressed in [RA, DEC] coordinates (degrees).
        - df -> [pandas.DataFrame] A pandas DataFrame containing detailed information about each star.
        - num_stars -> [int] Total number of the stars.
        - xy -> [2D array-like] Pixel coordinates of the stars.
        - radec -> [2D array-like] Celestial coordinates (Right Ascension and Declination) of the stars.

    Methods:
        - pixel_xy: Calculates pixel coordinates of the stars.
        - invariantfeatures: Calculates the geometric invariant features for stars.
    """   

    def __init__(self,info):
        """
        Initializes a new instance of the Stars class.

        Inputs:
            info -> [dict] A dictionary containing key-value pairs for the instance.
        """
        self.__dict__.update(info)

    def __repr__(self):
        """
        String representation of the Stars instance, typically used for debugging and logging.

        Returns:
            A formatted string that provides a summary of the Stars instance.
        """
    
        return "<Stars object: CATALOG = '{:s}' STARS_NUM = {:d} MCP = {:d}>".format(
            self.sc_name,self.num_stars,self.max_control_points)

    def pixel_xy(self,pixel_width,theta=0):
        """
        Calculates pixel coordinates of the stars in the collection.

        Usage:
            >>> stars = at_hyg24_simplified_stars.pixel_xy(0.01)
        Inputs:
            pixel_width -> [float] Pixel width in degrees, used to convert celestial coordinates to pixel coordinates.
            theta -> [float,optional,default=0] Rotation angle (in radians) to align the WCS frame(equivalent to ENU) with the image reference frame.
        Outputs:
            Stars instance with updated pixel coordinates.
        """

        # Convert celestial coordinates to pixel coordinates of stars
        x,y,wcs = xy_catalog(self.center,self.radec,pixel_width,theta)

        df = self.df
        df['pixelx'],df['pixely'] = x,y
        self.xy = np.stack([x,y]).T
        self.wcs = wcs
        return self

    def invariantfeatures(self,mode_invariants='triangles'):
        """
        Generates geometric invariant features based on the spatial configuration of stars, aiding in star pattern recognition.

        Inputs:
            mode_invariants -> [str, optional, default='triangles'] Mode of geometric invariants to use, e.g., 'triangles' or 'quads'.
        Outputs:
            Stars instance with updated invariant features and KDTree structure for efficient spatial queries.
        Steps:
            1. Derive unique geometric invariants for every possible triangle formed by groups of three stars or quads formed by groups of four stars.
            2. Construct a KDTree structure using these unique invariants, facilitating efficient spatial queries and pattern matching.
            3. Associate each set of invariants with indices of the stars forming the corresponding triangles or quads.
        Notes:
            This method should be called after pixel coordinates are calculated.
            The geometric invariants features are invariant to rotation, scaling, and translation.
        """
        if not hasattr(self,'xy'): 
            raise Exception("The pixel coordinates of stars should be calculated first by `.pixel_xy(pixel_width`)")
        inv_uniq,vrtx_uniq,inv_uniq_tree = calculate_invariantfeatures(self.xy, mode_invariants)
        self.invariants,self.asterisms,self.kdtree = inv_uniq,vrtx_uniq,inv_uniq_tree
        return self

    def tiles_draw(self):
        """
        Visualize the scope of the search area and the coverage of the corresponding tiles.

        Usage:
            >>> at_hyg24_simplified.tiles_draw()
        Outputs:
            An image illustrating the scope of the search area and the coverage of corresponding tiles.         
        """
        search_draw(self.nside,self.tiles_ids,self.search_area,self.radec,self._fov_min)
