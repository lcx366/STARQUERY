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

from .catalog_download import catalog_download,hyg_download
from .catalog_check import catalog_check
from .utils.starcatalog_statistic import tiles_statistic,starcatalog_info
from .utils.df2info import df2info
from .catalog_query import search_box,search_cone,search_box_magpm,search_cone_magpm,box2seqs,cone2seqs,_load_files
from .tiles_draw import search_draw
from .wcs import xy_catalog
from .invariantfeatures import unique_triangles

def solidangle_ratio(fov, r):
    """
    Calculate the ratio of solid angles for a square and a spherical cap, both centered at the same point on the sphere.

    Usage:
        >>> ratio = Solid_angles_ratio(8,10)
    Inputs:
        fov -> [float] Field of view of a camera, defining a square on the sphere, in degrees.
        r -> [float] Angular radius defining a spherical cap, in degrees.
    Outputs:
        ratio -> [float] Ratio of the solid angles of the square to the spherical cap.
    """
    fov_rad = np.deg2rad(fov) / 2
    r_rad = np.deg2rad(r)
    # Solid angle of the square
    Omega_square = 4*np.arcsin(np.tan(fov_rad/2)**2)
    # Solid angle of the spherical cap
    Omega_cone = 4*np.pi*np.sin(r_rad/2)**2
    return Omega_square / Omega_cap

class StarCatalog(object):
    """
    StarCatalog class for handling star catalogs with varying levels of details, including 'raw', 'reduced', and 'simplified' versions.
    Types of Catalogs:
        - 'raw': This is the original star catalog, encompassing comprehensive information about the stars, including their celestial coordinates, magnitudes, proper motions, and other astrophysical parameters.
        - 'reduced': The reduced star catalog retains key data such as the position, proper motion, apparent magnitude, and epoch of the stars, making it a more concise version of the raw catalog.
        - 'simplified': This version of the catalog is streamlined to include only essential information like the position and apparent magnitude of stars at a specific epoch. It is the most minimalist representation of star data in the catalog series.

    Methods:
        - get: Downloads star catalog data files from remote servers or loads them from a specified directory.
        - load: Load star catalog data from local storage.
        - read_h5_indices: Reads an h5-formatted star catalog index file generated using the healpix algorithm, which is essential for blind star pattern matching.
    """
    @staticmethod
    def get(sc_name,tile_size=2,dir_to=None):
        """
        Downloads star catalog data files from remote servers or loads them from a specified directory.

        Usage:
            >>> from starcatalogquery import StarCatalog
            >>> hygv37_raw = StarCatalog.get('hygv3.7',5)
        Inputs:
            sc_name -> [str] Name of the star catalog (e.g., 'hygv3.7', 'at-hygv2.4', 'gaiadr3'). 
            Available options include 'hygv3.7', 'at-hygv2.4', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc. 
            For more star catalogs, please refer to the Space Telescope Science Institute:
            https://outerspace.stsci.edu/display/GC/WebServices+for+Catalog+Access
            tile_size -> [int, optional, default=2] Size of the tile in degrees.
            dir_to -> [str, optional, default=None] Path to save the downloaded catalog. If None, it is assigned to a build-in directory by default.
        Outputs:
            sc_raw -> [StarCatalogRaw object] An instance of the StarCatalogRaw class with downloaded data.
        """

        # Handle specific catalogs that require different download methods
        if sc_name in ['hygv3.7', 'at-hygv2.4']:
            # Use a specialized method for HYG databases
            dir_to, dir_size, file_num, validity, tile_size = hyg_download(sc_name, tile_size, dir_to)

        elif sc_name in ['gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob']:
            # General case for other star catalogs
            if dir_to is None:
                dir_to = 'starcatalogs/raw/{:s}/res{:d}/'.format(sc_name, tile_size)

            if os.path.exists(dir_to):
                # Check existing catalog for validity
                dir_to,dir_size,file_num,validity = catalog_check(sc_name,tile_size,dir_to)
            else:
                # Download and check catalog if not already present
                dir_to,tile_size = catalog_download(sc_name,tile_size,dir_to)  
                dir_to,dir_size,file_num,validity = catalog_check(sc_name,tile_size,dir_to) 
        else:
            raise Exception("Star catalog '{:s}' is not supported.".format(sc_name))        

        # Gather catalog information
        stars_num,mag,description = starcatalog_info(sc_name)
        dict_values = [dir_to,dir_size,file_num,validity,sc_name,f'{tile_size} deg','raw',stars_num,mag,description]
        dict_keys = ['tiles_dir','sc_size','tiles_num','validity','sc_name','tile_size','_mode','stars_num','mag','description']
        info = dict(zip(dict_keys, dict_values))

        # Return an instance of the raw star catalog
        return StarCatalogRaw(info)  

    @staticmethod
    def load(dir_from=None):
        """
        Load star catalog data from local storage.

        Usage:
            >>> from starcatalogquery import StarCatalog
            >>> # load the raw star catalog HYG v3.7
            >>> dir_from_raw = '/Volumes/TOSHIBA/starcatalogs/raw/hygv3.7/res5/'
            >>> hygv37_raw = StarCatalog.load(dir_from_raw)
            >>>
            >>> # load the reduced star catalog HYG v3.7
            >>> dir_from_reduced = '/Volumes/TOSHIBA/starcatalogs/reduced/hygv3.7/res5/'
            >>> hygv37_reduced = StarCatalog.load(dir_from_reduced)
            >>>
            >>> # load the simplified star catalog HYG v3.7
            >>> dir_from_simplified = '/Volumes/TOSHIBA/starcatalogs/simplified/hygv3.7/res5/mag9.0/epoch2022.0/'
            >>> hygv37_simplified = StarCatalog.load(dir_from_simplified)
        Inputs:
            dir_from -> [str, optional, default=None] Path of the star catalog files. If None, a default directory is used.
        Outputs:
            StarCatalogRaw | StarCatalogReduced | StarCatalogSimplified: 
            An instance of the appropriate star catalog class, depending on the catalog type identified in the directory.
        """
        # Extract catalog type, name, and tile size from the directory path
        _mode,sc_name,tile_size = dir_from.split('starcatalogs/')[1].split('/')[:3]
        tile_size = int(tile_size[3:]) # Extracting numerical value from the tile size string

        # Based on the levels of details of catalogs, load the corresponding star catalog
        if _mode == 'raw':
            starcatalog = StarCatalogRaw.load(sc_name,tile_size,dir_from)
        elif _mode == 'reduced':
            starcatalog = StarCatalogReduced.load(sc_name,tile_size,dir_from)
        elif _mode == 'simplified':    
            mag_threshold,epoch = np.array(dir_from.split('mag')[1][:-1].split('/epoch'),dtype=float) 
            starcatalog = StarCatalogSimplified.load(sc_name,tile_size,mag_threshold,epoch,dir_from)

        return starcatalog 

    @staticmethod
    def read_h5_indices(infile):
        """
        Reads an h5-formatted star catalog index file generated using the healpix algorithm, which is essential for blind star pattern matching.
        This method extracts information from the file including center coordinates (Ra, Dec) of each sky area, 
        pixel coordinates of stars within these areas, triangle invariant features, and indices of stars forming these triangles.

        Usage:
            >>> from starcatalogquery import StarCatalog
            >>> infile_h5 = 'starcatalogs/indices/hygv3.7/k2_mag9.0_mcp30_2022.0.h5'
            >>> fp_radecs,stars_xy,stars_invariants,stars_asterisms = read_h5_indices(infile_h5)
        Inputs:
            infile -> [str] Path of the h5-formatted star catalog indices file.
        Outputs:
            fp_radecs -> [tuple] central coordinates of sky areas in form of [RA, DEC] in degrees.
            stars_xy -> [2D array-like,float] Pixel coordinates of stars in a sky area.
            stars_invariants -> [2d array-like,float] Triangle invariants for a sky area.
            stars_asterisms -> [2d array-like,int] Indices of stars forming these triangles.
        """
        with h5py.File(infile, 'r') as fin:
            # Extract data from the file
            fp_radecs = fin['fp_radecs'][:]
            stars_xy,stars_invariants, stars_asterisms = [],[],[]
            sort_index = np.argsort(np.abs(fp_radecs[:,1]))
            fp_radecs_sort = fp_radecs[sort_index]

            for j in sort_index:
                stars_xy.append(fin['stars_xy/'+str(j)][:])
                stars_invariants.append(fin['stars_invariants/'+str(j)][:])
                stars_asterisms.append(fin['stars_asterisms/'+str(j)][:])
    
        return fp_radecs_sort,stars_xy,stars_invariants,stars_asterisms  

class StarCatalogRaw(object):
    """
    This version of the catalog is the original star catalog, encompassing comprehensive information about the stars, including their celestial coordinates, magnitudes, proper motions, and other astrophysical parameters.

    Attributes:
        tiles_dir (str): Path where the star catalog files are stored.
        tiles_num (int): Total number of tile files in the catalog.
        tile_size (str): Size of each tile in the catalog, in degrees.
        sc_size (str): The overall size of the star catalog (e.g., in MB or GB).
        validity (bool): Indicates whether the star catalog is complete and valid for use.
        sc_name (str): The name of the star catalog, such as 'hygv3.7'.
        _mode (str): The mode of the catalog, set to 'raw' for this class.
        stars_num (str): The total number of stars included in the catalog.
        mag (str): The range of magnitudes of stars in the catalog.
        description (str): A brief description or summary of the catalog.

    Methods:
        load: Loads the raw star catalog from a specified directory.
        reduce: Reduces the raw star catalog to contain only essential information.
        search_box: Searches for stars within a specified rectangular area.
        search_cone: Searches for stars within a specified conical area.
        _search_draw: Visualizes the search area on a map and shows the corresponding tiles.  
    """   
    def __init__(self,info): 
        """
        Initializes a new instance of the StarCatalogRaw class.

        Inputs:
            info -> [dict] A dictionary containing key-value pairs of catalog attributes.
        """ 

        for key in info.keys():
            setattr(self, key, info[key])

    def __repr__(self):
        """
        Returns a string representation of the StarCatalogRaw instance.

        Returns:
            A string describing the StarCatalogRaw instance.
        """
    
        return "<StarCatalogRaw object: CATALOG_NAME = '{:s}' CATALOG_SIZE = '{:s}' TILES_NUM = {:d} TILE_SIZE = '{:s}' STARS_NUM = '{:s}' MAG = '{:s}'>".format(
            self.sc_name,self.sc_size,self.tiles_num,self.tile_size,self.stars_num,self.mag)

    @classmethod
    def load(cls, sc_name, tile_size, dir_from=None):
        """
        Load the raw star catalog files from a specified directory.

        Usage:
            >>> from starcatalogquery import StarCatalogRaw
            >>> # load the raw star catalog HYG v3.7
            >>> dir_from_raw = '/Volumes/TOSHIBA/starcatalogs/raw/hygv3.7/res5/'
            >>> hygv37_raw = StarCatalogRaw.load('hygv3.7',5,dir_from_raw)
        Inputs:
            sc_name -> [str] Name of the star catalog (e.g., 'hygv3.7').
            tile_size -> [int] Size of the tile in degrees.
            dir_from -> [str, optional, default=None] Directory containing the star catalog files. If None, defaults to a built-in directory.
        Outputs:
            hygv37_raw -> [StarCatalogRaw object] Instance of the class StarCatalogRaw with loaded data.
        Raises:
            Exception: If the specified directory does not exist.
        """
        # Assign default directory if none is provided
        if dir_from is None: dir_from = f'starcatalogs/raw/{sc_name}/res{tile_size}/'

        # Check if the directory exists, raise an exception if it does not
        if not os.path.exists(dir_from): raise Exception(f'Path of the star catalog {sc_name} does not exist.')

        # Calculate total size and number of tile files
        file_num, dir_size, validity = tiles_statistic(dir_from, tile_size)
        stars_num, mag, description = starcatalog_info(sc_name)

        info_keys = ['tiles_dir', 'sc_size', 'tiles_num', 'validity', 'sc_name', 'tile_size', '_mode', 'stars_num', 'mag', 'description']
        info_values = [dir_from, dir_size, file_num, validity, sc_name, f'{tile_size} deg', 'raw', stars_num, mag, description]
        info = dict(zip(info_keys, info_values))

        return cls(info)

    def reduce(self,dir_reduced=None):
        """
        Reduce the raw star catalog, retaining only essential star information.

        Usage:
            >>> from starcatalogquery import StarCatalogRaw
            >>> dir_from_raw = '/Volumes/TOSHIBA/starcatalogs/raw/hygv3.7/res5/'
            >>> hygv37_raw = StarCatalogRaw.load('hygv3.7',5,dir_from_raw)
            >>> hygv37_reduced = hygv37_raw.reduce()
        Inputs:
            dir_reduced -> [str, optional, default=None] Directory to store the reduced star catalog. Defaults to a built-in directory.
        Outputs:
            hygv37_reduced -> [StarCatalogReduced object] Instance of the reduced star catalog.
        Note:
            The reduction process simplifies the catalog by retaining only key details such as position, proper motion, and magnitude.
        """
        # Copying existing information
        info = self.__dict__.copy()

        # Extracting the tile size
        tile_size = int(self.tile_size.split(' ')[0])

        # Setting up the directory for the reduced catalog
        tiles_dir = self.tiles_dir
        if dir_reduced is None:
            dir_reduced = tiles_dir.replace('raw','reduced')
        Path(dir_reduced).mkdir(parents=True, exist_ok=True)    

        file_list = glob(tiles_dir+'*') 
        sc_name = self.sc_name 
        print(f'Reducing the star catalog {sc_name}, which may take a considerable amount of time')

        if sc_name in ['hygv3.7','at-hygv2.4']:
            # Processing each file in the raw star catalog directory
            for j, tile_file in enumerate(file_list, start=1):
                desc = f'Reducing {Fore.BLUE}{j}{Fore.RESET} of {self.tiles_num}'
                print(desc,end='\r')
                df = pd.read_csv(tile_file,skiprows=1,dtype=str)

                # Applying specific transformations based on the catalog name
                # units: ra->hourangle, dec->deg, pmra->mas/a, pmdec->mas/a, epoch->2000.0 
                if sc_name == 'at-hygv2.4':
                    columns_dict = {'pm_ra':'pmra', 'pm_dec':'pmdec'}
                    df.rename(columns=columns_dict, inplace=True)
                df_reduced = df.loc[:,['ra','dec','pmra','pmdec','mag']]
                df_reduced['epoch'] = '2000.0'
                df_reduced['ra'] = (df_reduced['ra'].astype(float)*15).round(6) # Convert hourangle to deg
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)
        elif sc_name == 'gsc12':
            for j, tile_file in enumerate(file_list, start=1):
                desc = f'Reducing {Fore.BLUE}{j}{Fore.RESET} of {self.tiles_num}'
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
        elif sc_name == 'gsc242':
            for j, tile_file in enumerate(file_list, start=1):
                desc = f'Reducing {Fore.BLUE}{j}{Fore.RESET} of {self.tiles_num}'
                print(desc,end='\r')
                df = pd.read_csv(tile_file,skiprows=1,dtype=str)
                # units: ra->deg, dec->deg, rapm->mas/a, decpm->mas/a, epoch->2012
                df_reduced = df.loc[:,['ra','dec','rapm','decpm','mag','epoch']]
                columns_dict = {'rapm':'pmra', 'decpm':'pmdec'}
                df_reduced.rename(columns=columns_dict, inplace=True)
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)  
        elif sc_name == 'gaiadr3':
            for j, tile_file in enumerate(file_list, start=1):  
                desc = f'Reducing {Fore.BLUE}{j}{Fore.RESET} of {self.tiles_num}'
                print(desc,end='\r')      
                df = pd.read_csv(tile_file,skiprows=1,dtype=str)
                # units: ra->deg, dec->deg, pmra->mas/a, pmdec->mas/a, epoch->2016
                df_reduced = df.loc[:,['ra','dec','pmra','pmdec','mag','epoch']]
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)           
        elif sc_name == 'ucac5':
            for j, tile_file in enumerate(file_list, start=1):        
                desc = f'Reducing {Fore.BLUE}{j}{Fore.RESET} of {self.tiles_num}'
                print(desc,end='\r')
                df = pd.read_csv(tile_file,skiprows=1,dtype=str)
                # units: ra->deg, dec->deg, pmur->mas/a, pmud->mas/a, epu->1998.754
                df_reduced = df.loc[:,['ra','dec','pmur','pmud','mag','epu']]
                columns_dict = {'pmur':'pmra', 'pmud':'pmdec','epu':'epoch'}
                df_reduced.rename(columns=columns_dict, inplace=True)
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)   
        elif sc_name == 'usnob':
            for j, tile_file in enumerate(file_list, start=1):        
                desc = f'Reducing {Fore.BLUE}{j}{Fore.RESET} of {self.tiles_num}'
                print(desc,end='\r')
                df = pd.read_csv(tile_file,skiprows=1,dtype=str)
                # units: ra->deg, dec->deg, pmRA->mas/a, pmDEC->mas/a, Epoch->1950,mag->unknown
                df_reduced = df.loc[:,['ra','dec','pmRA','pmDEC','mag','Epoch']]
                columns_dict = {'pmRA':'pmra', 'pmDEC':'pmdec','Epoch':'epoch'}
                df_reduced.rename(columns=columns_dict, inplace=True)
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)  
        elif sc_name == '2mass':
            for j, tile_file in enumerate(file_list, start=1):   
                desc = f'Reducing {Fore.BLUE}{j}{Fore.RESET} of {self.tiles_num}'
                print(desc,end='\r')
                df = pd.read_csv(tile_file,skiprows=1,dtype=str)
                # units: ra->deg, dec->deg, jdate->2451063.6417
                df_reduced = df.loc[:,['ra','dec','mag']]
                df_reduced['epoch'] = Time(df['jdate'].astype(float), format='jd').jyear.round(2) 
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)

        print('\nFinished')                                                   
                    
        file_num,dir_size,validity = tiles_statistic(dir_reduced,tile_size)  
        info['_mode'] = 'reduced'
        info['tiles_dir'] = dir_reduced
        info['sc_size'] = dir_size

        return StarCatalogReduced(info)  
    
    def search_box(self,radec_box,mag_threshold,t_pm,max_num=None,max_num_per_tile=None):
        """
        Perform a rectangle search of stars on raw star catalogs.

        Usage:
            >>> from starcatalogquery import StarCatalogRaw
            >>> dir_from_raw = '/Volumes/TOSHIBA/starcatalogs/raw/hygv3.7/res5/'
            >>> hygv37_raw = StarCatalogRaw.load('hygv3.7', 5, dir_from_raw)
            >>> stars = hygv37_raw.search_box([20, 30, 30, 40], 9, 2022.0)
        Inputs:
            radec_box -> [list] Rectangular search area in form of [ra_min, dec_min, ra_max, dec_max] in degrees.
            mag_threshold -> [float] Apparent magnitude limit.
            t_pm -> [float] Epoch to which the stars are unified.
            max_num -> [int,optional,default=None] Maximum number of stars to include in the search, sorted by brightness. If None, includes all stars meeting the magnitude criteria.
        Outputs:
            stars -> [Stars object] An instance of Stars with the search results.
        """
        # Calculating the center of the search box
        ra_min,dec_min,ra_max,dec_max = radec_box
        center = [(ra_min+ra_max)/2,(dec_min+dec_max)/2]

        # Determining the tile size and preparing to load catalog data
        tile_size = int(self.tile_size.split()[0])
        sc_path,sc_name = self.tiles_dir,self.sc_name

        # Calculate indices of star catalog tiles covering the rectangular search area.
        sc_indices = box2seqs(radec_box,tile_size) 

        # Loading data from relevant tiles and concatenating them into a DataFrame
        df = pd.concat(_load_files(sc_indices,sc_path,sc_name,self._mode,max_num_per_tile))

        if sc_name in ['hygv3.7','at-hygv2.4']:
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

        # Correct the proper motion
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

        info = df2info(self.sc_name,center,df,max_num,radec_box,'BOX') 
        return Stars(info)

    def search_cone(self,center,radius,mag_threshold,t_pm,max_num=None,max_num_per_tile=None):
        """
        Perform a cone search of stars on raw star catalogs.

        Usage:
            >>> from starcatalogquery import StarCatalogRaw
            >>> dir_from_raw = '/Volumes/TOSHIBA/starcatalogs/raw/hygv3.7/res5/'
            >>> hygv37_raw = StarCatalogRaw.load('hygv3.7', 5, dir_from_raw)
            >>> stars = hygv37_raw.search_cone([20, 30], 10, 9, 2022.0)
        Inputs:
            center -> [list] Center of the cone in form of [ra_c, dec_c] in degrees.
            radius -> [float] Angular radius of the cone.
            mag_threshold -> [float] Apparent magnitude limit.
            t_pm -> [float] Epoch to which the stars are unified.
            max_num -> [int,optional,default=None] Maximum number of stars to include in the search, sorted by brightness. If None, includes all stars meeting the magnitude criteria.
        Outputs:
            stars -> [Stars object] An instance of Stars with the search results.
        """
        # Extracting the center coordinates of the cone search
        ra_c,dec_c = center

        # Determining the tile size and preparing to load catalog data
        tile_size = int(self.tile_size.split()[0])
        sc_path,sc_name = self.tiles_dir,self.sc_name

        # Calculate indices of star catalog tiles covering the cone search area.
        sc_indices = cone2seqs(ra_c,dec_c,radius,tile_size) 

        # Loading data from relevant tiles and concatenating them into a DataFrame
        df = pd.concat(_load_files(sc_indices,sc_path,sc_name,self._mode,max_num_per_tile))

        if sc_name in ['hygv3.7','at-hygv2.4']:
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

        # Correct the proper motion
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

        info = df2info(self.sc_name,center,df,max_num,radius,'CONE') 
        return Stars(info)

    def _search_draw(self,search_area):
        """
        Visualize the scope of the search area and the coverage of corresponding tiles.

        Usage:
            >>> from starcatalogquery import StarCatalogRaw
            >>> dir_from_raw = '/Volumes/TOSHIBA/starcatalogs/raw/hygv3.7/res5/'
            >>> hygv37_raw = StarCatalogRaw.load('hygv3.7', 5, dir_from_raw)
            >>> cone_area = {'cone': [20, 30, 10]}
            >>> hygv37_raw._search_draw(cone_area)  
        Inputs:
            search_area -> [dict] Scope of the search area, in form of {'cone': [ra_c, dec_c, radius]} or {'box': [ra_min, dec_min, ra_max, dec_max]}, where
            where ra_min, ra_max, ra_c are Right Ascension and dec_min, dec_max, dec_c are Declination, all in degrees.
        Outputs:
            An image illustrating the scope of the search area and the coverage of corresponding tiles.
        """
        # Retrieve the tile size
        tile_size = int(self.tile_size.split()[0]) 
        search_draw(tile_size,search_area)       

class StarCatalogReduced(object):
    """
    This version of the catalog represents the reduced star catalog, which retains key data such as the position, proper motion, apparent magnitude, and epoch of the stars, making it a more concise version of the raw catalog.

    Attributes:
        tiles_dir (str): Path where the star catalog files are stored.
        tiles_num (int): Total number of tile files in the catalog.
        tile_size (str): Size of each tile in the catalog, in degrees.
        sc_size (str): The overall size of the star catalog (e.g., in MB or GB).
        validity (bool): Indicates whether the star catalog is complete and valid for use.
        sc_name (str): The name of the star catalog, such as 'hygv3.7'.
        _mode (str): The mode of the catalog, set to 'reduced' for this class.
        stars_num (str): The total number of stars included in the catalog.
        mag (str): The range of magnitudes of stars in the catalog.
        description (str): A brief description or summary of the catalog.

    Methods:
        load: Loads the reduced star catalog from a specified directory.
        simplify: Further refine the star catalog, resulting in a StarCatalogSimplified instance. 
        This process strips down the catalog to its most essential elements, typically including only the celestial positions and apparent magnitudes of stars at a specified epoch.
        search_box: Searches for stars within a specified rectangular area.
        search_cone: Searches for stars within a specified conical area.
        _search_draw: Visualizes the search area on a map and shows the corresponding tiles.  
    """ 
    def __init__(self,info):
        """
        Initializes a new instance of the StarCatalogReduced class.

        Inputs:
            info -> [dict] A dictionary containing key-value pairs of catalog attributes.
        """ 

        for key in info.keys():
            setattr(self, key, info[key])

    def __repr__(self):
        """
        Returns a string representation of the StarCatalogReduced instance.

        Returns:
            A string describing the StarCatalogReduced instance.
        """
    
        return "<StarCatalogReduced object: CATALOG_NAME = '{:s}' CATALOG_SIZE = '{:s}' TILES_NUM = {:d} TILE_SIZE = '{:s}' STARS_NUM = '{:s}' MAG = '{:s}'>".format(
            self.sc_name,self.sc_size,self.tiles_num,self.tile_size,self.stars_num,self.mag)

    @classmethod
    def load(cls,sc_name,tile_size,dir_from=None):
        """
        Load the reduced star catalog files from a specified directory.

        Usage:
            >>> from starcatalogquery import StarCatalogReduced
            >>> # load the reduced star catalog HYG v3.7
            >>> dir_from_reduced = '/Volumes/TOSHIBA/starcatalogs/reduced/hygv3.7/res5/'
            >>> hygv37_reduced = StarCatalogReduced.load('hygv3.7',5,dir_from_reduced)
        Inputs:
            sc_name -> [str] Name of the star catalog (e.g., 'hygv3.7').
            tile_size -> [int] Size of the tile in degrees.
            dir_from -> [str, optional, default=None] Directory containing the star catalog files. If None, defaults to a built-in directory.
        Outputs:
            hygv37_reduced -> [StarCatalogReduced object] Instance of the class StarCatalogReduced with loaded data.
        Raises:
            Exception: If the specified directory does not exist.
        """
        # Assign default directory if none is provided
        if dir_from is None: 
            dir_from = 'starcatalogs/reduced/{:s}/res{:d}/'.format(sc_name,tile_size)   

        # Check if the directory exists, raise an exception if it does not 
        if not os.path.exists(dir_from): raise Exception('Path of the star catalog {:s} does not exist.'.format(sc_name))  

        # Calculate total size and number of tile files
        file_num,dir_size,validity = tiles_statistic(dir_from,tile_size) 
        stars_num,mag,description = starcatalog_info(sc_name)

        dict_values = dir_from,dir_size,file_num,validity,sc_name,'{:d} deg'.format(tile_size),'reduced',stars_num,mag,description
        dict_keys = 'tiles_dir','sc_size','tiles_num','validity','sc_name','tile_size','_mode','stars_num','mag','description'
        info = dict(zip(dict_keys, dict_values))

        return cls(info)  

    def simplify(self,mag_threshold,t_pm,dir_simplified=None):
        """
        Simplify the reduced star catalog, making it a minimalist of the raw catalog.
        
        Usage:
            >>> from starcatalogquery import StarCatalogReduced
            >>> # load the reduced star catalog HYG v3.7
            >>> dir_from_reduced = '/Volumes/TOSHIBA/starcatalogs/reduced/hygv3.7/res5/'
            >>> hygv37_reduced = StarCatalogReduced.load('hygv3.7',5,dir_from_reduced)
            >>> hygv37_simplified = hygv37_reduced.simplify(9.0,2022.0)
        Inputs:
            mag_threshold -> [float] Apparent magnitude limit.
            t_pm -> [float] Epoch to which the stars are unified.
            dir_simplified -> [str, optional, default=None] Directory containing the star catalog files. If None, defaults to a built-in directory.
        Outputs:
            Instance of class StarCatalogSimplified      
        """
        # Copying existing information
        info = self.__dict__.copy()

        # Extracting the tile size
        tile_size = int(self.tile_size.split(' ')[0])

        # Setting up the directory for the reduced catalog
        tiles_dir = self.tiles_dir
        if dir_simplified is None:
            dir_simplified = tiles_dir.replace('reduced','simplified')+'mag{:.1f}/epoch{:.1f}/'.format(mag_threshold,t_pm)
        Path(dir_simplified).mkdir(parents=True, exist_ok=True)    

        file_list = glob(tiles_dir+'*') 
        sc_name = self.sc_name
        print('Simplifying the star catalog {:s}, which may take a considerable amount of time'.format(sc_name))  

        # Processing each file in the raw star catalog directory
        for j, tile_file in enumerate(file_list, start=1):
            desc = 'Simplifying {:s}{:d}{:s} of {:d}'.format(Fore.BLUE,j,Fore.RESET,self.tiles_num)
            print(desc,end='\r')

            df_simplified = pd.read_csv(tile_file)
            mag = df_simplified['mag']
            mag_flag = mag < mag_threshold
            df_simplified = df_simplified[mag_flag].sort_values(by=['mag'])    

            # Correct proper motion    
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

        print('\nFinished')

        file_num,dir_size,validity = tiles_statistic(dir_simplified,tile_size)  
        info['_mode'] = 'simplified'
        info['tiles_dir'] = dir_simplified
        info['sc_size'] = dir_size
        info['mag_threshold'] = mag_threshold
        info['epoch'] = t_pm

        return StarCatalogSimplified(info)    

    def search_box(self,radec_box,mag_threshold,t_pm,max_num=None,max_num_per_tile=None):
        """
        Perform a rectangle search of stars on reduced star catalogs.

        Usage:
            >>> from starcatalogquery import StarCatalogReduced
            >>> # load the reduced star catalog HYG v3.7
            >>> dir_from_reduced = '/Volumes/TOSHIBA/starcatalogs/reduced/hygv3.7/res5/'
            >>> hygv37_reduced = StarCatalogReduced.load('hygv3.7',5,dir_from_reduced)
            >>> stars = hygv37_reduced.search_box([20,30,30,40],9,2022.0)
        Inputs:
            radec_box -> [list] Rectangular search area in form of [ra_min, dec_min, ra_max, dec_max] in degrees.
            mag_threshold -> [float] Apparent magnitude limit.
            t_pm -> [float] Epoch to which the stars are unified.
            max_num -> [int,optional,default=None] Maximum number of stars to include in the search, sorted by brightness. If None, includes all stars meeting the magnitude criteria.
        Outputs:
            stars -> [Stars object] An instance of Stars with the search results.
        """
        # Calculating the center of the search box
        ra_min,dec_min,ra_max,dec_max = radec_box
        center = [(ra_min+ra_max)/2,(dec_min+dec_max)/2]

        # Determining the tile size and preparing to load catalog data
        tile_size = int(self.tile_size.split()[0])

        # Performs a rectangular search on the reduced star catalog considering magnitude and proper motion.
        df = search_box_magpm(radec_box,self.tiles_dir,self.sc_name,tile_size,mag_threshold,t_pm,max_num_per_tile)
        info = df2info(self.sc_name,center,df,max_num,radec_box,'BOX') 
        return Stars(info)

    def search_cone(self,center,radius,mag_threshold,t_pm,max_num=None,max_num_per_tile=None):   
        """
        Perform a cone search of stars on reduced star catalogs.

        Usage:
            >>> from starcatalogquery import StarCatalogReduced
            >>> # load the reduced star catalog HYG v3.7
            >>> dir_from_reduced = '/Volumes/TOSHIBA/starcatalogs/reduced/hygv3.7/res5/'
            >>> hygv37_reduced = StarCatalogReduced.load('hygv3.7',5,dir_from_reduced)
            >>> stars = hygv37_reduced.search_cone([20,30],10,9,2022.0)
        Inputs:
            center -> [list] Center of the cone in form of [ra_c, dec_c] in degrees.
            radius -> [float] Angular radius of the cone.
            mag_threshold -> [float] Apparent magnitude limit.
            t_pm -> [float] Epoch to which the stars are unified.
            max_num -> [int,optional,default=None] Maximum number of stars to include in the search, sorted by brightness. If None, includes all stars meeting the magnitude criteria.
        Outputs:
            stars -> [Stars object] Instance of class Stars containing search results.
        """
        # Determining the tile size
        tile_size = int(self.tile_size.split()[0])

        # Performs a conical search of stars on the reduced star catalog considering magnitude and proper motion.
        df = search_cone_magpm(center,radius,self.tiles_dir,self.sc_name,tile_size,mag_threshold,t_pm,max_num_per_tile)
        info = df2info(self.sc_name,center,df,max_num,radius,'CONE') 
        return Stars(info)

    def _search_draw(self,search_area):
        """
        Visualize the scope of the search area and the coverage of corresponding tiles.

        Usage:
            >>> from starcatalogquery import StarCatalogReduced
            >>> # load the reduced star catalog HYG v3.7
            >>> dir_from_reduced = '/Volumes/TOSHIBA/starcatalogs/reduced/hygv3.7/res5/'
            >>> hygv37_reduced = StarCatalogReduced.load('hygv3.7',5,dir_from_reduced)
            >>> cone_area = {'cone':[20,30,10]}
            >>> stars = hygv37_reduced._search_draw(cone_area)
        Inputs:
            search_area -> [dict] Scope of the search area, in form of {'cone': [ra_c, dec_c, radius]} or {'box': [ra_min, dec_min, ra_max, dec_max]}, where
            where ra_min, ra_max, ra_c are Right Ascension and dec_min, dec_max, dec_c are Declination, all in degrees.
        Outputs:
            An image illustrating the scope of the search area and the coverage of corresponding tiles.
        """
        tile_size = int(self.tile_size.split()[0]) 
        search_draw(tile_size,search_area)   

class StarCatalogSimplified(object):
    """
    This version of the catalog is streamlined to include only essential information like the position and apparent magnitude of stars at a specific epoch. 
    It is the most minimalist representation of star data in the catalog series.

    Attributes:
        tiles_dir (str): Path where the star catalog files are stored.
        tiles_num (int): Total number of tile files in the catalog.
        tile_size (str): Size of each tile in the catalog, in degrees.
        sc_size (str): The overall size of the star catalog (e.g., in MB or GB).
        validity (bool): Indicates whether the star catalog is complete and valid for use.
        sc_name (str): The name of the star catalog, such as 'hygv3.7'.
        _mode (str): The mode of the catalog, set to 'raw' for this class.
        stars_num (str): The total number of stars included in the catalog.
        mag (str): The range of magnitudes of stars in the catalog.
        description (str): A brief description or summary of the catalog.

    Methods:
        load: Loads the raw star catalog from a specified directory.
        search_box: Searches for stars within a specified rectangular area.
        search_cone: Searches for stars within a specified conical area.
        _search_draw: Visualizes the search area on a map and shows the corresponding tiles. 
        - h5_indices: Generate a h5-formatted star catalog indices file.
    """   
    def __init__(self,info):
        """
        Initializes a new instance of the StarCatalogSimplified class.

        Inputs:
            info -> [dict] A dictionary containing key-value pairs of catalog attributes.
        """ 

        for key in info.keys():
            setattr(self, key, info[key])

    def __repr__(self):
        """
        Returns a string representation of the StarCatalogSimplified instance.

        Returns:
            A string describing the StarCatalogSimplified instance.
        """
      
        return "<StarCatalogSimplified object: CATALOG_NAME = '{:s}' CATALOG_SIZE = '{:s}' TILES_NUM = {:d} TILE_SIZE = '{:s}' STARS_NUM = '{:s}' MAG = '{:s}' MAG_CUTOFF = {:.1f} EPOCH = {:.1f}>".format(
            self.sc_name,self.sc_size,self.tiles_num,self.tile_size,self.stars_num,self.mag,self.mag_threshold,self.epoch)   

    @classmethod
    def load(cls,sc_name,tile_size,mag_threshold,epoch,dir_from=None):
        """
        Load the simplified star catalog files from a specified directory.

        Usage:
            >>> from starcatalogquery import StarCatalogSimplified
            >>> # load the simplified star catalog HYG v3.7
            >>> dir_from_simplified = '/Volumes/TOSHIBA/starcatalogs/simplified/hygv3.7/res5/mag9.0/epoch2022.0/'
            >>> hygv37_simplified = StarCatalogSimplified.load('hygv3.7',5,9,2022,dir_from_simplified)
        Inputs:
            sc_name -> [str] Name of the star catalog (e.g., 'hygv3.7').
            tile_size -> [int] Size of the tile in degrees.
            dir_from -> [str, optional, default=None] Directory containing the star catalog files. If None, defaults to a built-in directory.
        Outputs:
            hygv37_simplified -> [StarCatalogSimplified object] Instance of the class StarCatalogSimplified with loaded data.
        Raises:
            Exception: If the specified directory does not exist.
        """
        # Assign default directory if none is provided
        if dir_from is None: 
            print(sc_name,tile_size,mag_threshold,epoch)
            dir_from = 'starcatalogs/simplified/{:s}/res{:d}/mag{:.1f}/epoch{:.1f}/'.format(sc_name,tile_size,mag_threshold,epoch)

        # Check if the directory exists, raise an exception if it does not
        if not os.path.exists(dir_from): raise Exception('Path of the star catalog {:s} does not exist.'.format(sc_name))  

        # Calculate total size and number of tile files   
        file_num,dir_size,validity = tiles_statistic(dir_from,tile_size) 
        stars_num,mag,description = starcatalog_info(sc_name)

        dict_values = dir_from,dir_size,file_num,validity,sc_name,'{:d} deg'.format(tile_size),'simplified',stars_num,mag,description,mag_threshold,epoch
        dict_keys = 'tiles_dir','sc_size','tiles_num','validity','sc_name','tile_size','_mode','stars_num','mag','description','mag_threshold','epoch'
        info = dict(zip(dict_keys, dict_values))

        return cls(info)   

    def search_box(self,radec_box,max_num=None,max_num_per_tile=None):
        """
        Perform a rectangle search of stars on simplified star catalogs.

        Usage:
            >>> from starcatalogquery import StarCatalogSimplified
            >>> # load the simplified star catalog HYG v3.7
            >>> dir_from_simplified = '/Volumes/TOSHIBA/starcatalogs/simplified/hygv3.7/res5/mag9.0/epoch2022.0/'
            >>> hygv37_simplified = StarCatalogSimplified.load('hygv3.7',5,9,2022,dir_from_simplified)
            >>> stars = hygv37_simplified.search_box([20,30,30,40])
        Inputs:
            radec_box -> [list] Rectangular search area in form of [ra_min, dec_min, ra_max, dec_max] in degrees.
            max_num -> [int,optional,default=None] Maximum number of stars to include in the search, sorted by brightness. If None, includes all stars meeting the magnitude criteria.
        Outputs:
            stars -> [Stars object] An instance of Stars with the search results.
        """
        # Calculating the center of the search box
        ra_min,dec_min,ra_max,dec_max = radec_box
        center = [(ra_min+ra_max)/2,(dec_min+dec_max)/2]

        # Determining the tile size
        tile_size = int(self.tile_size.split()[0])

        # Performs a rectangular search on the simplified star catalog without considering magnitude and proper motion.
        df = search_box(radec_box,self.tiles_dir,self.sc_name,tile_size,max_num_per_tile) 
        info = df2info(self.sc_name,center,df,max_num,radec_box,'BOX') 
        return Stars(info)

    def search_cone(self,center,radius,max_num=None,max_num_per_tile=None):   
        """
        Perform a cone search of stars on simplified star catalogs.

        Usage:
            >>> from starcatalogquery import StarCatalogSimplified
            >>> # load the simplified star catalog HYG v3.7
            >>> dir_from_simplified = '/Volumes/TOSHIBA/starcatalogs/simplified/hygv3.7/res5/mag9.0/epoch2022.0/'
            >>> hygv37_simplified = StarCatalogSimplified.load('hygv3.7',5,9,2022,dir_from_simplified)
            >>> stars = hygv37_simplified.search_cone([20,30],10)
        Inputs:
            center -> [list] Center of the cone in form of [ra_c, dec_c] in degrees.
            radius -> [float] Angular radius of the cone.
            max_num -> [int,optional,default=None] Maximum number of stars to include in the search, sorted by brightness. If None, includes all stars meeting the magnitude criteria.
        Outputs:
            stars -> [Stars object] Instance of class Stars containing search results.
        """
        # Extracting the center coordinates of the cone search
        tile_size = int(self.tile_size.split()[0])

        # Performs a cone search on the simplified star catalog without considering magnitude and proper motion.
        df = search_cone(center,radius,self.tiles_dir,self.sc_name,tile_size,max_num_per_tile)
        info = df2info(self.sc_name,center,df,max_num,radius,'CONE') 
        return Stars(info)

    def _search_draw(self,search_area):
        """
        Visualize the scope of the search area and the coverage of the corresponding tiles.

        Usage:
            >>> from starcatalogquery import StarCatalogSimplified
            >>> # load the simplified star catalog HYG v3.7
            >>> dir_from_simplified = '/Volumes/TOSHIBA/starcatalogs/simplified/hygv3.7/res5/mag9.0/epoch2022.0/'
            >>> hygv37_simplified = StarCatalogReduced.load('hygv3.7',5,9,2022,dir_from_simplified)
            >>> cone_area = {'cone':[20,30,10]}
            >>> stars = hygv37_simplified._search_draw(cone_area)
        Inputs:
            search_area -> [dict] Scope of the search area, in form of {'cone': [ra_c, dec_c, radius]} or {'box': [ra_min, dec_min, ra_max, dec_max]}, where
            where ra_min, ra_max, ra_c are Right Ascension and dec_min, dec_max, dec_c are Declination, all in degrees.
        Outputs:
            An image illustrating the scope of the search area and the coverage of corresponding tiles.         
        """
        # Retrieve the tile size
        tile_size = int(self.tile_size.split()[0]) 
        search_draw(tile_size,search_area)   

def h5_indices(self, fov, pixel_width, max_num=30):
    """
    Generate a h5-formatted star catalog indices file.

    Usage:
        >>> hygv37_simplified.h5_indices(8, 0.01)
    Inputs:
        fov -> [float] Field of View (FOV) of the camera in degrees.
        pixel_width -> [float] Pixel width in degrees.
        max_num -> [int, optional, default=30] Maximum number of stars to consider for each sky area.
    Outputs:
        outh5 -> h5-formatted file containing star catalog indices.
        ratio -> [float] Ratio of the solid angles spanned by the square and the cone derived from FOV.
    Note:
        The function divides the sky into equal areas using the healpix algorithm and generates indices for each area.
    """
    # Setting up the directory for h5 file
    dir_h5 = f'starcatalogs/indices/{self.sc_name}/'
    Path(dir_h5).mkdir(parents=True, exist_ok=True)  

    # Determine healpix level and search radius based on FOV
    k, nside, search_radius = self._healpix_level(fov)

    # Calculate the ratio of the solid angles spanned by the square and the cone derived from FOV.
    ratio = solidangle_ratio(fov, search_radius)   

    # Setting up the name for h5 file  
    outh5 = f'{dir_h5}k{k}_mag{self.mag_threshold}_mcp{max_num}_{self.epoch}.h5'
    
    # Check if file already exists
    if os.path.exists(outh5): return outh5, ratio 

    # Writing data to the h5 file
    fout = h5py.File(outh5, 'w')
    stars_xy_grp = fout.create_group("stars_xy")
    stars_invariants_grp = fout.create_group("stars_invariants")
    stars_asterisms_grp = fout.create_group("stars_asterisms")

    # Calculate the number of healpix polygons
    npix = hp.nside2npix(nside)
    fp_radecs = [] # Pointing of each healpix area

    # Process each healpix area
    for seq in range(npix):

        desc = f'Generating starcatalog sky area index {seq + 1} of {hp.nside2npix(nside)} for level k{k}'
        print(desc,end='\r')

        # Calculate the center of healpix polygon
        fp_radec = hp.pix2ang(nside,seq,lonlat=True)
        fp_radecs.append(fp_radec)

        # Perform a cone search of stars on the star catalogs.
        stars = self.search_cone(fp_radec,search_radius,max_num)
        # Calculates the pixel coordinates for each star
        stars.pixel_xy(pixel_width)
        # Calculates invariant features for the set of stars, used in pattern recognition.
        stars.invariantfeatures()

        stars_xy_grp.create_dataset(str(seq), data=stars.xy)
        stars_invariants_grp.create_dataset(str(seq), data=stars.invariants)
        stars_asterisms_grp.create_dataset(str(seq), data=stars.asterisms)  

    print('\nFinished')
        
    fout.create_dataset("fp_radecs", data=np.array(fp_radecs))
    fout.close()
    return outh5,ratio      

def _healpix_level(self, fov):
    """
    Determine the healpix level (k) and search radius based on the FOV.

    Usage:
        >>> _determine_healpix_level_and_radius(35)
    Inputs:
        fov -> [float] Field of View (FOV) of the camera in degrees.
    Outputs:
        k -> [int] Healpix level
        nside -> [int] 2 to the kth power, used to compute the number of healpix polygons
        search_radius -> [float] Derived search radius based on the FOV
    """
    if fov >= 29.3 and fov < 58.6: # solidangle_ratio(58.6,37.2) = 1, solidangle_ratio(28.4,37.2) = 1/5
        k,nside = 1,2
        search_radius = 37.2
    elif fov >= 14.7 and fov < 29.3:
        k,nside = 2,4
        search_radius = 17.0
    elif fov >= 7.33 and fov < 14.7:
        k,nside = 3,8
        search_radius = 8.3
    elif fov >= 3.66 and fov < 7.33:
        k,nside = 4,16
        search_radius = 4.1  
    elif fov >= 1.83 and fov < 3.66:
        k,nside = 5,32
        search_radius = 2.1
    else:
        raise ValueError("FOV should be within [1.83, 58.6] degrees.")
    return k, nside, search_radius
              
class Stars(object):
    """
    A class to represent a collection of stars.

    Attributes:
        - sc_name (str): The name of the star catalog from which the data is derived.
        - center (tuple of float): The center of the region of interest in the sky, expressed in [RA, DEC] coordinates (degrees).
        - df (Pandas DataFrame): A pandas DataFrame containing detailed information about each star.
        - num_stars (int): Total number of stars.
        - xy (2D array-like): Pixel coordinates of stars.
        - radec (2D array-like): Celestial coordinates (Right Ascension and Declination) of each star.
    Methods:
        - pixel_xy: Calculates the pixel coordinates of stars for a given pixel width and rotating orientation.
        - invariantfeatures: Generates invariant features for star pattern recognition.
    """   

    def __init__(self,info):
        """
        Initializes a new instance of the Stars class.

        Inputs:
            info -> [dict] A dictionary containing key-value pairs for the instance.
        """ 

        for key in info.keys():
            setattr(self, key, info[key])

    def __repr__(self):

        """
        String representation of the Stars instance, typically used for debugging and logging.

        Returns:
            A formatted string that provides a summary of the Stars instance.
        """
    
        return "<Stars object: CATALOG = '{:s}' CENTER = '{} in [RA,DEC]' SEARCH_AREA = '{} deg' STARS_NUM = {:d} MCP = {:d} MODE = '{:s}'>".format(
            self.sc_name,self.center,self.search_area,self.num_stars,self.max_control_points,self.mode)

    def pixel_xy(self,pixel_width,theta=0):
        """
        Calculates the pixel coordinates for each star in the collection.

        Usage:
            >>> from starcatalogquery import StarCatalogSimplified
            >>> # load the simplified star catalog HYG v3.7
            >>> dir_from_simplified = '/Volumes/TOSHIBA/starcatalogs/simplified/hygv3.7/res5/mag9.0/epoch2022.0/'
            >>> hygv37_simplified = StarCatalogReduced.load(hygv3.7',5,9,2022,dir_from_simplified)
            >>> hygv37_simplified_stars = hygv37_simplified.search_cone([20,30],10)
            >>> stars = hygv37_simplified_stars.pixel_xy(0.01)
        Inputs:
            pixel_width -> [float] The width of a pixel in degrees, used to convert celestial coordinates to pixel coordinates.
            theta -> [float,optional,default=0] Rotation angle (in radians) to align the WCS frame(equivalent to ENU) with the image reference frame.
        Outputs:
            The Stars instance with updated pixel coordinates for each star.
        """

        # Convert celestial coordinates to pixel coordinates of stars
        x,y,wcs = xy_catalog(self.center,self.radec,pixel_width,theta)

        df = self.df
        df['pixelx'],df['pixely'] = x,y
        self.xy = np.stack([x,y]).T
        self.wcs = wcs
        return self

    def invariantfeatures(self):
        """
        Generates geometrically invariant features based on the spatial configuration of stars, aiding in star pattern recognition.
        Steps:
            1. Derive unique geometric invariants (ratios L2/L1 and L1/L0) for every possible triangle formed by groups of three stars, where L2, L1, and L0 are the sides of each triangle sorted in descending order.
            2. Construct a KDTree structure using these unique invariants, facilitating efficient spatial queries and pattern matching.
            3. Associate each set of invariants with indices of the stars forming the corresponding triangle, preserving the relationship between the geometric features and the individual stars.
        Notes:
            This method should be called after pixel coordinates are calculated.
            It generates features that are invariant to rotation, scaling, and translation.
        """
        if not hasattr(self,'xy'): 
            raise Exception("The pixel coordinates of stars should be caiculated first by `.pixel_xy(pixel_width`)")
        # Derive unique geometric invariants    
        inv_uniq, triang_vrtx_uniq = unique_triangles(self.xy)
        # Construct a KDTree structure using the unique invariants
        inv_uniq_tree = KDTree(inv_uniq)
        self.invariants,self.asterisms,self.kdtree = inv_uniq,triang_vrtx_uniq,inv_uniq_tree
        return self
