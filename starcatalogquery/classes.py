# Load external packages
import os,warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from glob import glob
from astropy.time import Time
from tqdm import tqdm
from natsort import natsorted

# Load internal function library
from .catalog_index import estimate_mag_limit,estimate_level,find_healpix_level,build_catalog_indices,h5_hashes,read_h5_hashes
from .catalog_download import stsci_download,hyg_download
from .catalog_check import stsci_check,update_star_catalog,convert_csv_to_parquet
from .utils.starcatalog_statistic import tiles_statistic,starcatalog_info
from .utils.df2info import df2info
from .catalog_query import search_box_simplified,search_cone_simplified
from .tiles_draw import search_draw
from .wcs import xy_catalog
from .invariantfeatures import calculate_invariantfeatures
from .astrometry_corrections import parallax2dist

TILE_SIZE = 1.83  # Pixel size (in degrees) for HEALPix hierarchy level 'K5'

def format_columns(df, format_dict = {'ra': 8, 'dec': 8, 'pm_ra': 3, 'pm_dec': 3, 'dist': 8, 'mag': 3, 'epoch': 6}):
    """
    Formats the specified columns of a DataFrame according to the provided format dictionary.

    Inputs:
        df -> [pandas.DataFrame] The DataFrame to be formatted.
        format_dict -> [dict] A dictionary where the keys are column names and the values are the number of decimal places to round to.
    Returns:
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

    def get(sc_name,mag_range=None,dir_base=None):
        """
        Downloads 'raw' star catalog files from remote servers or loads them from a local directory.

        Usage:
            >>> from starcatalogquery import StarCatalog
            >>> sc_raw = StarCatalog.get('at-hyg32')
        Inputs:
            sc_name -> [str] Name of the star catalog.
            mag_range -> [tuple,optional,default=None] Range of magnitudes.
            Available star catalog include 'hyg41','at-hyg32','gaiadr3','gsc30','ucac5','usnob','2mass', etc.
            For more star catalogs, please refer to the Space Telescope Science Institute:
            https://outerspace.stsci.edu/display/GC/WebServices+for+Catalog+Access
            dir_base -> [str, optional, default=None] Base directory to save files.
            If None, a built-in directory is used by default.
        Returns:
            sc_raw -> [StarCatalogRaw object] An instance of the StarCatalogRaw class.
        """
        # List of supported star catalogs
        valid_catalogs = ['hyg41','at-hyg32','gaiadr3','gsc30','ucac5','usnob','2mass']
        if sc_name not in valid_catalogs:
            raise ValueError(f'Star catalog {sc_name} is not supported. Valid catalogs are: {valid_catalogs}')

        # Gather star catalog information
        stars_num, mag, description = starcatalog_info(sc_name)

        if dir_base is None:
            dir_base = os.path.expanduser(f'~/src/sc-data')
        dir_to = os.path.join(dir_base, f'starcatalogs/raw/{sc_name}/')

        # Handle specific catalogs that require different download methods
        if sc_name in valid_catalogs[:2]:
            # For HYG and AT-HYG star catalogs
            dir_size, file_num, stars_num, validity = hyg_download(sc_name,dir_to)
        else:
            # For GAIA DER3, GSC 30, UCAC5, USNOB, 2MASS star catalogs
            dir_url = os.path.join(dir_base, 'starcatalogs/url')
            url_file = os.path.join(dir_url, f'{sc_name}.txt')
            os.makedirs(dir_url, exist_ok=True)

            # Limit the star magnitude range, otherwise remote servers may cause data overflow issues, leading to download failures
            if mag_range is None:
                mag_range = (3, 17)
            else:
                mag_range = (max(mag_range[0], 3), mag_range[1])

            if os.path.exists(dir_to):
                # Check existing star catalog for validity
                dir_size,file_num,validity = stsci_check(sc_name,dir_to,url_file)
            else:
                # Download and check star catalog if not already present
                stsci_download(sc_name, mag_range, url_file, dir_to)
                dir_size,file_num,validity = stsci_check(sc_name,dir_to,url_file)

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
            'description': description
        }

        return StarCatalogRaw(info)  

    def load(dir_from):
        """
        Load star catalogs at various levels of hierarchy from local storage.

        Usage:
            >>> from starcatalogquery import StarCatalog
            >>> # load the raw star catalog AT-HYG v3.2
            >>> dir_from_raw = 'starcatalogs/raw/at-hyg32/'
            >>> at_hyg32_raw = StarCatalog.load(dir_from_raw)
            >>>
            >>> # load the reduced star catalog AT-HYG v3.2
            >>> dir_from_reduced = 'starcatalogs/reduced/at-hyg32/'
            >>> at_hyg32_reduced = StarCatalog.load(dir_from_reduced)
            >>>
            >>> # load the simplified star catalog AT-HYG v3.2
            >>> dir_from_simplified = 'starcatalogs/simplified/at-hyg32/mag12.0/epoch2019.5/'
            >>> at_hyg32_simplified = StarCatalog.load(dir_from_simplified)
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
        parse_dir = os.path.normpath(dir_from).split(os.sep)

        if 'epoch' in parse_dir[-1]:
            _mode, sc_name, level_mag_threshold, epoch = parse_dir[-4:]
            level,mag_threshold = level_mag_threshold.split('-')
        else:
            _mode, sc_name = parse_dir[-2:]

        # load the star catalog
        if _mode == 'raw':
            starcatalog = StarCatalogRaw._load(sc_name,dir_from)
        elif _mode == 'reduced':
            starcatalog = StarCatalogReduced._load(sc_name,dir_from)
        elif _mode == 'simplified':
            mag_threshold = float(mag_threshold[3:])
            epoch = float(epoch[5:])
            lvl = int(level[3:])
            starcatalog = StarCatalogSimplified._load(sc_name,mag_threshold,lvl,epoch,dir_from)
            starcatalog.build_indices()
        else:
            raise Exception(f"Invalid star catalog level: {_mode}.")

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
        - sc_name (str): Name of the star catalog, such as 'at-hyg32'.
        - _mode (str): Level of the star catalog, set to 'raw' for this class.
        - stars_num (str): Total number of the stars included in the catalog.
        - mag (str): Range of stars magnitudes in the catalog.
        - description (str): A brief description or summary of the catalog.

    Methods:
        - _load: Loads the raw star catalog from a specified directory.
        - reduce: Reduces the raw star catalog to contain only essential information.
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
    
        return "<StarCatalogRaw object: CATALOG_NAME = '{:s}' CATALOG_SIZE = '{:s}' TILES_NUM = {:d} TILE_SIZE = '{:s}' STARS_NUM = '{:}' MAG = '{:s}'>".format(
            self.sc_name,self.sc_size,self.tiles_num,self.tile_size,self.stars_num,self.mag)

    @classmethod
    def _load(cls, sc_name, dir_from):
        """
        Load the raw star catalog files from a specified directory.

        Usage:
            >>> from starcatalogquery import StarCatalogRaw
            >>> # load the raw star catalog AT-HYG v3.2
            >>> dir_from_raw = 'starcatalogs/raw/at-hyg32/'
            >>> at_hyg32_raw = StarCatalogRaw.load('at-hyg32',dir_from_raw)
        Inputs:
            sc_name -> [str] Name of the star catalog (e.g., 'at-hyg32').
            dir_from -> [str] Directory for the star catalog files.
        Outputs:
            at_hyg32_raw -> [StarCatalogRaw object] Instance of the class StarCatalogRaw.
        """
        # Calculate total file size and total number of the tile files
        stars_num, mag, description = starcatalog_info(sc_name)
        file_num, dir_size, total_lines, validity = tiles_statistic(dir_from)
        stars_num = total_lines - file_num*2

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
            >>> at_hyg32_reduced = at_hyg32_raw.reduce()
        Outputs:
            at_hyg32_reduced -> [StarCatalogReduced object] The reduced star catalog.
        """
        # Copying existing information
        info = self.__dict__.copy()

        # Setting up the directory for the reduced catalog
        tiles_dir = self.tiles_dir
        dir_reduced = tiles_dir.replace('raw','reduced')
        os.makedirs(dir_reduced, exist_ok=True)

        file_list = glob(f'{tiles_dir}*')
        sc_name = self.sc_name
        desc = "Reducing tile files"
        if sc_name in ['hyg41','at-hyg32']:
            if sc_name == 'hyg41':
                original_columns = ['ra','dec','pmra','pmdec','dist','mag']
            elif sc_name == 'at-hyg32':
                original_columns = ['ra','dec','pm_ra','pm_dec','dist','mag']
            rename_dict = {'pmra':'pm_ra', 'pmdec':'pm_dec'}    
            # Processing each file in the raw star catalog directory
            for tile_file in tqdm(file_list, desc=desc, unit="file"):
                df = pd.read_csv(tile_file,skiprows=1,dtype=str,na_values=[' ', ''],on_bad_lines='skip',usecols=original_columns)
                # Applying specific transformations based on the catalog name
                # units: ra->hourangle, dec->deg, pmra->mas/a, pmdec->mas/a, dist->parsecs, epoch->2000.0
                df_reduced = df.apply(pd.to_numeric,errors='coerce')
                df_reduced = df_reduced.reindex(columns=original_columns)
                df_reduced['ra'] = df_reduced['ra']*15 # Convert hourangle to deg
                df_reduced['dist'] = df_reduced['dist'] / 1e3  # Convert pc to kpc
                df_reduced['epoch'] = 2000.0
                df_reduced.rename(columns=rename_dict, inplace=True)
                df_reduced = format_columns(df_reduced)
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)
        elif sc_name == 'gaiadr3':
            for tile_file in tqdm(file_list, desc=desc, unit="file"):
                # Define the columns to be read and the renamed mapping after reading
                original_columns = ['ra', 'dec', 'pmra', 'pmdec', 'parallax', 'mag', 'epoch']
                rename_dict = {'pmra': 'pm_ra', 'pmdec': 'pm_dec'}
                # Read star catalog tile files, select only the required columns
                df = pd.read_csv(tile_file,skiprows=1,dtype=str,na_values=[' ', ''],on_bad_lines='skip',usecols=original_columns)
                # units: ra->deg, dec->deg, pmra->mas/a, pmdec->mas/a, epoch->2016
                df_reduced = df.apply(pd.to_numeric,errors='coerce')
                df_reduced = df_reduced.reindex(columns=original_columns)
                # Calculate distance from parallax
                df_reduced['dist'] = parallax2dist(df_reduced['parallax'])
                parallax_index = df_reduced.columns.get_loc('parallax')
                # Remove the 'parallax' column from the DataFrame
                df_reduced.drop(columns=['parallax'], inplace=True)
                # Insert the 'dist' column in the position of the original 'parallax' column
                df_reduced.insert(parallax_index, 'dist', df_reduced.pop('dist'))
                # Rename columns
                df_reduced.rename(columns=rename_dict, inplace=True)
                df_reduced = format_columns(df_reduced)
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)
        elif sc_name == 'gsc30':
            for tile_file in tqdm(file_list, desc=desc, unit="file"):
                # Define the columns to be read and the renamed mapping after reading
                original_columns = ['ra','dec','rapm','decpm','parallax','mag','epoch']
                rename_dict = {'rapm':'pm_ra', 'decpm':'pm_dec'}
                df = pd.read_csv(tile_file,skiprows=1,dtype=str,na_values=[' ', ''],on_bad_lines='skip',usecols=original_columns)
                # units: ra->deg, dec->deg, rapm->mas/a, decpm->mas/a, epoch->2012
                df_reduced = df.apply(pd.to_numeric,errors='coerce')
                df_reduced = df_reduced.reindex(columns=original_columns)
                # Calculate distance from parallax
                df_reduced['dist'] = parallax2dist(df_reduced['parallax'])
                # Get the position of the 'parallax' column
                parallax_index = df_reduced.columns.get_loc('parallax')
                # Remove the 'parallax' column from the DataFrame
                df_reduced.drop(columns=['parallax'], inplace=True)
                # Insert the 'dist' column in the position of the original 'parallax' column
                df_reduced.insert(parallax_index, 'dist', df_reduced.pop('dist'))
                # Rename columns
                df_reduced.rename(columns=rename_dict, inplace=True)
                df_reduced = format_columns(df_reduced)
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)
        elif sc_name == 'ucac5':
            for tile_file in tqdm(file_list, desc=desc, unit="file"):
                # Define the columns to be read and the renamed mapping after reading
                original_columns = ['ra','dec','pmur','pmud','mag','epu']
                rename_dict = {'pmur':'pm_ra', 'pmud':'pm_dec','epu':'epoch'}
                df = pd.read_csv(tile_file,skiprows=1,dtype=str,na_values=[' ', ''],on_bad_lines='skip',usecols=original_columns)
                # units: ra->deg, dec->deg, pmur->mas/a, pmud->mas/a, epu->1998.754
                df_reduced = df.apply(pd.to_numeric,errors='coerce')
                df_reduced = df_reduced.reindex(columns=original_columns)
                # Rename columns
                df_reduced.rename(columns=rename_dict, inplace=True)
                df_reduced = format_columns(df_reduced)
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)
        elif sc_name == 'usnob':
            for tile_file in tqdm(file_list, desc=desc, unit="file"):
                # Define the columns to be read and the renamed mapping after reading
                original_columns = ['ra','dec','pmRA','pmDEC','mag','Epoch']
                rename_dict = {'pmRA':'pm_ra', 'pmDEC':'pm_dec','Epoch':'epoch'}
                df = pd.read_csv(tile_file,skiprows=1,dtype=str,na_values=[' ', ''],on_bad_lines='skip',usecols=original_columns)
                # units: ra->deg, dec->deg, pmRA->mas/a, pmDEC->mas/a, Epoch->1950,mag->unknown
                df_reduced = df.apply(pd.to_numeric,errors='coerce')
                df_reduced = df_reduced.reindex(columns=original_columns)
                # Rename columns
                df_reduced.rename(columns=rename_dict, inplace=True)
                df_reduced = format_columns(df_reduced)
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)  
        elif sc_name == '2mass':
            for tile_file in tqdm(file_list, desc=desc, unit="file"):
                # Define the columns to be read and the renamed mapping after reading
                original_columns = ['ra','dec','mag','jdate']
                df = pd.read_csv(tile_file,skiprows=1,dtype=str,na_values=[' ', ''],on_bad_lines='skip',usecols=original_columns)
                # units: ra->deg, dec->deg, jdate->2451063.6417
                df_reduced = df.apply(pd.to_numeric,errors='coerce')
                df_reduced = df_reduced.reindex(columns=original_columns)
                # Convert Julian date to epoch year
                df_reduced['epoch'] = Time(df_reduced['jdate'], format='jd').jyear
                df_reduced.pop('jdate')
                df_reduced = format_columns(df_reduced)
                df_reduced.to_csv(tile_file.replace('raw','reduced'),index=False)
                    
        file_num,dir_size,total_lines,validity = tiles_statistic(dir_reduced)
        stars_num = total_lines - file_num

        info['_mode'] = 'reduced'
        info['tiles_dir'] = dir_reduced
        info['sc_size'] = dir_size
        info['stars_num'] = stars_num

        return StarCatalogReduced(info)

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
        - sc_name (str): Name of the star catalog, such as 'at-hyg32'.
        - _mode (str): Level of the star catalog, set to 'reduced' for this class.
        - stars_num (str): Total number of the stars included in the catalog.
        - mag (str): Range of star magnitudes in the catalog.
        - description (str): A brief description or summary of the catalog.

    Methods:
        - _load: Loads the reduced star catalog from a specified directory.
        - simplify: Simplifies the reduced star catalog to a more basic version.
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
    
        return "<StarCatalogReduced object: CATALOG_NAME = '{:s}' CATALOG_SIZE = '{:s}' TILES_NUM = {:d} TILE_SIZE = '{:s}' STARS_NUM = '{:d}' MAG = '{:s}'>".format(
            self.sc_name,self.sc_size,self.tiles_num,self.tile_size,self.stars_num,self.mag)

    @classmethod
    def _load(cls,sc_name,dir_from):
        """
        Load the reduced star catalog files from a specified directory.

        Usage:
            >>> from starcatalogquery import StarCatalogReduced
            >>> # load the reduced star catalog HYG v3.7
            >>> dir_from_reduced = 'starcatalogs/reduced/at-hyg32/'
            >>> at_hyg32_reduced = StarCatalogReduced.load('at-hyg32',dir_from_reduced)
        Inputs:
            sc_name -> [str] Name of the star catalog (e.g., 'at-hyg32').
            dir_from -> [str] Directory for the star catalog files.
        Outputs:
            at_hyg32_reduced -> [StarCatalogReduced object] Instance of the class StarCatalogReduced.
        """ 

        # Calculate total file size and total number of the tile files
        stars_num, mag, description = starcatalog_info(sc_name)
        file_num, dir_size, total_lines, validity = tiles_statistic(dir_from)
        stars_num = total_lines - file_num

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

    def simplify(self, t_pm, fov=None, mag_threshold=None):
        """
        Simplify the reduced star catalog, applying magnitude truncation and proper motion correction，
        making it a basic version suitable for quick lookups.

        Usage:
            >>> at_hyg32_simplified = at_hyg32_reduced.simplify(9.0,2025.5)
        Inputs:
            t_pm -> [float] Epoch to which the stars are unified.
            fov -> [tuple,optional,default=None] FOV of camera, such as (7,9)
            mag_threshold -> [float,optional,default=None] Magnitude threshold to cutoff star catalog entries.
        Outputs:
            Instance of class StarCatalogSimplified
        """
        # Copying existing information
        info = self.__dict__.copy()

        if fov is not None:
            # Estimate the magnitude threshold
            level, nside, npix, pixel_size = find_healpix_level(fov_min=min(fov))
            lvl = int(level[1:])
            mag_threshold = estimate_mag_limit(lvl)
        else:
            if mag_threshold is None:
                raise AttributeError('No magnitude threshold specified.')
            mag_threshold = round(mag_threshold, 1)
            #  Compute the theoretical highest HEALPix order such that a star catalogue
            #  complete to magnitude *mag_threshold* still provides at least *30*
            #  stars **on average** in every pixel.
            lvl_max = estimate_level(mag_threshold)
            lvl = int(f'99{lvl_max}')  # Indicates that the healpix order has not been determined yet

        # Setting up the directory for the reduced catalog
        tiles_dir = self.tiles_dir
        dir_simplified = os.path.join(tiles_dir.replace('reduced', 'simplified'), f'lvl{lvl}-mag{mag_threshold:.1f}',
                                      f'epoch{t_pm:.1f}')
        os.makedirs(dir_simplified, exist_ok=True)

        # Construct file pattern
        sc_name = self.sc_name
        file_pattern = os.path.join(tiles_dir, f'{sc_name}-*.csv')
        file_list = natsorted(glob(file_pattern))

        desc = 'Simplifying tile files'
        # Processing each file in the reduced star catalog directory
        for tile_file in tqdm(file_list, desc=desc, unit="file"):
            df_simplified = pd.read_csv(tile_file).dropna()
            mag = df_simplified['mag']
            mag_flag = mag < mag_threshold
            df_simplified = df_simplified[mag_flag]

            # Correct proper motion
            dt = t_pm - df_simplified['epoch']

            if {'pm_ra', 'pm_dec'}.issubset(df_simplified.columns):
                df_simplified['ra'] += df_simplified['pm_ra'] / 3.6e6 * dt
                df_simplified['dec'] += df_simplified['pm_dec'] / 3.6e6 * dt
                df_simplified['pm_ra'] = df_simplified['pm_ra']
                df_simplified['pm_dec'] = df_simplified['pm_dec']
                df_simplified['epoch'] = t_pm
            else:
                warnings.warn(f'Proper motion data for catalog {sc_name} is not found.')

            if 'dist' in df_simplified.columns:
                dist_flag = df_simplified['dist'] > 0
                df_simplified = df_simplified[dist_flag]

            df_simplified = format_columns(df_simplified)
            csv_path = os.path.join(dir_simplified, os.path.basename(tile_file))
            parquet_path = csv_path.replace('.csv', '.parquet')

            # Write DataFrame to Parquet file
            df_simplified.to_csv(csv_path, index=False)

        update_star_catalog(dir_simplified,sc_name)
        convert_csv_to_parquet(dir_simplified)

        file_num, dir_size, total_lines, validity = tiles_statistic(dir_simplified, 'parquet')
        stars_num = total_lines

        info['_mode'] = 'simplified'
        info['tiles_dir'] = dir_simplified
        info['sc_size'] = dir_size
        info['mag_threshold'] = mag_threshold
        info['epoch'] = t_pm
        info['lvl'] = lvl
        info['stars_num'] = stars_num

        starcatalog_simplified = StarCatalogSimplified(info)
        starcatalog_simplified.build_indices()

        return starcatalog_simplified

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

    Methods:
        - _load: Loads the simplified star catalog from a specified directory.
        - search_box: Searches for stars within a specified rectangular area.
        - search_cone: Searches for stars within a specified conical area.
        - build_indices: Aggregates indices from all tile files and generates a comprehensive catalog index file.
        - check_lvl: Check consistency between the camera-derived HEALPix level and the star catalog level.
        - hashes: Generates a h5-formatted star catalog geometric invariants hashed file.
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

        return "<StarCatalogSimplified object: CATALOG_NAME = '{:s}' CATALOG_SIZE = '{:s}' TILES_NUM = {:d} TILE_SIZE = '{:s}' STARS_NUM = '{:d}' LVL = '{:d}' MAG_CUTOFF = {:.1f} EPOCH = {:.1f}>".format(
            self.sc_name,self.sc_size,self.tiles_num,self.tile_size,self.stars_num,self.lvl,self.mag_threshold,self.epoch)

    @classmethod
    def _load(cls,sc_name,mag_threshold,lvl,epoch,dir_from):
        """
        Load the simplified star catalog files from a specified directory.

        Usage:
            >>> from starcatalogquery import StarCatalogSimplified
            >>> # load the simplified star catalog AT-HYG v3.2
            >>> dir_from_simplified = 'starcatalogs/simplified/at_hyg32/mag12.0/epoch2019.5/'
            >>> at_hyg32_simplified = StarCatalogSimplified.load('at-hyg32',5,9,2022,dir_from_simplified)
        Inputs:
            sc_name -> [str] Name of the star catalog (e.g., 'at-hyg32').
            dir_from -> [str] Directory for the star catalog files.
        Outputs:
            at_hyg32_simplified -> [StarCatalogSimplified object] Instance of the class StarCatalogSimplified with loaded data.
        """  

        # Calculate total size and number of tile files
        stars_num,mag,description = starcatalog_info(sc_name)
        file_num, dir_size, total_lines, validity = tiles_statistic(dir_from,'parquet')
        stars_num = total_lines

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
            'lvl': lvl,
            'epoch': epoch,
        }

        return cls(info) 

    def build_indices(self):
        """
        Aggregate indices from all tile files for each star across different HEALPix hierarchy levels,
        generating a comprehensive csv-formatted catalog index file.

        Usage:
            >>> at_hyg32_simplified.build_indices()

        Outputs:
            - A sorted Parquet file located at:
              <project_root>/starcatalogs/indices/<sc_name>.parquet
            - The index file includes one row per star, with columns:
              K6, K7, K8, K9, K10, K5_SUB
        Notes:
            - * lvl ≤ 5  : no index is built (too coarse for practical use).
        """
        lvl = self.lvl
        tb_name = f'{self.sc_name}_lvl{lvl}-mag{self.mag_threshold}_epoch{self.epoch}'
        indices_dir = os.path.join(self.tiles_dir.split('starcatalogs')[0], 'starcatalogs/indices')
        indices_path = os.path.join(indices_dir, f'{tb_name}.parquet')

        #  Decide k_list and k_range
        if lvl > 99:  # encoded form, e.g. 998 ⇒ lvl_max = 8
            lvl = int(str(lvl)[2:])
            k_list = range(6, lvl + 1)  # index 6 … lvl_max
            k_range = (1, lvl)
        else:  # normal case
            k_list = range(lvl, lvl + 1)
            k_range = (k_list.start, k_list.stop - 1)

        #  Build index only when it makes sense and is not yet present
        if lvl > 5 and not os.path.exists(indices_path):
            os.makedirs(indices_dir, exist_ok=True)
            build_catalog_indices(self.tiles_dir,self.sc_name,indices_path,k_list)

        self._indices_path = indices_path
        self.tb_name = tb_name
        self.k_range = k_range

    def check_lvl(self,camera_params):
        """
        Check consistency between the camera-derived HEALPix level and the star catalog level.

        This function verifies whether the HEALPix level determined from camera parameters
        (such as field of view, pixel size, and resolution) matches the level information
        specified in the simplified star catalog (`sc_simplified`).

        - If the star catalog explicitly specifies a level, it must match the camera's expected level.
        - If the catalog does not specify a level (level > 99), then:
            - The catalog must support at least the required level (source level <= catalog max level).
            - Otherwise, an error is raised.
        - If the camera parameters are not provided, but the catalog specifies a level, raise an error.

        Inputs:
            camera_params -> [dict]
                Dictionary containing camera information, must include key 'fov' (field of view in degrees).
                Example: {'fov': [fov_width, fov_height], ...}
        """
        self._lvl_tmp = None
        lvl_catalog = self.lvl
        fov = camera_params.get('fov',0)

        if fov:
            # Calculate the expected HEALPix level based on the minimum FOV
            lvl_source,_nside,_npix, _pixel_size = find_healpix_level(fov_min=min(fov))
            lvl_source = int(lvl_source[1:])

            if lvl_catalog > 99:
                # Catalog does not specify a level, parse maximum level
                lvl_max = int(str(lvl_catalog)[2:])
                if lvl_source > lvl_max:
                    raise Exception(
                        f'Given FOV {fov} deg, the star catalog maximum level {lvl_max} is insufficient (required: {lvl_source}).')
                else:
                    self.case = '1-2'
                    warnings.warn('FOV is specified, but star catalog level is unspecified.')
            else:
                # Catalog specifies a definite level, must match exactly
                if lvl_source != lvl_catalog:
                    raise Exception(
                        f'Camera-derived HEALPix level ({lvl_source}) does not match catalog level ({lvl_catalog}).')
                self.case = '1-1'
            self._lvl_tmp = lvl_source
        else:
            if lvl_catalog < 99:
                # Camera FOV missing, but catalog level specified: inconsistent
                raise Exception('Camera FOV is not provided, but star catalog specifies a level.')
            self.case = '0-2'

    def search_box(self,radec_box,max_num=None,max_num_per_tile=None,lvl=None,astrometry_corrections={}):
        """
        Perform a rectangle search of stars on simplified star catalogs.

        Usage:
            >>> stars = at_hyg32_simplified.search_box([20,30,30,40])
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
        sc_name = self.sc_name
        indices_path = self._indices_path
        if lvl is None: lvl = self.lvl

        # Performs a rectangular search on the simplified star catalog
        df,level,nside,ids,pixel_size = search_box_simplified(radec_box, dir_sc, sc_name, indices_path, lvl, max_num_per_tile,astrometry_corrections)
        info = df2info(sc_name,center,df,max_num,level,nside,ids,pixel_size,search_area)
        return Stars(info)

    def search_cone(self,center,radius,max_num=None,max_num_per_tile=None,lvl=None,astrometry_corrections={}):
        """
        Perform a conical search of stars on simplified star catalogs.

        Usage:
            >>> stars = at_hyg32_simplified.search_cone([20,30],10)
        Inputs:
            center -> [list] Center of the cone in form of [ra_c, dec_c] in degrees.
            radius -> [float] Angular radius of the cone.
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
            stars -> [Stars object] An instance of class Stars with the search results.
        """
        search_area = {'cone':(center,radius)}

        dir_sc = self.tiles_dir
        sc_name = self.sc_name
        indices_path = self._indices_path
        if lvl is None: lvl = self.lvl

        # Performs a cone search on the simplified star catalog
        df,level,nside,ids,pixel_size = search_cone_simplified(center, radius, dir_sc, sc_name, indices_path, lvl, max_num_per_tile,astrometry_corrections)
        info = df2info(sc_name,center,df,max_num,level,nside,ids,pixel_size,search_area)
        return Stars(info)

    def hashes(self, mode_invariants='quads', n_stars=5, num_nearest_neighbors=10):
        """
        Loads a h5-formatted hashed file containing geometric invariants such as triangles or quads.
        If not exists, generate it.

        Inputs:
            mode_invariants -> [str] Type of invariant to load, must be one of ['triangles', 'quads'].
            n_stars -> [int,optional,default=9] Number of stars per tile for each level.
            num_nearest_neighbors -> [int,optional,default=15] Number of nearest neighbors to consider for each point.
        Returns:
            H5HashesData: Parsed hash data loaded from .h5 file.
        """
        # Validate invariant mode
        if mode_invariants not in ['triangles', 'quads']:
            raise ValueError(f"Unrecognized mode invariants type: {mode_invariants}")

        # Determine hashes directory based on tiles_dir
        hashes_dir = os.path.join(self.tiles_dir.split('starcatalogs')[0], 'starcatalogs/hashes')

        # Construct expected file pattern
        pattern = os.path.join(hashes_dir, f"{self.tb_name}_{mode_invariants}_*.h5")
        matched_files = glob(pattern)

        # If matching file exists, use it; otherwise generate it
        if matched_files:
            hash_file = matched_files[0]
            hashed_data = read_h5_hashes(hash_file) # Read from the hash file and wrap in H5HashesData
        else:
            k_min, k_max = self.k_range
            hash_file, hashed_data = h5_hashes(self._indices_path, self.tiles_dir, self.sc_name, self.tb_name, k_min, k_max, n_stars, num_nearest_neighbors, mode_invariants)
        return H5HashesData(self, hash_file, hashed_data, mode_invariants)

class H5HashesData:
    """
    A class to manage hashed data related to simplified star catalogs and mode invariants.

    Attributes:
        - sc_simplified -> [StarCatalogSimplified object] The simplified star catalog data.
        - hashed_data -> [dict] The hashed representation of the geometric invariants data.
        - mode_invariants -> [str] Mode of geometric invariants
    """

    def __init__(self, sc_simplified, hash_file, hashed_data, mode_invariants):
        """
        Initializes a new instance of the H5HashesData class.

        Inputs:
            - sc_simplified -> [StarCatalogSimplified object] The simplified star catalog data.
            - hashed_data -> [dict] The hashed representation of the geometric invariants data.
            - mode_invariants -> [str] Mode of geometric invariants
        """
        self.sc_simplified = sc_simplified
        self.hash_file = hash_file
        self.hashed_data = hashed_data
        self.mode_invariants = mode_invariants

    def __repr__(self):
        """
        String representation of the H5HashesData instance.

        Returns:
            A formatted string that provides a summary of the H5HashesData instance.
        """

        return f"<H5HashesData object: mode_invariants={self.mode_invariants}>"


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
        - tiles_draw: Visualize the scope of the search area and the coverage of the corresponding tiles.
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
        Calculates pixel coordinates of the stars.

        Usage:
            >>> stars = at_hyg32_simplified_stars.pixel_xy(0.01)
        Inputs:
            pixel_width -> [float] Pixel width in degrees, used to convert celestial coordinates to pixel coordinates.
            theta -> [float,optional,default=0] Rotation angle (in radians) to align the WCS frame(equivalent to ENU) with the image reference frame.
        """

        # Convert celestial coordinates to pixel coordinates of stars
        x,y,wcs = xy_catalog(self.center,self.radec,pixel_width,theta)

        df = self.df
        df['pixelx'],df['pixely'] = x,y
        self.xy = np.stack([x,y]).T
        self.wcs = wcs

    def invariantfeatures(self,num_nearest_neighbors,mode_invariants='triangles'):
        """
        Generates geometric invariant features based on the spatial configuration of stars, aiding in star pattern recognition.

        Inputs:
            num_nearest_neighbors -> [int] Number of nearest neighbors to consider for each point.
            mode_invariants -> [str, optional, default='triangles'] Mode of geometric invariants to use, e.g., 'triangles' or 'quads'.

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
        inv_uniq,vrtx_uniq,inv_uniq_tree = calculate_invariantfeatures(self.xy, num_nearest_neighbors, mode_invariants)
        self.invariants,self.asterisms,self.kdtree = inv_uniq,vrtx_uniq,inv_uniq_tree

    def tiles_draw(self):
        """
        Visualize the scope of the search area and the coverage of the corresponding tiles.

        Usage:
            >>> at_hyg32_simplified.tiles_draw()
        Outputs:
            An image illustrating the scope of the search area and the coverage of corresponding tiles.
        """
        search_draw(self.nside,self.tiles_ids,self.search_area,self.radec,self.level)

    def plot_scatter(self, width, height, output_path='stars.png'):
        """
        Plot a scatter diagram for the given stars, with the x and y axes scaled equally.

        Inputs:
            xy -> [numpy.ndarray] An array of shape (n, 2) containing the (x, y) coordinates of n stars.
            width -> [int] The width resolution of the output image.
            height -> [int] The height resolution of the output image.
            output_path -> [str,optional,default='stars.png'] Path to the output image.
        Outputs:
            Scatter diagram of stars
        """
        # Extract x and y coordinates
        x,y = self.xy.T

        # Filter points within the specified range
        x_min, x_max = -width / 2, width / 2
        y_min, y_max = -height / 2, height / 2
        mask = (x >= x_min) & (x <= x_max) & (y >= y_min) & (y <= y_max)
        x_filtered = x[mask]
        y_filtered = y[mask]

        # Plot the scatter diagram
        plt.figure(figsize=(6, 6), dpi=300)
        plt.scatter(x_filtered, y_filtered, s=5)
        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)
        plt.gca().set_aspect('equal', adjustable='box')  # Set the aspect of the axes to be equal
        plt.savefig(output_path,bbox_inches='tight')
        plt.close()
