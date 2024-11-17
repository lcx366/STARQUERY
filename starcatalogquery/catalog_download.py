import os
import numpy as np
import healpy as hp
import pandas as pd
from gzip import GzipFile

from .utils.try_download import wget_download
from .utils.starcatalog_statistic import tiles_statistic

NSIDE = 32  # Default HEALPix partition level, corresponding to K=5 (2^K=NSIDE)

def stsci_download(sc_name, mag_range, dir_to, url_file):
    """
    Downloads star catalog tile files, each representing a region of the sky.
    This function supports downloading files for different catalogs such as GAIA DR3, GSC30, UCAC5, USNOB, or 2MASS.

    Usage:
        >>> sc_name = 'gaiadr3'
        >>> dir_to = 'starcatalogs/raw/gaiadr3/'
        >>> url_file = 'starcatalogs/url/gaiadr3.txt'
        >>> stsci_download(sc_name, dir_to, url_file)
    Inputs:
        sc_name -> [str] Star catalog to download (e.g., 'gaiadr3', 'gsc30', 'ucac5', 'usnob', '2mass').
        mag_range -> [tuple] Range of magnitude, such as [3,17].
        dir_to -> [str] Directory for storing the star catalog files.
        url_file -> [str] URL file for downloading star catalog.
    Outputs:
        star catalog tile files
    """
    mag_min,mag_max = mag_range

    # Create the directory for tiles if it doesn't exist
    os.makedirs(dir_to, exist_ok=True)

    # Calculate the number of pixels for the given NSIDE
    npix = hp.nside2npix(NSIDE)

    print(f'Generating URL list file with {npix} entries.')

    url_prefix = 'https://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?'

    with open(url_file, 'w') as f:
        for k in range(npix):
            corners_vec = hp.boundaries(NSIDE, k)  # Get vector boundaries of each HEALPix pixel
            corners_lonlat = hp.vec2ang(np.transpose(corners_vec),lonlat=True)  # Convert vectors to longitudes and latitudes
            polygon_coords = ','.join(f"{lon} {lat}" for lon, lat in zip(*corners_lonlat))
            url = f'{url_prefix}STCS=POLYGON {polygon_coords}&format=csv&catalog={sc_name}&magrange={mag_min},{mag_max}'
            file_path = os.path.join(dir_to, f"{sc_name}-{k}.csv")
            f.write(f"'{file_path}' '{url}'\n")  # Write URL to file

    print(f'Downloading the catalog tile files for {sc_name} ...')
    os.system(f"cat '{url_file}' | xargs -P 16 -n 2 wget -c -O")
    print('Download complete.')


def hyg_download(sc_name, dir_to):
    """
    Downloads and converts HYG and AT-HYG databases into a tile-based format for easier handling.

    Inputs:
        sc_name -> [str] Name of the star catalog to download (e.g., 'hyg37' or 'at-hyg24').
        dir_to -> [str] Directory for storing the star catalog files.
    Outputs:
        dir_size -> [str] Total size of the star catalog.
        file_num -> [int] Number of the tile files.
        validity -> [bool] Validity of the star catalog.
    """
    # Create the directory for tiles if it doesn't exist
    os.makedirs(dir_to, exist_ok=True)

    catalog_info = {
        'hyg37': ('https://astronexus.com/downloads/catalogs/hygdata_v37.csv.gz', 'hygdata_v37.csv'),
        'at-hyg24': ('https://www.astronexus.com/downloads/catalogs/athyg_v24.csv.gz', 'athyg_v24.csv')
    }

    url, raw_file = catalog_info[sc_name]
    raw_dir_to = os.path.expanduser('~/src/sc-data/hyg/')
    raw_file_path = os.path.join(raw_dir_to, raw_file)
    raw_file_gz = os.path.join(raw_dir_to, os.path.basename(url))

    os.makedirs(raw_dir_to, exist_ok=True)  # Ensure the directory exists

    if not os.path.exists(raw_file_path):
        desc = f"Downloading {sc_name} from {url}"
        wget_out = wget_download(url, raw_file_gz, desc)
        with GzipFile(raw_file_gz) as gz:
            with open(raw_file_path, "wb") as f:
                f.write(gz.read())
        os.remove(raw_file_gz)
    else:
        print(f"Catalog file {raw_file} already present in {raw_dir_to}")

    # Dividing the downloaded star catalog into tiles
    print('Dividing the star catalog into tiles ...')
    df = pd.read_csv(raw_file_path, skiprows=[1], dtype=str)  # Skip the sun
    ra, dec = df['ra'].astype(float) * 15, df['dec'].astype(float)  # Convert RA to degrees and get DEC
    pixels = hp.ang2pix(NSIDE, ra, dec, lonlat=True)  # Calculate HEALPix pixel number for each star
    df['pixel'] = pixels  # Add pixel numbers as a new column in the DataFrame

    # Group stars by HEALPix pixel and save to separate CSV files
    for pixel_number, group in df.groupby('pixel'):
        tile_file = dir_to + '{:s}-{:d}.csv'.format(sc_name, pixel_number)
        with open(tile_file, 'w') as fn:
            fn.write('#Objects found: {:d}\n'.format(len(group)))
        group.to_csv(tile_file, mode='a', index=False, header=True)

    print('Finished processing the catalog.')

    # Make statistics for the tiles
    file_num, dir_size, validity = tiles_statistic(dir_to)

    return dir_size, file_num, validity
