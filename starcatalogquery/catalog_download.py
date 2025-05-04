import os
import numpy as np
import healpy as hp
import pandas as pd

from tqdm import tqdm
from gzip import GzipFile

from .utils.try_download import wget_download
from .utils.starcatalog_statistic import tiles_statistic

NSIDE = 32  # Default HEALPix partition level, corresponding to K=5 (2^K=NSIDE)

def stsci_download(sc_name, mag_range, url_file, dir_to):
    """
    Downloads tile-based star catalog files from the STScI VO service.

    For each HEALPix pixel (at the given NSIDE), this function:
      - Constructs a polygonal STCS query for the sky region
      - Writes the corresponding API request URL to a file
      - Downloads each tile in parallel using wget

    Supported catalogs: GAIA DR3, GSC30, UCAC5, USNOB, 2MASS (and others via STScI VO)

    Usage:
        >>> sc_name = 'gaiadr3'
        >>> mag_range = (3, 17)
        >>> url_file = 'starcatalogs/url/gaiadr3.txt'
        >>> dir_to = 'starcatalogs/raw/gaiadr3'
        >>> stsci_download(sc_name, mag_range, url_file, dir_to)

    Inputs:
        sc_name -> [str] Name of the star catalog (e.g., 'gaiadr3', 'gsc30', '2mass', etc.)
        mag_range -> [tuple] Magnitude range filter as (min_mag, max_mag)
        url_file -> [str] URL file listing the API request URL
        dir_to -> [str] Directory to save tile files.

    Outputs:
        - URL list file: <url_file>
        - Tile file per HEALPix pixel (e.g., gaiadr3-1234.csv)
    """
    mag_min,mag_max = mag_range
    os.makedirs(dir_to, exist_ok=True)

    # Calculate the number of pixels for the given NSIDE
    npix = hp.nside2npix(NSIDE)

    print(f'Generating URL list file with {npix} entries.')

    url_prefix = 'https://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?'

    with open(url_file, 'w') as f:
        for k in range(npix):
            corners_vec = hp.boundaries(NSIDE, k, nest=True)  # Get vector boundaries of each HEALPix pixel
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
    Downloads and processes HYG or AT-HYG star catalogs into HEALPix tile-based CSV files.

    This function supports both single-file and multi-part catalogs. It downloads the catalog data,
    merges parts if necessary, and splits the full catalog into smaller, tile-based CSV files.

    Usage:
        >>> hyg_download("hyg41",'starcatalogs/raw')

    Inputs:
        sc_name -> [str] Catalog name. Supported values: 'hyg41' or 'at-hyg24'.
        dir_to -> [str] Directory to save tile files.

    Returns:
        dir_size -> [str] Total size of all generated tile files (e.g., "24.51 MB").
        file_num -> [int] Number of tile files generated.
        validity -> [bool] True if tile files were successfully generated and non-empty.

    Output:
        - All tile files are written to the given `dir_to` directory.
        - Each tile corresponds to one HEALPix pixel number (NSIDE = 32 by default).
        - File naming convention: `<sc_name>-<pixel_id>.csv`, for example:
              at-hyg32-231.csv
              hyg41-1765.csv
        - Each tile CSV begins with a header line:
              #Objects found: <number>
        - These files can be used for spatially-indexed catalog fast lookups.
    """
    os.makedirs(dir_to, exist_ok=True)

    # Catalog metadata: URLs and final merged file name
    catalog_info = {
        'hyg41': {
            'urls': ['https://astronexus.com/downloads/catalogs/hygdata_v41.csv.gz'],
            'output': 'hygdata_v41.csv'
        },
        'at-hyg32': {
            'urls': [
                'https://codeberg.org/astronexus/athyg/media/branch/main/data/athyg_v32-1.csv.gz',
                'https://codeberg.org/astronexus/athyg/media/branch/main/data/athyg_v32-2.csv.gz'
            ],
            'output': 'athyg_v32.csv'
        }
    }

    # Directory where raw catalog files will be stored
    dir_base = os.path.expanduser(f'~/src/sc-data')
    raw_dir_to = os.path.join(dir_base,'hyg')
    os.makedirs(raw_dir_to, exist_ok=True)

    urls = catalog_info[sc_name]['urls']
    raw_file = catalog_info[sc_name]['output']
    raw_file_path = os.path.join(raw_dir_to, raw_file)

    # If catalog file does not exist, download it
    if not os.path.exists(raw_file_path):
        if sc_name == 'hyg41':
            # Single file decompression
            url = urls[0]
            gz_path = os.path.join(raw_dir_to, os.path.basename(url))
            wget_download(url, gz_path, f"Downloading {sc_name} ...")
            with GzipFile(gz_path) as gz:
                with open(raw_file_path, "wb") as f_out:
                    f_out.write(gz.read())
            os.remove(gz_path)
        elif sc_name == 'at-hyg32':
            # Multi file concatenation, retaining the header of the first file
            with open(raw_file_path, "wb") as fout:
                for i, url in enumerate(urls):
                    gz_path = os.path.join(raw_dir_to, os.path.basename(url))
                    wget_download(url, gz_path, f"Downloading part {i + 1}/{len(urls)} ...")
                    with GzipFile(gz_path) as gz:
                        lines = gz.readlines()
                        fout.writelines(lines)
                    os.remove(gz_path)
            print(f"Merged parts saved to: {raw_file_path}")
    else:
        print(f"Catalog file {raw_file} already exists at {raw_dir_to}")
    # Load catalog and convert coordinates
    df = pd.read_csv(raw_file_path, skiprows=[1], dtype=str,na_values = [' ', ''], on_bad_lines = 'skip')  # Skip Sun entry
    ra = df['ra'].astype(float) * 15  # RA in degrees
    dec = df['dec'].astype(float)
    pixels = hp.ang2pix(NSIDE, ra, dec, nest=True, lonlat=True)
    df['pixel'] = pixels

    # Write one CSV per HEALPix pixel
    desc = 'Dividing the star catalog into HEALPix tiles'
    for pixel_number, group in tqdm(df.groupby('pixel'), desc=desc, unit="tile"):
        tile_file = os.path.join(dir_to, f"{sc_name}-{pixel_number}.csv")
        with open(tile_file, 'w') as fn:
            fn.write(f'#Objects found: {len(group)}\n')

        group.drop(columns='pixel', inplace=True)
        group['mag'] = pd.to_numeric(group['mag'], errors='coerce')
        group = group.sort_values(by='mag')
        group.to_csv(tile_file, mode='a', index=False, header=True)

    # Return directory statistics
    file_num, dir_size, total_lines, validity = tiles_statistic(dir_to)
    stars_num = total_lines - file_num*2
    return dir_size, file_num, stars_num, validity
