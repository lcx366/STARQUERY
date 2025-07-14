import os,re
import pandas as pd
import healpy as hp

from glob import glob
from natsort import natsorted
from tqdm import tqdm

NUM_TILES = 12288  # Default number of HEALPix pixels for (K=5, NSIDE=32)

def find_smallest_csv_file(dir_from):
    """
    Identify the smallest CSV file in a given directory.

    Inputs:
        dir_from -> [str] Path to the directory containing star catalog tiles.
        Example: '/path/to/starcatalogs/raw/<catalog_name>/'
    Returns:
        smallest_file -> [str] Full path to the smallest CSV file in the directory.
    """
    csv_files = glob(os.path.join(dir_from, '*.csv'))
    smallest_file = min(csv_files, key=os.path.getsize, default=None)

    return smallest_file

def read_urls(file_path):
    """
    Extract the URL from each line in a text file.
    Each line should contain two quoted strings:
        1. The relative path to a catalog CSV file
        2. The associated download URL

    Usage:
        >>> file_path = 'starcatalogs/url/gaiadr3.txt'
        >>> urls = read_urls(file_path)
    Inputs:
        file_path -> [str] Path to the file containing lines with quoted CSV paths and URLs.
    Returns:
        urls -> [list of str] A list of extracted URLs.
    """
    urls = []
    with open(file_path, 'r') as file:
        for line in file:
            url = re.findall(r"'(.*?)'", line)[1]
            urls.append(url)

    return urls

def stsci_check(sc_name, dir_from, url_file):
    """
    Checks the validity of star catalog tile files downloaded from STScI, and re-fetches invalid files from a remote server.

    Usage:
        >>> dir_from = 'starcatalogs/raw/gaiadr3/'
        >>> dir_size, file_num, validity = stsci_check('gaiadr3', dir_from)
    Inputs:
        sc_name -> [str] Name of the star catalog.
        dir_from -> [str] Directory of the star catalog tile files. Expected format: '<your_path>/starcatalogs/raw/<sc_name>/'
        url_file -> [str] The path to the URL file.
    Returns:
        - dir_size -> [str] Total size of the star catalog tile files with appropriate units (KB, MB).
        - file_num -> [int] Total number of the star catalog tile files.
        - validity -> [bool] Indicates whether the star catalog is complete and safe to use.
    """
    urls = read_urls(url_file)

    # List all CSV files in the directory and sort them
    file_pattern = os.path.join(dir_from, f'{sc_name}-*.csv')
    files_list = natsorted(glob(file_pattern))
    file_num = len(files_list)

    # Initialize directory size and issue flag
    dir_size = 0
    issue_flag = False

    print(f'Checking and validating the star catalog tile files for {sc_name} ...')

    # Prepare a temporary file to log URLs of invalid files for re-download
    temp_url = f'temp_{sc_name}.txt'

    with open(temp_url, 'w') as outputf:
        # Iterate over each expected tile and check its validity
        for i in tqdm(NUM_TILES, desc="Checking .csv files", unit="file"):
            file_path = os.path.join(dir_from, f'{sc_name}-{i}.csv')

            if file_path in files_list:
                file_size = os.path.getsize(file_path)
                dir_size += file_size / 1024  # Convert bytes to KB

                # Output the total number of rows in the star catalog file
                file_lines = 0
                with open(file_path, 'rb') as f:
                    while chunk := f.read(chunk_size):
                        file_lines += chunk.count(b'\n')
                    # Subtract the file description line and header line
                    star_count = file_lines -2

                # Check the first line for object count consistency
                with open(file_path, 'r') as fn:
                    firstline = fn.readline().strip()

                expected_count = int(firstline.split(' : ')[1]) if 'Objects found' in firstline else 1
                if star_count != expected_count:
                    # Log invalid or missing files for re-download
                    issue_flag = True
                    outputf.write(f"'{file_path}' '{urls[i]}'\n")  # Write URL to file
                    os.remove(file_path)
            else:
                issue_flag = True
                outputf.write(f"'{file_path}' '{urls[i]}'\n")

    # Re-download invalid files if necessary
    if issue_flag:
        cmd = f"cat '{temp_url}' | xargs -P 16 -n 2 wget -c -O"
        os.system(cmd)
        # Extract the smallest star catalog tile file
        smallest_tile_file = find_smallest_csv_file(dir_from)
        # Verify the validity of the star catalog tile file
        with open(smallest_tile_file, 'r') as csvfile:
            first_line = csvfile.readline().strip()

        if "Objects found" not in first_line:
            print('The star catalog tile files are incomplete, which need to be confirmed manually. You may try downloading again.')
    else:
        print('All star catalog tile files are valid.')

    os.remove(temp_url)  # Clean up the temporary file

    # Format the directory size into a human-readable format
    dir_size_str = f'{dir_size:.2f} KB' if dir_size < 1024 else f'{dir_size / 1024:.2f} MB'
    validity = not issue_flag

    return dir_size_str, file_num, validity

def _get_pixel_id_from_bak(path: str) -> int:
    """
    Extracts the HEALPix pixel ID from a backup star catalog filename.

    Usage:
        Input: "/path/to/gaiadr3-123.csv.bak"
        Output: 123

    Inputs:
        path -> [str] Full path to a .csv.bak file.

    Returns:
        pixel_id -> [int] Pixel ID extracted from the filename.
    """
    pixel_id = os.path.basename(path).split('-')[-1].split(".")[0]
    return int(pixel_id)

def update_star_catalog(dir_sc: str, sc_name: str, nside: int = 32):
    """
    Verifies and corrects star placement errors in HEALPix-tiled star catalog files.

    This function scans all star catalog tiles (CSV files), identifies stars that were assigned
    to incorrect HEALPix pixels, and redistributes them into the correct tiles. It also re-sorts
    all tile files by magnitude from bright to faint.

    Usage:
        >>> update_star_catalog("sc-data/starcatalogs/raw/gaiadr3/", "gaiadr3", nside=32)

    Inputs:
        dir_sc -> [str] Directory containing the star catalog tiles (CSV files).
        sc_name -> [str] Filename of the star catalog (e.g., "gaiadr3").
        nside -> [int,optional,default=32] HEALPix resolution parameter.

    Output:
        - All original tile files are renamed to `.csv.bak` during processing.
        - Corrected tiles are written back to `.csv` with the same name format:
              <sc_name>-<pixel_id>.csv
        - Each file contains a line at the top: `#Objects found: <N>`
        - All tile files are sorted by 'mag' column from brightest to faintest.
        - If no star data is present, no output file will be created for that tile.
    """
    # Step 1: Rename all .csv files to .csv.bak for safe processing
    pattern = os.path.join(dir_sc, f"{sc_name}-*.csv")
    all_csv = glob(pattern)

    for fpath in all_csv:
        os.rename(fpath, fpath + ".bak")

    # Initialize header tracking for each pixel
    npix = 12 * (nside ** 2)
    header_written = [False] * npix

    # Process each .csv.bak file
    all_bak = glob(os.path.join(dir_sc, f"{sc_name}-*.bak"))
    for bak_file in tqdm(all_bak, desc="Redistributing stars into correct tile files", unit="file"):
        old_pix = _get_pixel_id_from_bak(bak_file)

        df = pd.read_csv(bak_file, comment='#')

        if df.empty:
            # Empty file => simply restore to .csv
            os.rename(bak_file, bak_file[:-4])
            header_written[old_pix] = True
            continue

        # Compute the correct pixel for each star
        correct_pix = hp.ang2pix(nside, df['ra'], df['dec'], nest=True, lonlat=True)
        df['pixel'] = correct_pix

        # Check if any stars are misassigned
        mismatch = (correct_pix != old_pix)
        if not mismatch.any():
            # No misplaced stars => rename back to .csv
            os.rename(bak_file, bak_file[:-4])
            header_written[old_pix] = True
        else:
            # Some stars need to be reassigned
            wrong_subdf = df.loc[mismatch].copy()
            # Group misplaced stars by target pixel
            for pix_id, subd in wrong_subdf.groupby('pixel'):
                out_file = os.path.join(dir_sc, f"{sc_name}-{pix_id}.csv")
                write_header = not header_written[pix_id]
                subd.drop(columns=['pixel'], inplace=True)
                subd.to_csv(out_file, mode='a', index=False, header=write_header)
                if write_header:
                    header_written[pix_id] = True

            # Write correctly placed stars back to their original tile
            correct_subdf = df.loc[~mismatch].drop(columns=['pixel'])
            out_file = os.path.join(dir_sc, f"{sc_name}-{old_pix}.csv")
            write_header = not header_written[old_pix]
            correct_subdf.to_csv(out_file, mode='a', index=False, header=write_header)
            if write_header:
                header_written[old_pix] = True

            # Remove the .bak file after correction
            os.remove(bak_file)

    # Step 2: Re-sort each output CSV by magnitude (brightest to faintest)
    for csv_file in tqdm(all_csv, desc="Sorting by increasing magnitude (from brightest to faintest)", unit="file"):
        df_all = pd.read_csv(csv_file, comment='#')

        df_all.sort_values(by='mag', inplace=True, ignore_index=True)
        nrows = len(df_all)
        tmp_file = csv_file + ".tmp"

        # Write header and sorted content to a temporary file
        with open(tmp_file, 'w', encoding='utf-8') as fw:
            fw.write(f"#Objects found: {nrows}\n")
        df_all.to_csv(tmp_file, mode='a', index=False)

        # Replace old file with sorted version
        os.remove(csv_file)
        os.rename(tmp_file, csv_file)

    print("Catalog update complete. All files sorted by magnitude (brightest to faintest).")

def convert_csv_to_parquet(dir_sc):
    """
    Convert all CSV tiles in a directory to Parquet format.

    Inputs:
        dir_sc -> [str] Path to the directory containing star catalog CSV tile files.
    Outputs:
        The newly created Parquet files.
    Notes:
        - CSV files are permanently removed after successful conversion.
    """
    # List all CSV files in the specified directory
    csv_files = glob(os.path.join(dir_sc, '*.csv'))

    for csv_file in tqdm(csv_files, desc="Converting CSV to Parquet", unit="file"):
        parquet_file = csv_file.replace(".csv", ".parquet")

        # Load CSV data
        df = pd.read_csv(csv_file,comment='#')

        # Write data to Parquet format with compression
        df.to_parquet(parquet_file, index=False, compression="zstd")

        # Delete the original CSV file
        os.remove(csv_file)