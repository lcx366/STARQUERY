import os,re,subprocess
from glob import glob
from colorama import Fore
from natsort import natsorted

NUM_TILES = 3072  # Default number of HEALPix pixels for (K=4, NSIDE=16)

def read_urls(file_path):
    """
    Reads the URL file.
    Each line of the file is like the following:
    'starcatalogs/raw/sky2000/sky2000-0.csv' 'https://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?STCS=POLYGON 0.0 90.0,0.0 87.1,45.0 84.1,90.0 87.1&format=csv&catalog=sky2000'

    Usage:
        >>> file_path = 'starcatalogs/url/gaiadr3.txt'
        >>> urls = read_urls(file_path)
    Inputs:
        file_path -> [str] The path to the URL file.
    Outputs:
        urls -> [list] A list of URLs extracted from the file.
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
    Outputs:
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
        for i in range(NUM_TILES):
            file_path = os.path.join(dir_from, f'{sc_name}-{i}.csv')

            # Display progress
            print(f'Checking {Fore.BLUE}{i + 1}{Fore.RESET} of {NUM_TILES}', end='\r')

            if file_path in files_list:
                file_size = os.path.getsize(file_path)
                dir_size += file_size / 1024  # Convert bytes to KB

                # Output the total number of rows in the star catalog file
                result = subprocess.run(['wc', '-l', file_path], capture_output=True, text=True)
                star_count = int(result.stdout.split()[0]) - 2  # Subtract the file description line and header line

                # Check the first line for object count consistency
                with open(file_path, 'r') as fn:
                    firstline = fn.readline().strip()

                expected_count = int(firstline.split(' : ')[1]) if 'Objects found' in firstline else 1
                if star_count == expected_count:
                    continue  # Valid file, skip to next

            # Log invalid or missing files for re-download
            issue_flag = True
            url = urls[i]
            outputf.write(f"'{file_path}' '{url}'\n")  # Write URL to file

    print()
    # Re-download invalid files if necessary
    if issue_flag:
        cmd = f"cat '{temp_url}' | xargs -P 16 -n 2 wget -c -O"
        os.system(cmd)
        print('Warnings: There may still be invalid files, which need to be confirmed manually.')
    else:
        print('All star catalog tile files are valid.')

    os.remove(temp_url)  # Clean up the temporary file

    # Format the directory size into a human-readable format
    dir_size_str = f'{dir_size:.2f} KB' if dir_size < 1024 else f'{dir_size / 1024:.2f} MB'
    validity = not issue_flag

    return dir_size_str, file_num, validity
