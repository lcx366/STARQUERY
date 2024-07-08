from datetime import datetime, timedelta
from os import path, makedirs, remove

from .try_download import wget_download

def download_iers(update_interval_days=7, dir_to='~/src/iers/'):
    """
    Downloads or updates the Earth Orientation Parameters (EOP) file and Leap Seconds(LS) file from IERS.

    Usage:
        >>> dir_to, dir_eop_file, dir_leapsecond_file = download_iers()
    Inputs:
        update_interval_days -> [int, optional, default = 7] Number of days after which the IERS files should be updated if they were last modified.
        dir_to -> [str, optional, default = '~/src/iers/'] Directory for storing IERS files.
    Outputs:
        dir_to -> [str] Directory of the IERS files.
        dir_eop_file -> [str] Path of the EOP file.
        dir_leapsecond_file -> [str] Path of the Leap Second file.
    """
    # Expanding the user directory
    dir_to = path.expanduser(dir_to)

    # Define filenames and URLs
    eop_file = 'finals2000A.all'
    leapsecond_file = 'Leap_Second.dat'
    dir_eop_file = path.join(dir_to, eop_file)
    dir_leapsecond_file = path.join(dir_to, leapsecond_file)

    url_eop = 'https://datacenter.iers.org/products/eop/rapid/standard/finals2000A.all'
    url_leapsecond = 'https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat'

    # Create the directory if it does not exist
    if not path.exists(dir_to):
        makedirs(dir_to)

    # Download or update the EOP file
    if not path.exists(dir_eop_file):
        desc = f"Downloading the latest EOP file '{eop_file}' from IERS."
        wget_download(url_eop, dir_eop_file, desc)
    else:
        # Check if update is needed based on the modification time
        modified_time = datetime.fromtimestamp(path.getmtime(dir_eop_file))
        if datetime.now() > modified_time + timedelta(days=update_interval_days):
            remove(dir_eop_file)
            desc = f"Updating the EOP file '{eop_file}' from IERS."
            wget_download(url_eop, dir_eop_file, desc)
        else:
            print(f"The EOP file '{eop_file}' in {dir_to} is already the latest.")

    # Download or update the Leap Second file
    if not path.exists(dir_leapsecond_file):
        desc = f"Downloading the latest Leap Second file '{leapsecond_file}' from IERS."
        wget_download(url_leapsecond, dir_leapsecond_file, desc)
    else:
        # Check if update is needed based on the modification time
        modified_time = datetime.fromtimestamp(path.getmtime(dir_leapsecond_file))
        if datetime.now() > modified_time + timedelta(days=update_interval_days):
            remove(dir_leapsecond_file)
            desc = f"Updating the Leap Second file '{leapsecond_file}' from IERS."
            wget_download(url_leapsecond, dir_leapsecond_file, desc)
        else:
            print(f"The Leap Second file '{leapsecond_file}' in {dir_to} is already the latest.")

    return dir_to, dir_eop_file, dir_leapsecond_file

def download_sspe(jpleph, dir_to='~/src/ephem/'):
    """
    Downloads the JPL Solar System Planetary Ephemeris(SSPE) file from NAIF.

    Usage:
        >>> ephem_file = download_sspe('DE440S')
    Inputs:
        jpleph -> [str] Version of the JPL ephemeris, such as 'DE430' or 'DE440S'.
        dir_to -> [str, optional, default = '~/src/ephem/'] Directory for storing the ephemeris file.
    Outputs:
        ephem_file -> [str] Path of the downloaded JPL ephemeris file.
    """
    # Expanding the user directory
    dir_to = path.expanduser(dir_to)

    file = jpleph + '.bsp'
    dir_file = path.join(dir_to, file)
    url = 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/' + file

    # Create the directory if it does not exist
    if not path.exists(dir_to):
        makedirs(dir_to)

    # Download the JPL ephemeris file if it does not exist
    if not path.exists(dir_file):
        desc = f"Downloading the JPL Solar System Planetary Ephemeris file '{file}' from NAIF."
        try:
            wget_download(url, dir_file, desc)
        except Exception as e:
            print(f"Failed to download from NAIF, error message: {str(e)}.\n"
                  f"Please try a backup URL and manually place a pre-prepared {file} file into '{dir_to}'.")
    else:
        print(f"The JPL Solar System Planetary Ephemeris file '{file}' is already in {dir_to}.")

    return dir_file
