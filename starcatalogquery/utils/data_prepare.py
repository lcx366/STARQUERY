from skyfield.api import Loader, load
from skyfield.data import iers as iers_skyfield
from astropy.utils import iers as iers_astropy

from .data_download import download_iers, download_sspe

def iers_load():
    """
    Loads the Earth Orientation Parameters(EOP) and Leap Seconds(LS) files from IERS.
    These files are essential for accurate time and coordinate transformations in astronomical calculations.
    This function downloads the necessary files if they are not found locally, and then sets up both Skyfield and Astropy libraries to use this data.

    Outputs -> Global variable to store the IERS time system
    """
    global ts  # Global variable to store the IERS time system

    # Download the EOP and LS files
    dir_iers, eop_file, leapsecond_file = download_iers()

    # Load the IERS data for Skyfield
    load_iers = Loader(dir_iers)
    ts = load_iers.timescale(builtin=False)  # Load the IERS time system
    with load.open(eop_file) as f:
        pm_data = iers_skyfield.parse_x_y_dut1_from_finals_all(f)  # Parse the EOP file
    iers_skyfield.install_polar_motion_table(ts, pm_data)  # Load the EOP data into Skyfield

    # Load the IERS data for Astropy
    iers_astropy.conf.auto_download = False  # Prevent automatic IERS data download by Astropy
    iers_a = iers_astropy.IERS_A.open(eop_file)  # Load the EOP data into Astropy
    iers_astropy.LeapSeconds.from_iers_leap_seconds(leapsecond_file)  # Load the Leap Seconds data into Astropy
    iers_astropy.earth_orientation_table.set(iers_a)  # Configure Astropy to use the IERS data

def sspe_load(jpleph):
    """
    Loads the Solar System Planetary Ephemeris (SPPE) file from NAIF.

    Usage:
        >>> sspe_load('DE440S')
    Inputs:
        jpleph -> [str] Name of the JPL ephemeris, e.g. 'DE440S'.
    Outputs:
        earth, sun -> Global variables to store the ephemerides of the Earth and Sun
    """
    global earth, sun  # Global variables to store the ephemerides of the Earth and Sun

    jpleph = jpleph.lower() # Convert to lowercase for consistency

    if jpleph in ['de430', 'de440', 'de440s']:
        ephem_file = download_sspe(jpleph)  # Download the SPPE file
        eph = load(ephem_file)  # Load the SPPE file into Skyfield
        earth,sun = eph['earth'],eph['sun']
    else:
        raise Exception(
            f"The specified ephemeris '{jpleph}' is not supported. Available ephemerides are:\n- DE430\n- DE440\n- DE440S")
