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
    Loads the Solar System Planetary Ephemeris (SSPE) file from NAIF (NASA's Navigation and Ancillary Information Facility)
    using the Skyfield library. This function is responsible for downloading and loading a specific JPL ephemeris into
    the global variable `eph`, which is later used to compute positions and motions of celestial objects, particularly
    the Earth and the Sun.

    Usage:
        >>> sspe_load('DE440S')

    Inputs:
        jpleph -> [str]: The name of the desired JPL ephemeris. Common values include:
                          'DE430' - JPL Developmental Ephemeris 430
                          'DE440' - JPL Developmental Ephemeris 440
                          'DE440S' - A simplified version of DE440

    Outputs:
        This function does not explicitly return any value. Instead, it sets the global variable `eph`,
        which contains the ephemeris data for computing celestial positions.
    """
    global eph  # Global variable to store the ephemeris data

    jpleph = jpleph.lower()  # Convert input to lowercase for case-insensitive matching

    # Check if the input ephemeris is one of the supported types
    if jpleph in ['de430', 'de440', 'de440s']:
        ephem_file = download_sspe(jpleph)  # Download the ephemeris file
        eph = load(ephem_file)  # Load the downloaded ephemeris file into the global variable `eph`
    else:
        # Raise an exception if the ephemeris name is not supported
        raise Exception(
            f"The specified ephemeris '{jpleph}' is not supported. Available ephemerides are:\n- DE430\n- DE440\n- DE440S")
