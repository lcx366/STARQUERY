import numpy as np
from datetime import datetime, timezone

from .utils import data_prepare, Const
from .utils.math import spherical_to_cartesian, unit_vector, matrix_dot_vector, axes_to_antisymmetric_matrices

# Constants used for astrometric corrections
ANGLE_SUN_STAR = 5  # [degrees] Light deflection due to general relativity is considered when the angular separation between the Sun and the star is less than this value.
SOLAR_ANGULAR_RADIUS = 16 / 60  # [degrees] The apparent angular radius of the Sun.
DEFLECTION_COEFF = np.deg2rad(1.75 / 3600) * np.deg2rad(16 / 60)  # Light deflection coefficient. 1.75 arcseconds is the deflection at the Sun's limb.

def parallax2dist(theta_mas):
    """
    Convert stellar parallax (in milliarcseconds) to distance (in kiloparsecs).

    Usage:
        >>> dist_kpc = parallax2dist(1000)
    Inputs:
        theta_mas -> [float] Parallax angle in milliarcseconds (mas).
    Returns:
        dist_kpc -> [float] Distance to the object in kiloparsecs (kpc).

    Notes:
        The conversion uses the small-angle approximation:
        distance (kpc) = 1 / (parallax (radians) × AU per kiloparsec),
        where 1 kpc ≈ 206265 × 10^3 AU.
    """
    theta_rad = np.deg2rad(theta_mas / 3.6e6)
    dist_kpc = 1 / (Const.kpc_in_au * theta_rad)

    return dist_kpc

def apply_astrometry_corrections(df, astrometry_corrections, ra_rad, dec_rad):
    """
    This function applies corrections for proper motion, parallax, aberration, and light deflection to the
    right ascension (RA) and declination (Dec) of stars in a catalog based on specified parameters.

    Usage:
        >>> corrections = {'t':'2019-02-26T20:11:14.347','proper-motion':None,'aberration':(0.55952273, -1.17780654,  7.50324956),'parallax':None}
        >>> corrected_df = apply_astrometry_corrections(df, corrections, ra_rad, dec_rad)
    Inputs:
        df -> [DataFrame] Star catalog data.
        astrometry_corrections -> [dict] Dictionary specifying the types of corrections to apply.
            - 't' -> [str] Observation time in UTC, such as '2019-02-26T20:11:14.347'.
               It specifies the time at which corrections are to be applied.
            - 'proper-motion' -> [None] If present, apply proper motion correction.
               This term corrects for the motion of stars across the sky due to their velocities.
            - 'aberration' -> [tuple] Aberration correction parameters. Observer's velocity relative to Earth's center (vx, vy, vz) in km/s.
               This term corrects for the apparent shift in star positions due to the motion of the observer.
            - 'parallax' -> [None] If present, apply parallax correction.
               This term corrects for the apparent shift in star positions due to the change in observer's viewpoint as the Earth orbits the Sun.
            - 'deflection' -> [None] If present, apply light deflection correction.
               This term corrects for the bending of light from stars due to the gravitational field of the Sun, based on general relativity.
        ra_rad -> [array-like] Right Ascension values in radians.
        dec_rad -> [array-like] Declination values in radians.
    Returns:
        df -> [DataFrame] Corrected star catalog data.
    """
    # Load time system and JPL ephemeris for Earth and Sun
    ts, eph = data_prepare.ts, data_prepare.eph
    earth, sun = eph['earth'],eph['sun']
    dist = df['dist'].values * Const.kpc_in_au  # Convert distance from kpc to AU

    # Parse the time string into a datetime object and specify it as the UTC time zone
    t_ap = datetime.fromisoformat(astrometry_corrections['t']).replace(tzinfo=timezone.utc)
    t_sf = ts.utc(t_ap) # Observation time in skyfield
    epoch = ts.J(df['epoch'].values) # Epoch of star positions
    delta_year = (t_sf - epoch) / 365.25 # Time difference in years

    delta_uxyz = np.zeros((len(dist), 3))  # Initialize astrometry correction vector
    earth_t = earth.at(t_sf)
    earth_pos = earth_t.position.au  # Earth position in AU
    earth_vel = earth_t.velocity.km_per_s  # Earth velocity in km/s

    sun_t = sun.at(t_sf)
    sun_pos = sun_t.position.au  # Sun position in AU

    # Compute the direction of the Sun relative to Earth
    sun_relative_earth_pos = sun_pos - earth_pos
    sun_relative_earth_uec = unit_vector(sun_relative_earth_pos)

    # Compute the direction of the stars relative to Earth
    star_uec = spherical_to_cartesian(ra_rad, dec_rad, 1, False)
    star_relative_earth_uec = star_uec - earth_pos / dist[:, None]

    # Compute the cosine of the angle between the Sun and the stars
    sun_star_cos = np.dot(star_relative_earth_uec, sun_relative_earth_uec)

    # Apply proper motion correction
    if 'proper-motion' in astrometry_corrections:
        pm_ra = np.radians(df['pm_ra'].values / 3.6e6) # Convert proper motion from mas/yr to rad/yr
        pm_dec = np.radians(df['pm_dec'].values / 3.6e6)
        ra_rad += pm_ra * delta_year
        dec_rad += pm_dec * delta_year

    # Apply parallax correction
    if 'parallax' in astrometry_corrections:
        delta_uxyz += -earth_pos / dist[:, None]

    # Apply aberration correction
    if 'aberration' in astrometry_corrections:
        # Observer's velocity relative to Earth
        sat_relative_earth_vel = np.array(astrometry_corrections['aberration'])
        # Compute observer's velocity relative to the Solar System barycenter
        delta_uxyz += (earth_vel + sat_relative_earth_vel) / Const.light_speed

    # Apply light deflection correction
    if 'deflection' in astrometry_corrections:
        # Check if the angle between the Sun and the stars is less than the threshold (e.g., 5 degrees) and greater than the solar angular radius
        deflection_flags = (sun_star_cos > np.cos(np.deg2rad(ANGLE_SUN_STAR))) & (sun_star_cos < np.cos(np.deg2rad(SOLAR_ANGULAR_RADIUS)))
        if np.any(deflection_flags):
            _star_relative_earth_uec = star_relative_earth_uec[deflection_flags]
            _sun_star_cos = sun_star_cos[deflection_flags]
            _star_uec = star_uec[deflection_flags]
            _rotation_axis = np.cross(_star_relative_earth_uec, sun_relative_earth_uec)  # Rotation axis vector
            _rotation_uec = unit_vector(_rotation_axis)  # Rotation axis unit vector
            _antisym_matrix = axes_to_antisymmetric_matrices(_rotation_uec)  # Antisymmetric matrix from rotation axis
            _delta_theta = DEFLECTION_COEFF / np.arccos(_sun_star_cos)  # Light deflection angle
            delta_uxyz[deflection_flags] += -matrix_dot_vector(_antisym_matrix, _star_uec) * _delta_theta[:, None]

    # Compute corrections to RA and Dec
    delta_ra = (-np.sin(ra_rad) * delta_uxyz[:, 0] + np.cos(ra_rad) * delta_uxyz[:, 1]) / np.cos(dec_rad)
    delta_dec = -np.sin(dec_rad) * (np.cos(ra_rad) * delta_uxyz[:, 0] + np.sin(ra_rad) * delta_uxyz[:, 1]) + np.cos(dec_rad) * delta_uxyz[:, 2]

    # Apply corrections to DataFrame
    df['ra'] = np.rad2deg(ra_rad + delta_ra)
    df['dec'] = np.rad2deg(dec_rad + delta_dec)
    return df
