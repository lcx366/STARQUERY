import numpy as np

def df2info(sc_name, center, df, max_control_points, level, nside, pixels_ids, pixel_size, search_area):
    """
    Extracts information from a star catalog dataFrame.

    Usage:
        >>> info = df2info(sc_name, center, df, max_control_points, level, nside, pixels_ids, pixel_size, search_area)
    Inputs:
        sc_name -> [str] Name of the star catalog, e.g. 'gaiadr3'.
        center -> [tuple] The central coordinates of the search area in form of (Ra, Dec), where Ra is Right Ascension and Dec is Declination, both in degrees.
        df -> [Pandas dataframe] The star catalog data frame.
        max_control_points -> [int] Maximum number of stars to consider in the search, e.g., the top 100 brightest stars.
        level -> [str] Partition level, such as 'K4'. It indicates the hierarchical level of the HEALPix partitioning.
        nside -> [int] The number of divisions along the side of a HEALPix base pixel. It determines the resolution of the HEALPix grid.
        pixels_ids -> [array-like] The HEALPix pixel IDs covering the search area.
        pixel_size -> [float] Approximate size of each pixel in degrees.
        search_area -> [dict] The area of the sky being searched. For cone search, it denotes the center and search radius; for box search, it denotes the rectangular range.
    Outputs:
        info -> [dict] A dictionary containing the metadata and star data(Ra,Dec,Mag) extracted from the star catalog data frame.
    """

    # Determine the number of stars to consider
    num_stars = len(df)
    if max_control_points is None or max_control_points > num_stars:
        max_control_points = num_stars
    else:
        df = df[:max_control_points].copy()  # Limit the DataFrame to the specified maximum number of stars

    # Remove duplicate stars based on 'ra' and 'dec' to avoid considering variable stars
    df.drop_duplicates(subset=['ra', 'dec'], keep=False, inplace=True)

    # Extract RA, Dec, and magnitude values from the DataFrame
    ra, dec = df['ra'].values, df['dec'].values
    radec = np.stack([ra, dec]).T  # Combine RA and Dec into a single array
    mag = df['mag'].values

    # Compile the extracted data into a dictionary
    info = {
        'sc_name': sc_name,
        'num_stars': num_stars,
        'df': df,
        'radec': radec,
        'mag': mag,
        'max_control_points': max_control_points,
        'center': center,
        'level': level,
        'nside': nside,
        'tiles_ids': pixels_ids,
        'tile_size': '{:.2f} deg'.format(pixel_size),
        'search_area': search_area
    }

    return info
