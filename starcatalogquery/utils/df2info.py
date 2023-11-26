import numpy as np

def df2info(sc_name, center, df, max_control_points, search_area, mode):
    """
    Extracts information from a specific DataFrame and generates a dictionary with star catalog data.

    Usage:
        >>> info = df2info(sc_name, center, df, max_control_points, search_area, mode)
    Inputs:
        sc_name -> [str] Name of the star catalog.
        center -> [tuple] The central coordinates of the search area in the form (Ra, Dec), where Ra is Right Ascension and Dec is Declination, both in degrees.
        df -> The DataFrame containing star data.
        max_control_points -> [int] Maximum number of stars to consider in the search, e.g., the top 100 brightest stars.
        search_area -> [float or tuple] The area of the sky being searched. For cone search, it represents the search radius; for box search, it denotes the rectangular range.
        mode -> [str] The search mode, either 'cone' or 'box'.
    Outputs:
        info -> [dict] A dictionary containing extracted star data and metadata.
    """
    
    # Determine the number of stars to consider
    num_stars = len(df)
    if max_control_points is None or max_control_points > num_stars: 
        max_control_points = num_stars
    else:
        df = df[:max_control_points]  # Limit the DataFrame to the specified maximum number of stars

    # Extract RA, Dec, and magnitude values from the DataFrame
    ra, dec = df['ra'].values, df['dec'].values
    radec = np.stack([ra, dec]).T  # Combine RA and Dec into a single array
    mag = df['mag'].values

    # Compile the extracted data into a dictionary
    info = {
        'sc_name': sc_name,
        'num_stars': num_stars,
        'df': df,
        'center': center,
        'radec': radec,
        'mag': mag,
        'max_control_points': max_control_points,
        'search_area': search_area,
        'mode': mode
    }   

    return info