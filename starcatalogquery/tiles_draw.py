import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.patches import Rectangle
from astropy.visualization.wcsaxes import SphericalCircle
from astropy import units as u

from .catalog_query import cone2seqs,box2seqs,seq2radec  

def search_draw(tile_size, kwargs):
    """
    Visualizes the search area and corresponding tile coverage on a map.

    Usage:
        >>> kwargs_box = {'box':[20,30,30,40]}
        >>> kwargs_cone = {'cone':[20,30,10]}
        >>> search_draw(3,kwargs_box)
        >>> search_draw(3,kwargs_cone)
    Inputs:
        tile_size -> [int] Size of each tile in degrees.
        kwargs -> [dict] Keyword arguments defining the search area. Can be either:
            {'cone': [ra_c, dec_c, radius]} for conical search areas, or
            {'box': [ra_min, dec_min, ra_max, dec_max]} for rectangular search areas.
    """
    # Handle conical search area visualization
    if 'cone' in kwargs:
        # Extract cone parameters and calculate corresponding tile sequences
        ra_c, dec_c, r = kwargs['cone']
        seqs = cone2seqs(ra_c, dec_c, r, tile_size)
        radec_boxes, box_center = seq2radec(seqs, tile_size)
        plot_circle_PlateCarree(radec_boxes, ra_c, dec_c, r)

    # Handle rectangular search area visualization
    elif 'box' in kwargs:
        # Extract box parameters and calculate corresponding tile sequences
        ra_min, dec_min, ra_max, dec_max = kwargs['box']
        # Check for valid latitude range
        if dec_min < -90 or dec_max > 90:
            raise Exception('Rectangle search contains poles, please replace it with cone search')
        seqs = box2seqs(kwargs['box'], tile_size)
        radec_boxes, box_center = seq2radec(seqs, tile_size)
        plot_rectangle_PlateCarree(radec_boxes, ra_min, dec_min, ra_max, dec_max) 

def plot_circle_PlateCarree(radec_boxes, ra_c, dec_c, r):
    """
    Plots the coverage of tiles corresponding to a cap (conical) search area on a map.

    Usage:
        >>> plot_circle_PlateCarree(radec_boxes, ra_c, dec_c, r)
    Inputs:
        radec_boxes -> [2D array-like] Array of [ra_min, dec_min, ra_max, dec_max] for each tile.
        ra_c -> [float] Right ascension of the center of the cap search area in degrees.
        dec_c -> [float] Declination of the center of the cap search area in degrees.
        r -> [float] Radius of the cap search area in degrees.
    Outputs:
        A matplotlib plot showing the cap search area and the tile coverages.
    """
    # Define the extended region around the search area to display
    region = left_ra, right_ra, lower_dec, upper_dec = ra_c - 1.5 * r, ra_c + 1.5 * r, dec_c - 1.5 * r, dec_c + 1.5 * r

    # Define spans for RA and DEC based on the radius
    ra_span = dec_span = 3 * r

    # Set the map projection with Plate Carree and centering longitude at ra_c
    proj = ccrs.PlateCarree(central_longitude=ra_c)

    # Initialize and clear the figure
    plt.clf()
    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(1, 1, 1, projection=proj)

    # Set the map to cover the defined region
    ax.set_extent(region, crs=ccrs.PlateCarree())

    # Calculate grid line locations, adjust for right ascension exceeding 180 degrees
    xlocs = np.linspace(left_ra,right_ra,find_divisors(int(ra_span))+1)
    ylocs = np.linspace(lower_dec,upper_dec,find_divisors(int(dec_span))+1)
    if right_ra > 180:
        block_flags = xlocs > 180
        xlocs[block_flags] -= 360

    # Add grid lines and label styles
    gl = ax.gridlines(xlocs=xlocs, ylocs=ylocs, draw_labels=True, linestyle='--', alpha=0.7)
    gl.xlabel_style = {'size': 7, 'color': 'gray'}
    gl.ylabel_style = {'size': 7, 'color': 'gray'}

    # Draw the cap search area
    circle = SphericalCircle((ra_c * u.deg, dec_c * u.deg), r * u.deg, fc='none', ec='m', lw=1.5, transform=ccrs.Geodetic())
    ax.add_patch(circle)

    # Draw each tile as a rectangle on the map
    for radec_box in radec_boxes:
        ra_min, dec_min, ra_max, dec_max = radec_box
        width = ra_max - ra_min
        height = dec_max - dec_min
        ax.add_patch(Rectangle((ra_min, dec_min), width, height, ec='g', lw=1, transform=ccrs.PlateCarree(), alpha=0.3))

    # Display the plot
    plt.show() 

def plot_rectangle_PlateCarree(radec_boxes, ra_min, dec_min, ra_max, dec_max):
    """
    Plots a rectangular search area and its tile coverages on a map.

    Usage:
        >>> lot_rectangle_PlateCarree(radec_boxes, ra_min, dec_min, ra_max, dec_max)
    Inputs:
        radec_boxes -> [array-like] Array of [ra_min, dec_min, ra_max, dec_max] for each tile.
        ra_min -> [float] Minimum right ascension of the search area in degrees.
        dec_min -> [float] Minimum declination of the search area in degrees.
        ra_max -> [float] Maximum right ascension of the search area in degrees.
        dec_max -> [float] Maximum declination of the search area in degrees.
    Outputs:
        A matplotlib plot displaying the rectangular search area and tile coverages.
    """
    # Calculate center and span for the search area
    ra_c, dec_c = (ra_min + ra_max) / 2, (dec_min + dec_max) / 2
    width, height = ra_max - ra_min, dec_max - dec_min
    ra_span, dec_span = width * 1.5, height * 1.5

    # Define the extended region to display
    region = ra_c - ra_span / 2, ra_c + ra_span / 2, dec_c - dec_span / 2, dec_c + dec_span / 2

    # Set map projection with central longitude at the center of the search area
    proj = ccrs.PlateCarree(central_longitude=ra_c)

    # Initialize and clear the plot
    plt.clf()
    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(1, 1, 1, projection=proj)

    # Set the extent of the map to cover the defined region
    ax.set_extent(region, crs=ccrs.PlateCarree())

    # Calculate grid line locations
    xlocs = np.linspace(region[0], region[1], find_divisors(int(ra_span)) + 1)
    ylocs = np.linspace(region[2], region[3], find_divisors(int(dec_span)) + 1)

    # Adjust for right ascension exceeding 180 degrees
    if region[1] > 180:
        xlocs[xlocs > 180] -= 360

    # Add grid lines and set label styles
    gl = ax.gridlines(xlocs=xlocs, ylocs=ylocs, draw_labels=True, linestyle='--', alpha=0.7)
    gl.xlabel_style = {'size': 7, 'color': 'gray'}
    gl.ylabel_style = {'size': 7, 'color': 'gray'}

    # Draw the main rectangle representing the search area
    ax.add_patch(Rectangle((ra_min, dec_min), width, height, fc='none', ec='m', lw=1.5, transform=ccrs.PlateCarree()))

    # Draw each tile as a rectangle on the map
    for radec_box in radec_boxes:
        ra_min, dec_min, ra_max, dec_max = radec_box
        width = ra_max - ra_min
        height = dec_max - dec_min
        ax.add_patch(Rectangle((ra_min, dec_min), width, height, ec='g', lw=1, transform=ccrs.PlateCarree(), alpha=0.3))    

    # Display the plot
    plt.show()

def find_divisors(x):
    """
    Calculates the optimal grid line spacing for regional maps based on the area span.

    Usage:
        >>> print(find_divisors(20)) 
        >>> # 4
        >>> print(find_divisors(25)) 
        >>> # 5
    Inputs:
        x -> [int] Area span in degrees (RA or DEC). Should not be a prime number.
    Outputs:
        y -> [int] Optimal grid line spacing.
    """

    # Compute all common divisors of x
    y = np.unique(np.gcd(np.arange(x), x))
    n = len(y)

    # Return the middle value for an odd number of divisors, or the one before middle for even
    return y[n // 2] if n % 2 == 1 else y[n // 2 - 1]  