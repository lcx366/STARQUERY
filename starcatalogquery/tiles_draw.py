import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import cartopy.crs as ccrs
from matplotlib.patches import Rectangle
from astropy.visualization.wcsaxes import SphericalCircle
from astropy import units as u 

def search_draw(nside, pixels, search_area, points, fov_min):
    """
    Visualizes the search area and corresponding tile coverage on a map.

    Inputs:
        nside -> [int] HEALPix Nside parameter.
        pixels -> [list] List of HEALPix pixel indices.
        search_area -> [dict] Keyword arguments defining the search area. Can be either:
            {'cone': [ra_c, dec_c, radius]} for conical search areas, or
            {'box': [ra_min, dec_min, ra_max, dec_max]} for rectangular search areas.
        points -> [numpy array] Array of points to plot.
        fov_min -> [float] Field of view parameters in degrees.
    Outputs:
        A matplotlib plot showing the search area and tile coverages.
    """
    zoom = 400/fov_min # Adjust zoom based on the field of view
    # Handle conical search area visualization
    if 'cone' in search_area:
        [ra_c, dec_c], r = search_area['cone']
        plot_circle_PlateCarree(nside, pixels, ra_c, dec_c, r,points,zoom)
    # Handle rectangular search area visualization
    elif 'box' in search_area:
        ra_min, dec_min, ra_max, dec_max = search_area['box']
        # Check for valid latitude range
        if dec_min < -90 or dec_max > 90:
            raise Exception('Rectangle search contains poles, please replace it with cone search')
        plot_rectangle_PlateCarree(nside, pixels, ra_min, dec_min, ra_max, dec_max,points,zoom) 

def plot_circle_PlateCarree(nside, pixels, ra_c, dec_c, r,points,zoom):
    """
    Plots the coverage of tiles corresponding to a cap (conical) search area on a map.

    Inputs:
        nside -> [int] HEALPix Nside parameter.
        pixels -> [list] List of HEALPix pixel indices.
        ra_c -> [float] Right ascension of the center of the cap search area in degrees.
        dec_c -> [float] Declination of the center of the cap search area in degrees.
        r -> [float] Radius of the cap search area in degrees.
        points -> [numpy array] Array of points to plot.
        zoom -> [float] Zoom level for the map projection.
    Outputs:
        A matplotlib plot showing the cap search area and the tile coverages.
    """

    # Set the map projection
    proj = ccrs.NearsidePerspective(ra_c, dec_c, satellite_height=6.5e7/zoom)

    # Initialize and clear the figure
    plt.clf()
    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(1, 1, 1, projection=proj)

    ax.set_global()

    # Draw the cap search area
    circle = SphericalCircle((ra_c * u.deg, dec_c * u.deg), r * u.deg, fc='none', ec='m', lw=1, transform=ccrs.Geodetic())
    ax.add_patch(circle)

    # Draw each healpix tile on the map

    # Get boundary points, return theta and phi coordinates, set step=1 to get four vertices
    corners = hp.boundaries(nside, pixels)
    s1,s2,s3 = corners.shape

    lons_lats = hp.vec2ang(corners.transpose(0,2,1),lonlat=True)
    quadrangles = np.stack(lons_lats).T.reshape(s1,s3,2)

    # Loop to add quadrilaterals to the plot
    for quad in quadrangles:
        # Create a polygon, where quad is the list of vertices, closed=True means the vertices form a closed polygon
        polygon = patches.Polygon(quad, closed=True, ec='g',lw=1, transform=ccrs.Geodetic(), alpha=0.3)
        ax.add_patch(polygon)

    # Add the star distribution
    lons, lats = points.T
    ax.plot(lons, lats, marker='o', linestyle='', color='red', markersize=1, transform=ccrs.Geodetic())    

    plt.show() 

def plot_rectangle_PlateCarree(nside, pixels, ra_min, dec_min, ra_max, dec_max,points,zoom):
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

    # Set the map projection
    proj = ccrs.NearsidePerspective(ra_c, dec_c, satellite_height=6.5e7/zoom)

    # Initialize and clear the plot
    plt.clf()
    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(1, 1, 1, projection=proj)

    ax.set_global()

    # Draw the main rectangle representing the search area
    ax.add_patch(Rectangle((ra_min, dec_min), width, height, fc='none', ec='m', lw=1, transform=ccrs.Geodetic())) # ccrs.PlateCarree()

    # Draw each healpix tile on the map

    # Get boundary points, return theta and phi coordinates, set step=1 to get four vertices
    corners = hp.boundaries(nside, pixels)
    s1,s2,s3 = corners.shape

    lons_lats = hp.vec2ang(corners.transpose(0,2,1),lonlat=True)
    quadrangles = np.stack(lons_lats).T.reshape(s1,s3,2)

    # Loop to add quadrilaterals to the plot
    for quad in quadrangles:
        # Create a polygon, where quad is the list of vertices, closed=True means the vertices form a closed polygon
        polygon = patches.Polygon(quad, closed=True, ec='g',lw=1, transform=ccrs.Geodetic(), alpha=0.3)
        ax.add_patch(polygon)

    # Add the star distribution
    lons, lats = points.T
    ax.plot(lons, lats, marker='o', linestyle='', color='red', markersize=1, transform=ccrs.Geodetic())        

    plt.show()