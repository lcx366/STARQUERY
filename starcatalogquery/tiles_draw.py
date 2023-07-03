import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.patches import Rectangle
from astropy.visualization.wcsaxes import SphericalCircle
from astropy import units as u

from .catalog_query import cone2seqs,box2seqs,seq2radec  

def search_draw(tile_size,kwargs):
    """
    Visualize the scope of the search area and the coverage of the corresponding tiles.

    Usage:
        >>> kwargs_box = {'box':[20,30,30,40]}
        >>> kwargs_cone = {'cone':[20,30,10]}
        >>> search_draw(3,kwargs_box)
        >>> search_draw(3,kwargs_cone)

    Inputs:
        tile_size -> [int] Size of the tile in [deg]
        kwargs -> [dict] Keyword arguments, which defines the scope of the search area, such as {'cone':[20,30,10]} or {'box':[20,30,30,40]}, where
            for {'box':[ra_min,dec_min,ra_max,dec_max]}:
                ra_min -> [float] Left border of RA in [deg].
                dec_min -> [float] Lower border of DEC in [deg].
                ra_max -> [float] Right border of RA in [deg].
                dec_max -> [float] Upper border of DEC in [deg].
            for {'cone':[ra_c,dec_c,radius]}:
                ra_c -> [float] RA, in [deg].
                dec_c -> [float] DEC, in [deg].
                radius -> [float] Angular radius of the cap, in [deg].  

        Outputs:
            An image that shows the scope of the search area and the coverage of the corresponding tiles.          
    """
    if 'cone' in kwargs:
        ra_c,dec_c,r = kwargs['cone']
        seqs = cone2seqs(ra_c,dec_c,r,tile_size)
        radec_boxes,box_center = seq2radec(seqs,tile_size)
        plot_circle_PlateCarree(radec_boxes,ra_c,dec_c,r)
    elif 'box' in kwargs:    
        ra_min,dec_min,ra_max,dec_max = kwargs['box']
        if dec_min < -90 or dec_max > 90: raise Exception('Rectangle search contains poles, please replace it with cone search')
        seqs = box2seqs(kwargs['box'],tile_size)
        radec_boxes,box_center = seq2radec(seqs,tile_size)
        plot_rectangle_PlateCarree(radec_boxes,ra_min,dec_min,ra_max,dec_max)    

def plot_circle_PlateCarree(radec_boxes,ra_c,dec_c,r):
    """
    Plot the scope of the cap search area and the coverage of the corresponding tiles.         
    """
    region = left_ra, right_ra, lower_dec, upper_dec = ra_c-1.5*r, ra_c+1.5*r, dec_c-1.5*r, dec_c+1.5*r 
    ra_span = dec_span = 3*r
        
    # set the map projection
    proj = ccrs.PlateCarree(central_longitude=ra_c)
   
    # plot
    plt.clf()
    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(1, 1, 1,projection = proj)

    # set the range of the area of interest
    ax.set_extent(region,crs = ccrs.PlateCarree())
    xlocs = np.linspace(left_ra,right_ra,med(int(ra_span))+1)
    ylocs = np.linspace(lower_dec,upper_dec,med(int(dec_span))+1)
        
    if right_ra > 180:
        block_flags = xlocs > 180
        xlocs[block_flags] -= 360

    gl = ax.gridlines(xlocs = xlocs, ylocs = ylocs,draw_labels=True,linestyle='--',alpha=0.7)
    gl.xlabel_style = {'size': 7, 'color': 'gray'}
    gl.ylabel_style = {'size': 7, 'color': 'gray'}
    circle = SphericalCircle((ra_c*u.deg, dec_c*u.deg), r*u.deg,fc='none',ec='m',lw=1.5,transform=ccrs.Geodetic())
    ax.add_patch(circle)
    
    for radec_box in radec_boxes:
        ra_min,dec_min,ra_max,dec_max = radec_box
        width = ra_max - ra_min
        height = dec_max - dec_min
    
        ax.add_patch(Rectangle((ra_min,dec_min),width, height, ec ='g',lw = 1,transform=ccrs.PlateCarree(),alpha=0.3))

    plt.show()   

def plot_rectangle_PlateCarree(radec_boxes,ra_min,dec_min,ra_max,dec_max):
    """
    Plot the scope of the rectangle search area and the coverage of the corresponding tiles.         
    """
    ra_c = (ra_min + ra_max)/2
    dec_c = (dec_min + dec_max)/2
    width = ra_max - ra_min
    height = dec_max - dec_min
    ra_span = width*1.5
    dec_span = height*1.5
    region = left_ra, right_ra, lower_dec, upper_dec = ra_c - ra_span/2, ra_c + ra_span/2, dec_c - dec_span/2, dec_c + dec_span/2

    # set the map projection
    proj = ccrs.PlateCarree(central_longitude=ra_c)
   
    # plot
    plt.clf()
    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(1, 1, 1,projection = proj)

    # set the range of the area of interest
    ax.set_extent(region,crs = ccrs.PlateCarree())
    xlocs = np.linspace(left_ra,right_ra,med(int(ra_span))+1)
    ylocs = np.linspace(lower_dec,upper_dec,med(int(dec_span))+1)
        
    if right_ra > 180:
        block_flags = xlocs > 180
        xlocs[block_flags] -= 360

    gl = ax.gridlines(xlocs = xlocs, ylocs = ylocs,draw_labels=True,linestyle='--',alpha=0.7)
    gl.xlabel_style = {'size': 7, 'color': 'gray'}
    gl.ylabel_style = {'size': 7, 'color': 'gray'}
    ax.add_patch(Rectangle((ra_min,dec_min),width,height, fc ='none',ec ='m',lw = 1.5,transform=ccrs.PlateCarree()))

    for radec_box in radec_boxes:
        ra_min,dec_min,ra_max,dec_max = radec_box
        width = ra_max - ra_min
        height = dec_max - dec_min
    
        ax.add_patch(Rectangle((ra_min,dec_min),width,height, ec ='g',lw = 1,transform=ccrs.PlateCarree(),alpha=0.3))

    plt.show()      

def med(x):
    """
    Calculate the optimal grid line spacing for regional maps.

    Usage: 
        >>> y = med(x)
    
    Inputs:
        x -> [int] Area span in RA or DEC in [deg]. It cannot be a prime number.

    Outputs:
        y -> [int] The optimal grid line spacing

    Examples:
        >>> print(med(20)) # 4
        >>> print(med(25)) # 5
    """
    y = np.unique(np.gcd(np.arange(x),x))
    n = len(y)
    if n%2 == 1:
        return y[n//2]
    else:
        return y[n//2-1]    