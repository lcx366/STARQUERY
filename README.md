# Welcome to the STARQUERY package

[![PyPI version shields.io](https://img.shields.io/pypi/v/starcatalogquery.svg)](https://pypi.python.org/pypi/starcatalogquery/) [![PyPI pyversions](https://img.shields.io/pypi/pyversions/starcatalogquery.svg)](https://pypi.python.org/pypi/starcatalogquery/) [![PyPI status](https://img.shields.io/pypi/status/starcatalogquery.svg)](https://pypi.python.org/pypi/starcatalogquery/) [![GitHub contributors](https://img.shields.io/github/contributors/lcx366/STARQUERY.svg)](https://GitHub.com/lcx366/STARQUERY/graphs/contributors/) [![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/lcx366/STARQUERY/graphs/commit-activity) [![GitHub license](https://img.shields.io/github/license/lcx366/STARQUERY.svg)](https://github.com/lcx366/STARQUERY/blob/master/LICENSE) [![Documentation Status](https://readthedocs.org/projects/starcatalogquery/badge/?version=latest)](http://starcatalogquery.readthedocs.io/?badge=latest) [![Build Status](https://travis-ci.org/lcx366/starcatalogquery.svg?branch=master)](https://travis-ci.org/lcx366/starcatalogquery)

STARQUERY is an archive of scientific routines for establishing an offline star catalog query database.

## Key Features

1. **Offline Star Catalog Database Creation:** Allows users to import extensive star catalog data from the [STScI Outerspace](https://outerspace.stsci.edu/display/GC/WebServices+for+Catalog+Access) and [The Astronomy Nexus](https://www.astronexus.com/hyg) astronomical sources, creating a detailed and accessible local offline database.
2. **Data Simplification**: Condenses the star catalog data, ensuring minimal storage space usage while retaining essential information.
3. **Advanced Search and Query**: Offers sophisticated search tools, including rectangle search and spherical cap (cone) search, for precise star filtering based on location, magnitude, and more.
4. **Search Area Visualization**: Visualizes the scope of the search area and the coverage of corresponding catalog tiles.
5. **Pixel Coordinate Calculation**: Given a pixel width, calculates the pixel coordinates of filtered stars.
6. **Feature Libraries Construction**: Generates geometrically invariant features based on the spatial configuration of stars, aiding in star pattern recognition.
7. **HEALPix-based Sky Area Division**: Divide the celestial sphere into multiple sky areas with the HEALPix algorithm.
8. **Local Offline Star Catalog Database Loading**: Enables loading of the local offline star catalog database for efficient data access and manipulation.

## How to Install

Install STARQUERY with ease using the following pip command in your terminal:

```
pip install starcatalogquery
pip install starcatalogquery --upgrade # to upgrade a pre-existing installation
```

### Dependencies

STARQUERY relies on several external libraries and packages to function optimally. Below is a list of these dependencies with a brief description of their role in the project:

#### healpy

It is a Python package to handle pixelated data on the sphere. It is based on HEALPix (Hierarchical Equal Area isoLatitude Pixelization), a versatile scheme used for the pixelization of data on the sphere. In STARQUERY, healpy is primarily used for dividing the celestial sphere into multiple equal-area sky regions, a critical step in the blind matching of star maps.

*healpy* can be easily installed via pip. Run the following command to install:

```
pip install healpy
```

For more information on healpy and its functionalities, refer to the [healpy documentation](https://healpy.readthedocs.io/en/latest/).

## How to use

### Build an Offline Star Catalog Database

Creating a local offline star catalog with STARQUERY is straightforward and efficient. You can obtain star catalog data from renowned astronomical data repositories like the [STScI Outerspace](https://outerspace.stsci.edu/display/GC/WebServices+for+Catalog+Access) and [The Astronomy Nexus](https://www.astronexus.com/hyg). STARQUERY supports a wide range of star catalogs including 'hygv3.7', 'at-hygv2.4', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', and more.

To start building your database, use the `StarCatalog` class in the STARQUERY package. For example, to download the HYGv3.7 star catalog with a tile size of 5 degrees, you can use the following Python script:

```python
>>> from starcatalogquery import StarCatalog
>>> # Get the star catalog HYGv3.7 with the tile size of 5 deg
>>> hyg_raw = StarCatalog.get('hygv3.7',5)
>>> print(hyg_raw) # Display basic information about the downloaded catalog
>>> # <StarCatalogRaw object: CATALOG_NAME = 'hygv3.7' CATALOG_SIZE = '31.3 MB' TILES_NUM = 2592 TILE_SIZE = '5 deg' STARS_NUM = '~ 0.12 Million' MAG = '< 9'>
```

This script initiates the download and construction of your offline database. You can customize the tile size as per your requirements, and STARQUERY will handle the rest, seamlessly integrating the data into your local storage for quick and easy access.

### Reduce the Raw Star Catalogs

The original star catalog, while comprehensive, often contains extensive information that results in substantial storage requirements and can impede the speed of star queries. To address this, STARQUERY provides a feature to streamline the catalog by extracting only the essential star information. This is achieved using the .reduce() method, which constructs a more compact version of the star catalog database.

```python
>>> # Reduce the raw star catalog to its essential components
>>> hyg_reduced = hyg_raw.reduce()
>>> print(hyg_reduced) # Display basic information about the reduced catalog
>>> # <StarCatalogReduced object: CATALOG_NAME = 'hygv3.7' CATALOG_SIZE = '5.2 MB' TILES_NUM = 2592 TILE_SIZE = '5 deg' STARS_NUM = '~ 0.12 Million' MAG = '< 9'>
```

The reduced star catalog will include only the celestial position, proper motion, apparent magnitude, and epoch of stars, significantly reducing the storage space and improving query performance.

### Simplify the Reduced Star Catalogs

To further optimize the size of the star catalog and enhance its query efficiency, STARQUERY allows users to filter the reduced star catalog based on the limit magnitude of the detector. In addition, it facilitates proper motion corrections to produce a minimalist star catalog. This process significantly trims down the data, focusing on the most relevant astronomical information.

```python
>>> # Set the magnitude threshold and target epoch for proper motion correction
>>> mag_threshold,t_pm = 9.0,2022.0 
>>> # Simplify the reduced star catalog
>>> hyg_simplified = hyg_reduced.simplify(mag_threshold,t_pm)
>>> print(hyg_simplified) # Display basic information about the simplified catalog
>>> # <StarCatalogSimplified object: CATALOG_NAME = 'hygv3.7' CATALOG_SIZE = '2.0 MB' TILES_NUM = 2592 TILE_SIZE = '5 deg' STARS_NUM = '~ 0.12 Million' MAG = '< 9' MAG_CUTOFF = 9.0 EPOCH = 2022.0>
```

This minimalist star catalog provides just the essential data: the celestial position and apparent magnitude of stars at a specific epoch. This streamlined dataset is ideal for applications where speed and efficiency are paramount.

### Query Information About Stars in a Specific Sky Area

STARQUERY allows users to perform a cone search in the raw, reduced, or simplified star catalog to find stars within a specific area of the sky. This feature is particularly useful for astronomers who need to focus on a particular region for their observational studies or data analysis.

```python
>>> # Set the center pointing in form of [Ra, Dec] in degrees and search radius in degrees
>>> center,radius = [20,30],15 
>>> # Set the cutoff magnitude and the observation epoch
>>> mag_threshold,t_pm = 9.0,2022.0
>>> # Optionally, set the maximum number of brightest stars to output
>>> max_num = 100
>>> # Perform a cone search on the raw star catalog
>>> hyg_raw_stars = hyg_raw.search_cone(center,radius,mag_threshold,t_pm,max_num)
>>> # Perform a cone search on the reduced star catalog
>>> hyg_reduced_stars = hyg_reduced.search_cone(center,radius,mag_threshold,t_pm,max_num)
>>> hyg_simplified_stars = hyg_simplified.search_cone(center,radius,max_num)
>>> print(hyg_simplified_stars)
>>> # <Stars object: CATALOG = 'hygv3.7' CENTER = '[20, 30] in [RA,DEC]' SEARCH_AREA = '15 deg' STARS_NUM = 1345 MCP = 100 MODE = 'CONE'>
```

This functionality enables the retrieval of detailed information about stars within a defined radius from a given center point in the celestial coordinates. The results can be further refined by specifying a magnitude threshold, observation epoch, and the number of brightest stars to be included in the output.

STARQUERY also provides the capability to perform a rectangle search in the raw, reduced, and simplified star catalogs. This function is useful for users who wish to focus on a specific rectangular area in the celestial coordinates, enabling targeted astronomical analysis and data collection.

```python
>>> # Set a rectangular search area in form of [ra_min, dec_min, ra_max, dec_max] in degrees
>>> radec_box = [5,15,35,45]
>>> # Set the cutoff magnitude and the observation epoch
>>> mag_threshold,t_pm = 9.0,2022.0
>>> # Optionally, set the maximum number of brightest stars to output
>>> max_num = 100
>>> # Perform a rectangle search on the raw star catalog
>>> hyg_raw_stars = hyg_raw.search_box(radec_box,mag_threshold,t_pm,max_num)
>>> # # Perform a rectangle search on the reduced star catalog
>>> hyg_reduced_stars = hyg_reduced.search_box(radec_box,mag_threshold,t_pm,max_num)
>>> # Perform a rectangle search on the simplified star catalog
>>> hyg_simplified_stars = hyg_simplified.search_box(radec_box,max_num)
>>> print(hyg_simplified_stars)
>>> # <Stars object: CATALOG = 'hygv3.7' CENTER = '[20.0, 30.0] in [RA,DEC]' SEARCH_AREA = '[5, 15, 35, 45] deg' STARS_NUM = 1455 MCP = 100 MODE = 'BOX'>
```

This functionality allows for the retrieval of star data within a specified rectangular area, offering flexibility in terms of magnitude cutoff, epoch, and number of stars to include. 

### Calculate the Pixel Coordinates of the Filtered Stars

STARQUERY enables users to calculate the pixel coordinates of filtered stars, leveraging the TAN projection in World Coordinate System (WCS) transformations. This feature is particularly useful in astronomical imaging and analysis, where precise pixel-level data is crucial.

```python
>>> # Specify the pixel width of the detector in degrees
>>> pixel_width = 0.01
>>> # Calculate the pixel coordinates using the TAN projection
>>> hyg_simplified_stars.pixel_xy(0.01)
>>> # Retrieve the calculated pixel coordinates
>>> xy = hyg_simplified_stars.xy
>>> print(xy)
```

This functionality is essential for converting celestial coordinates into pixel coordinates, which is a key step in many astronomical applications, particularly in the field of astrophotography and celestial navigation.

### Calculate the Triangle Invariants and the Asterism Indices of the Filtered Stars

This functionality in STARQUERY, which focuses on the calculation of triangle invariants and asterism indices, draws inspiration from the Astroalign software package developed by Beroiz, M. I. (https://astroalign.quatrope.org/en/latest/). We acknowledge  the work done in creating Astroalign, a Python package that has advanced the field of star map matching.

Steps:

1. **Derive Geometric Invariants:** For each possible triangle formed by groups of three stars, the software calculates unique geometric invariants (ratios L2/L1 and L1/L0), where L2, L1, and L0 are the sides of each triangle sorted in descending order.
2. **Construct a 2D-Tree Structure:** Utilizing these invariants, a 2D-Tree structure is constructed to facilitate efficient spatial queries and pattern matching in star maps.
3. **Associate Invariants with Star Indices:** Each invariant set is linked to the indices of the stars forming the corresponding triangle, thereby maintaining a connection between the geometric features and the stars.

Note that this method should be employed after the calculation of pixel coordinates.

```python
>>> # Generate the invariant features for star map matching
>>> hyg_simplified_stars.invariantfeatures()
>>> # Retrieve the generated invariants, asterisms, and 2D-Tree structure
>>> invariants = hyg_simplified_stars.invariants
>>> asterisms = hyg_simplified_stars.asterisms
>>> kdtree = hyg_simplified_stars.kdtree
>>> print(invariants,asterisms,kdtree)
```

This advanced feature of STARQUERY greatly enhances the capability to match and analyze star patterns.

### Visualization

STARQUERY provides a visualization feature to display the scope of the search area and the coverage of the corresponding catalog tiles. This is particularly useful for understanding the spatial extent of your search and the distribution of tiles within the catalog.

```python
>>> # Define a rectangular search area
>>> box_area = {'box':[5,15,35,45]} # Format: {'box': [ra_min, dec_min, ra_max, dec_max]}
>>> # Define a cone search area
>>> cone_area = {'cone':[20,30,15]} # Format: {'cone': [ra_c, dec_c, radius]}
>>> # Visualize the rectangle search area
>>> hyg_simplified._search_draw(box_area)
>>> # Visualize the cone search area
>>> hyg_simplified._search_draw(cone_area)
>>> # Note: ._search_draw is also applicable for hyg_raw and hyg_reduced
```

<p align="middle">
  <img src="readme_figs/box.png" width="500" />
</p>

<p align="middle">
  <img src="readme_figs/cone.png" width="500" />
</p>

### Divide the Sky into Multiple Equal-Area Sky Regions

To facilitate blind matching of star maps, STARQUERY pre-divides the celestial sphere into multiple equal-area sky regions using the HEALPix algorithm.

HEALPix Division Visualization:

<p align="middle">
  <img src="readme_figs/healpix_list.png" width="400" />
</p>

<p align="middle">
  <img src="readme_figs/healpix_table.png" width="400" />
</p>

The default strategy for dividing the sky area is based on the field of view (FOV):

- For 58.6 > FOV >= 29.3, k = 1, radius of cone search = 37.2; 

- For 29.3 > FOV >= 14.7, k = 2, radius of cone search = 17.0;

- For 14.7 > FOV >= 7.33, k = 3, radius of cone search = 8.3;

- For 7.33 > FOV >= 3.66, k = 4, radius of cone search = 4.1;

- For 3.66 > FOV >= 1.83, k = 5, radius of cone search = 2.1;

```python
>>> # Set the field of view and pixel width
>>> fov,pixel_width = 8,0.01 # In degrees
>>> # Set the maximum number of brightest stars in each sky area
>>> max_num = 30 
>>> # Generate a h5-formatted star catalog indices file
>>> outh5,ratio = hyg_simplified.h5_indices(fov,pixel_width,max_num)
>>> # ratio: Ratio of the solid angles spanned by the square and the cone derived from FOV.
```

A h5-formatted star catalog indices file is generated, recording the center pointing, pixel coordinates of the stars, triangle invariants, and asterism indices of each sky area.

### Read and Parse the h5-formatted Star Catalog Indices File

STARQUERY enables users to read and parse h5-formatted star catalog indices files efficiently. This feature is essential for retrieving detailed star data, including their celestial coordinates, pixel positions, and geometric invariants.

```python
>>> from starcatalogquery import StarCatalog
>>> # Specify the path of the h5-formatted star catalog indices file
>>> infile_h5 = 'starcatalogs/indices/hygv3.7/k3_mag9.0_mcp30_2022.0.h5'
>>> # Read and parse the h5 file
>>> fp_radecs,stars_xy,stars_invariants,stars_asterisms = StarCatalog.read_h5_indices(infile_h5)
```

This command extracts the necessary information from the h5 file.

### Load the Local Offline Star Catalog Database

STARQUERY provides a straightforward method to load local offline star catalog databases. This feature supports loading raw, reduced, and simplified versions of the star catalog.

```python
>>> from starcatalogquery import StarCatalog
>>> # Path of the raw star catalog
>>> dir_from_raw = 'starcatalogs/raw/hygv3.7/res5/'
>>> hyg_raw = StarCatalog.load(dir_from_raw)
>>> # Uncomment the following lines to load the reduced or simplified star catalog
>>> # dir_from_reduced = 'starcatalogs/reduced/hygv3.7/res5/' # Path of the reduced starcatalog
>>> # hyg_reduced = StarCatalog.load(dir_from_reduced)
>>> # dir_from_simplified = 'starcatalogs/simplified/hygv3.7/res5/mag9.0/epoch2022.0/' # Path of the simplified starcatalog
>>> # hyg_simplified = StarCatalog.load(dir_from_simplified)
```

This process ensures quick and efficient access to the star catalog data.

## Change log

- **0.1.14 — Nov 29, 2023**
  - Minor bugs fixed.

- **0.1.13 — Nov 26, 2023**
  
  - Added the capability to download, load, and simplify the *AT-HYGv2.4* star catalog.
  - Implemented a Rotation angle (in radians) in the `pixel_xy` method to align the WCS frame (equivalent to ENU) with the image reference frame.
  - Enhanced the `__repr__` method to provide a formatted summary of the instance.

- **0.1.12 — Sep 23, 2023**
  
  - Minor bugs fixed.

- **0.1.11 — Sep 04, 2023**
  
  - Adjusted parameters for dividing the celestial sphere into multiple equal-area sky regions using the HEALPix algorithm and the corresponding radius of cone search for blind matching of star maps.
  - Reorganized the equal-area sky regions to gradually transition from the celestial equator to the poles.
  - Minor bugs fixed.

- **0.1.10 — Jul 23, 2023**
  
  - Introduced the HEALPix algorithm for dividing the celestial sphere into multiple equal-area sky regions.

- **0.1.8 — Jul 03, 2023**
  
  - Removed dependency on **pyshtools**
    
    - Added function `from_cap` to replace `pysh.SHGrid.from_cap` in *catalog_query.py*
    
    - Using `SphericalCircle` in *astropy* to replace `MakeCircleCoord` in *pyshtools.utils* 
    
    - Added function `from_cap` in *catalog_query.py* to replace `pysh.SHGrid.from_cap`.
    
    - Replaced `MakeCircleCoord` in *pyshtools.utils* with `SphericalCircle` from *Astropy*.

- **0.1.7 — Jun 16, 2023**
  
  - Simplified parameter input for `StarCatalog.load` to facilitate star catalog loading.

- **0.1.5 — May 13, 2023**
  
  - Added the `.invariantfeatures()` method to the `Stars` class, which calculates triangle invariants, constructs a 2D-Tree, and records asterism indices for each triangle.

- **0.1.0 — May 10,  2023**
  
  - Release of the ***starcatalogquery*** package.

## Next Release

- Introduce filtering based on spectral types, allowing for even more refined searches in the star catalog.  

# Contributing

We welcome contributions to the STARQUERY project and are grateful for every bit of help. 

- **Bug Reports**: If you find a bug, please create an issue in our issue tracker. Be sure to include detailed information about the bug and steps to reproduce it.
- **Feature Requests**: If you have ideas for new features or improvements, feel free to open an issue to discuss them.

# Reference

The development of STARQUERY has been influenced and supported by a variety of external resources and tools. Below are some of the key references:

- [**STScI Outerspace**](https://outerspace.stsci.edu/display/GC/WebServices+for+Catalog+Access): The Space Telescope Science Institute provides web services for catalog access, which have been instrumental in the development of STARQUERY.
- [**The Astronomy Nexus**](https://www.astronexus.com/hyg): A valuable resource providing comprehensive star catalog data, such as  the HYG and AT-HYG database. STARQUERY utilizes this resource for accessing detailed astronomical data.
- [**Astroalign**](https://astroalign.quatrope.org/en/latest/index.html): This Python package by Beroiz, M. I. is a significant reference, especially in the development of star pattern recognition features within STARQUERY.
- [**HEALPix**](https://healpix.sourceforge.io): The Hierarchical Equal Area isoLatitude Pixelization (HEALPix) tool has been a reference for implementing the division of the celestial sphere into equal-area sky regions. Learn more about HEALPix at their website.
