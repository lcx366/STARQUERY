# Welcome to the STARQUERY package

[![PyPI version shields.io](https://img.shields.io/pypi/v/starcatalogquery.svg)](https://pypi.python.org/pypi/starcatalogquery/) [![PyPI pyversions](https://img.shields.io/pypi/pyversions/starcatalogquery.svg)](https://pypi.python.org/pypi/starcatalogquery/) [![PyPI status](https://img.shields.io/pypi/status/starcatalogquery.svg)](https://pypi.python.org/pypi/starcatalogquery/) [![GitHub contributors](https://img.shields.io/github/contributors/lcx366/STARQUERY.svg)](https://GitHub.com/lcx366/STARQUERY/graphs/contributors/) [![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/lcx366/STARQUERY/graphs/commit-activity) [![GitHub license](https://img.shields.io/github/license/lcx366/STARQUERY.svg)](https://github.com/lcx366/STARQUERY/blob/master/LICENSE) [![Documentation Status](https://readthedocs.org/projects/starcatalogquery/badge/?version=latest)](http://starcatalogquery.readthedocs.io/?badge=latest) [![Build Status](https://travis-ci.org/lcx366/starcatalogquery.svg?branch=master)](https://travis-ci.org/lcx366/starcatalogquery)

STARQUERY is an archive of scientific routines for establishing an offline star catalog query database.

## Key Features

1. **Offline Star Catalog Database Creation:** Import star catalog data from astronomical sources like [STScI Outerspace](https://outerspace.stsci.edu/display/GC/WebServices+for+Catalog+Access) and [The Astronomy Nexus](https://www.astronexus.com/hyg) to create a comprehensive local offline database.
2. **Data Simplification**: Extract necessary star information from huge star catalog files.
3. **Star Query**: Utilize the rectangle and spherical cap (cone) query for star filtering.
4. **Visualization**: Visualize the search area and the corresponding catalog tiles.
5. **Pixel Coordinate Calculation**: Calculate the pixel coordinates of filtered stars given a specific pixel width.
6. **Invariant Features Computation**: Compute geometrically invariant features based on the spatial configuration of stars to aid in star pattern recognition.
7. **HEALPix-based Sky Area Division**: Divide the celestial sphere into multi-level equal-area sky areas with the HEALPix algorithm.
8. **Star Catalog Database Loading**: Load the local offline star catalog database for efficient data access and manipulation.
9. **Astronomical Corrections**: Apply corrections for proper motion, aberration, parallax, and deflection to star positions.

## How to Install

Install STARQUERY using the pip command in your terminal:

```
pip install starcatalogquery
pip install starcatalogquery --upgrade # to upgrade a pre-existing installation
```

If an error message similar to "`ERROR: Could not build wheels for cartopy, which is required to install pyproject.toml-based projects`" is displayed, a good solution is

```bash
mamba install h5py
mamba install cartopy 
```

## How to use

### Build an Offline Star Catalog Database

Creating a local offline star catalog database with STARQUERY is straightforward and efficient. The star catalog data can be obtained from the renowned astronomical data repositories like [STScI Outerspace](https://outerspace.stsci.edu/display/GC/WebServices+for+Catalog+Access) and [The Astronomy Nexus](https://www.astronexus.com/hyg). STARQUERY supports a wide range of star catalogs, including 'HYG v3.7', 'AT-HYG v2.4', 'GAIA DR3', 'GSC 30', 'UCAC5', 'USNOB1.0', '2MASS', and more.

To start building your database, for example, just download the AT-HYG v2.4 star catalog.

```python
>>> from starcatalogquery import StarCatalog
>>> sc_raw = StarCatalog.get('at-hyg24') # Get the raw star catalog AT-HYG v2.4
>>> print(sc_raw) # Display basic information about the downloaded catalog
```

### Simplify the Raw Star Catalog

The raw star catalog usually contains extensive information about stars that result in substantial storage requirements, thus impeding the speed of star queries. To address this, by extracting only the essential star information, a more compact version of the star catalog is created.

```python
>>> sc_reduced = sc_raw.reduce()
>>> print(sc_reduced)
```

The reduced star catalog will only contains: celestial position, proper motion, apparent magnitude, parallax, and epoch of stars.

### Filter the Reduced Star Catalogs

For further enhancing the query efficiency, the reduced star catalog is filtered based on the limit magnitude of the detector. Simultaneously, the proper motion corrections is applied to make the star catalog close to the observation time. 

```python
>>> mag_threshold = 13.0 # Set the magnitude threshold
>>> t_pm = 2019.5 # Set the observation time for proper motion correction
>>> sc_simplified = sc_reduced.simplify(mag_threshold,t_pm)
>>> print(sc_simplified) 
```

### Query Stars over a Specific Sky Area

STARQUERY impliments a conical and rectangular query on the **raw**, **reduced**, or **simplified** star catalog. Before conducting a catalog query, an index database needs to be prepared. The **raw** catalog is too large to directly build an index database, and it shares the same index database with the **reduced**  catalog. If you do not query the **raw or reduced** catalog, there is no need to build their index database.

#### Build a catalog index database for the first time

```python
>>> from starcatalogquery import CatalogDB
>>> sc_database = CatalogDB(sc_reduced._db_path) # Linking databases
>>> sc_reduced.build_indices() # Establish catalog index data table
>>> sc_database.add_table(sc_reduced._indices_path) # Add the data table to the database
>>> print(sc_database)
```

Add the index data table for the simplified star catalog to the database.

```python
>>> sc_simplified.build_indices() # Establish catalog index data table for simplified star catalog
>>> sc_database.add_table(sc_simplified._indices_path) # Add the data table to the database
>>> print(sc_database)
```

Delete the star index data tables from the the database.

```python
>>> sc_database.del_table(sc_reduced._tb_name)
>>> print(sc_database)
```

#### Perform a star catalog query

Based on the size of the search area, the star catalog query adaptively select the hierarchy level and corresponding tiles with the HEALPix.

##### Conical star query on raw/reduced catalog

Extract all stars located in the search area through cutoff magnitude and search boundaries.

```python
>>> center = [20,30] # Set the center pointing in form of [Ra, Dec] in degrees
>>> radius = 5 # Set the search radius in degrees
>>> mag_threshold = 13.0 # Set the cutoff magnitude
>>> t_pm = 2019.5 # Set the observation time
>>> max_num = 100 # Optionally, set the maximum number of brightest stars to output
>>> sc_raw_stars = sc_raw.search_cone(center,radius,mag_threshold,t_pm,max_num=max_num) # Cone search on the raw star catalog
>>> print(sc_raw_stars)
```

Alternatively, by filtering out the brightest stars in each tile, extract the stars located in the search area.

```python
>>> max_num_per_tile = 5
>>> sc_raw_stars = sc_raw.search_cone(center,radius,mag_threshold,t_pm,max_num_per_tile=max_num_per_tile)
```

##### Conical star query on simplified catalog

Star query on simplified catalog no longer requires magnitude truncation and prior proper motion correction. 

```python
>>> sc_simplified_stars = sc_simplified.search_cone(center,radius)
>>> print(sc_simplified_stars)
```

During the star query, astronomical corrections can be applied:

    - 'proper-motion': Corrections for the motion of stars across the sky due to their velocities.
    - 'aberration': Corrections for the apparent shift in star positions due to the motion of the observer.
    - 'parallax': Corrections for the apparent shift in star positions due to the change in observer's viewpoint as the Earth orbits the Sun.
    - 'deflection': Corrections for the bending of light from stars due to the gravitational field of the Sun, based on general relativity.

```python
>>> astrometry_corrections = {'t':'2019-02-26T20:11:14.347','proper-motion':None,'aberration':(0.55952273, -1.17780654,  7.50324956),'parallax':None}
>>> sc_simplified_stars = sc_simplified.search_cone(center,radius,astrometry_corrections=astrometry_corrections)
```

Details on the *astrometry_corrections* dict:
    - 't': Observation time in UTC, which specifies the time at which corrections are to be applied.
    - 'proper-motion': If present, apply proper motion correction.
    - 'aberration': Observer's velocity relative to Earth's center, (vx, vy, vz) in km/s.
    - 'parallax': If present, apply parallax correction.
    - 'deflection': If present, apply light deflection correction.

##### Rectangular star query on raw/reduced/Simplified catalog

The usage is similar to that of the cone search.

```python
>>> radec_box = [5,15,15,25] # Set the rectangular search area in form of [ra_min, dec_min, ra_max, dec_max] in degrees
>>> mag_threshold = 13.0 # Set the cutoff magnitude
>>> t_pm = 2019.5 # Set the observation time
>>> max_num = 100 # Optionally, set the maximum number of brightest stars to output
>>> # Perform a rectangle search on the raw star catalog
>>> sc_raw_stars = sc_raw.search_box(radec_box,mag_threshold,t_pm,max_num=max_num)
>>> # Perform a rectangle search on the reduced star catalog
>>> sc_reduced_stars = sc_reduced.search_box(radec_box,mag_threshold,t_pm,max_num=max_num)
>>> # Perform a rectangle search on the simplified star catalog
>>> sc_simplified_stars = sc_simplified.search_box(radec_box)
```

### Calculate the Pixel Coordinates of the Filtered Stars

Calculate the pixel coordinates of filtered stars, using the *TANGENT* projection in World Coordinate System (WCS) transformations. 

```python
>>> pixel_width = 0.001 # Set the pixel width in degrees
>>> # Calculate the pixel coordinates using the TAN projection
>>> sc_simplified_stars.pixel_xy(pixel_width)
>>> # Retrieve the calculated pixel coordinates
>>> xy = sc_simplified_stars.xy
>>> print(xy)
```

### Calculate the Geometric Invariants and the Asterism Indices of the Filtered Stars

- Method for calculating the triangle invariants draws inspiration from the Astroalign software package developed by Beroiz, M. I. (https://astroalign.quatrope.org/en/latest/). 
- Method for calculating the quads invariants refers to Astrometry.net(https://astrometry.net).

Steps for the calculation of triangle invariants:

1. **Derive Geometric Invariants:** For each possible triangle formed by groups of three stars, calculates unique geometric invariants (ratios L2/L1 and L1/L0), where L2, L1, and L0 are the sides of each triangle sorted in descending order.
2. **Construct a 2D-Tree Structure:** Utilizing these invariants, a 2D-Tree structure is constructed to facilitate efficient spatial queries.
3. **Associate Invariants with Star Indices:** Each invariant set is linked to the indices of the stars forming the corresponding triangle.

Steps for the calculation of quad invariants:

1. **Derive Geometric Invariants:** For each possible quad formed by groups of four stars, use the most widely separated pair of stars in the quad to define a local coordinate system, labeling them "A" and "B." The remaining two stars, "C" and "D," have their positions (xC, yC) and (xD, yD) respectively in this coordinate system. The geometric hash code is the 4-vector (xC, yC, xD, yD).
2. **Break symmetries** By requiring that xC ≤ xD and xC + xD ≤ 1, only consider the permutation that meets these conditions.
3. **Construct a 2D-Tree Structure:** Utilizing these invariants, a 4D-Tree structure is constructed to facilitate efficient spatial queries.
4. **Associate Invariants with Star Indices:** Each invariant set is linked to the indices of the stars forming the corresponding quad.

```python
>>> mode_invariants = 'triangles' # quads
>>> sc_simplified_stars.invariantfeatures(mode_invariants)
>>> # Retrieve the generated invariants, asterisms, and 2D-Tree structure
>>> invariants = sc_simplified_stars.invariants
>>> asterisms = sc_simplified_stars.asterisms
>>> kdtree = sc_simplified_stars.kdtree
>>> print(invariants,asterisms,kdtree)
```

### Visualization

Display the scope of the search area and the coverage of the corresponding catalog tiles. 

```python
>>> sc_simplified_stars.tiles_draw() # Visualize the cone search area
```

<p align="middle">
  <img src="readme_figs/box.png" width="500" />
</p>

```python
>>> sc_simplified_stars.tiles_draw() # Visualize the box search area
```

<p align="middle">
  <img src="readme_figs/cone.png" width="500" />
</p>

### Divide the Sky into multi-level Equal-Area Sky Regions

To facilitate blind matching of star maps, STARQUERY divides the celestial sphere into multi-level equal-area sky regions using the HEALPix algorithm in advance. For each level(from K1 to K11), the brightest stars (such as top 5) are extracted and integrated into the K1 level. The geometric invariants are then calculated and stored into a h5-formatted hash file.

<p align="middle">
  <img src="readme_figs/healpix_list.png" width="400" />
</p>

<p align="middle">
  <img src="readme_figs/healpix_table.png" width="400" />
</p>

```python
>>> k_min,k_max = 1,6
>>> sc_simplified.h5_hashes(k_min,k_max,mode_invariants)
>>> print(sc_simplified.hashed_h5)
```

The h5-formatted hash data records the center pointing, pixel coordinates of the stars, geometric invariants, and asterism indices for multi-level sky region.

### Read the h5-formatted Hash File

```python
>>> infile_h5 = 'starcatalogs/indices/at-hyg24_mag13.0_epoch2019.5_triangles_K1_K6.h5'
>>> sc_hashed = sc_simplified.read_h5_hashes(infile_h5)
>>> print(sc_hashed.hashed_data)
```

### Load the Local Offline Star Catalog Database

```python
>>> from starcatalogquery import StarCatalog
>>> dir_from_raw = 'starcatalogs/raw/at-hyg24/' # Path of the raw star catalog
>>> sc_raw = StarCatalog.load(dir_from_raw)
>>> dir_from_reduced = 'starcatalogs/reduced/at-hyg24/' # Path of the reduced starcatalog
>>> sc_reduced = StarCatalog.load(dir_from_reduced)
>>> dir_from_simplified = 'starcatalogs/simplified/at-hyg24/mag13.0/epoch2019.5/' # Path of the simplified starcatalog
>>> sc_simplified = StarCatalog.load(dir_from_simplified)
```

## Change log

- **1.1.3 — Oct 14, 2024**
  
  - Fixed an error in parallel processing when applying astronomical corrections (primarily caused by JPL ephemeris memory mapping).

- **1.1.2 — Sep 29, 2024**
  
  - Normalize RA to [0, 360) when converting cartesian coordinates to spherical coordinates.

- **1.1.1 — Sep 04, 2024**
  
  - The number of nearest neighbor stars for constructing geometric invariants has been increased to *nine*.
  - Added functions `vectorized_unique_quads` and `vectorized_unique_triangles` for vectorizing the calculation of geometric invariants.
  - The pixel scale used for conversion between pixel coordinates and astronomical coordinates via WCS (World Coordinate System) transformations has been refined from *0.01 degrees/pixel* to *0.001 degrees/pixel*.
  - The number of stars extracted per tile at all levels has been reduced to *five* during the generation of hash files for geometric invariants.
  - Adjusted the calculation for determining healpix level based on field of view size.
  - A new class, `H5HashesData`, has been introduced.

- **1.0.5 — Aug 08, 2024**
  
  - Fixed the memory overflow issue caused by processing giant star catalog tile files.

- **1.0.3 — Jul 27, 2024**
  
  - Raise the healpix level of the star catalog tiles from K4 to K5 to avoid the trouble of downloading large files from remote servers that are prone to failure.

- **1.0.2 — Jul 17, 2024**
  
  - In checking the validity of star catalog files, use the `wc -l` command with the subprocess module to quickly count the number of lines in the CSV file, which is more efficient than a pure Python implementation.
  - Fixed a bug in the validity check of star catalog files

- **1.0.1 — Jul 16, 2024**
  
  - Added parameters that limit the magnitude range to avoid the problem of remote server data overflow and download failure.

- **1.0.0 — Jul 04, 2024**
  
  - Replaced the spherical rectangular partitioning with a multi-level equal-area partitioning strategy based on the HEALPix algorithm.
  - Established a star catalog index database system that adaptively selects the appropriate level based on the size of the search area or FOV.
  - Added functionality to generate geometrically invariant features using four stars and create hash tables.
  - Enhanced star position corrections with proper motion, annual parallax, aberration, and light deflection adjustments.
  - Added visualization capabilities for multi-level equal-area sky region partitioning.

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
