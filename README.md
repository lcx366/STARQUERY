# Welcome to the STARQUERY package

[![PyPI version shields.io](https://img.shields.io/pypi/v/starcatalogquery.svg)](https://pypi.python.org/pypi/starcatalogquery/) [![PyPI pyversions](https://img.shields.io/pypi/pyversions/starcatalogquery.svg)](https://pypi.python.org/pypi/starcatalogquery/) [![PyPI status](https://img.shields.io/pypi/status/starcatalogquery.svg)](https://pypi.python.org/pypi/starcatalogquery/) [![GitHub contributors](https://img.shields.io/github/contributors/lcx366/STARQUERY.svg)](https://GitHub.com/lcx366/STARQUERY/graphs/contributors/) [![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/lcx366/STARQUERY/graphs/commit-activity) [![GitHub license](https://img.shields.io/github/license/lcx366/STARQUERY.svg)](https://github.com/lcx366/STARQUERY/blob/master/LICENSE) [![Documentation Status](https://readthedocs.org/projects/starcatalogquery/badge/?version=latest)](http://starcatalogquery.readthedocs.io/?badge=latest) [![Build Status](https://travis-ci.org/lcx366/starcatalogquery.svg?branch=master)](https://travis-ci.org/lcx366/starcatalogquery)

STARQUERY is an archive of scientific routines for establishing an offline star catalog query database.

## Key Features

1. **Offline Star Catalog Database Creation:** Import star catalog data from extensive astronomical sources like [STScI Outerspace](https://outerspace.stsci.edu/display/GC/WebServices+for+Catalog+Access) and [The Astronomy Nexus](https://www.astronexus.com/hyg) to create a comprehensive local offline database.
2. **Data Simplification**: Condense star catalog data, minimizing storage while retaining essential information.
3. **Advanced Search and Query**: Utilize sophisticated search tools, including rectangle and spherical cap (cone) searches, for precise star filtering based on location, magnitude, and other criteria.
4. **Search Area Visualization**: Visualize the scope of the search area and the coverage of corresponding catalog tiles.
5. **Pixel Coordinate Calculation**: Calculate the pixel coordinates of filtered stars given a specific pixel width.
6. **Feature Libraries Construction**: Generate geometrically invariant features based on the spatial configuration of stars to aid in star pattern recognition.
7. **HEALPix-based Sky Area Division**: Divide the celestial sphere into multi-level equal-area sky areas with the HEALPix algorithm.
8. **Local Offline Star Catalog Database Loading**: Load the local offline star catalog database for efficient data access and manipulation.
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

Creating a local offline star catalog with STARQUERY is straightforward and efficient. You can obtain star catalog data from renowned astronomical data repositories like the [STScI Outerspace](https://outerspace.stsci.edu/display/GC/WebServices+for+Catalog+Access) and [The Astronomy Nexus](https://www.astronexus.com/hyg). STARQUERY supports a wide range of star catalogs including 'HYG v3.7', 'AT-HYG v2.4', 'GAIA DR3', 'GSC 30', 'UCAC5', 'NSNOB', '2MASS', and more.

To start building your database, use the `StarCatalog` class in the STARQUERY package. For example, to download the AT-HYG v2.4 star catalog, you can use the following Python script:

```python
>>> from starcatalogquery import StarCatalog
>>> # Get the star catalog AT-HYG v2.4
>>> sc_raw = StarCatalog.get('at-hyg24')
>>> print(sc_raw) # Display basic information about the downloaded catalog
```

### Reduce the Raw Star Catalog

The original star catalog, while comprehensive, often contains extensive information that results in substantial storage requirements and can impede the speed of star queries. To address this, STARQUERY provides a feature to streamline the catalog by extracting only the essential star information. This is achieved using the .reduce() method, which constructs a more compact version of the star catalog.

```python
>>> # Reduce the raw star catalog to its essential components
>>> sc_reduced = sc_raw.reduce()
>>> print(sc_reduced) # Display basic information about the reduced catalog
```

The reduced star catalog will include only the celestial position, proper motion, apparent magnitude, parallax, and epoch of stars, significantly reducing the storage space and improving query performance.

### Simplify the Reduced Star Catalogs

To further enhance its query efficiency, STARQUERY allows users to filter the reduced star catalog based on the limit magnitude of the detector. In addition, it facilitates proper motion corrections to produce a star catalog close to observation time. 

```python
>>> # Set the magnitude threshold and observation epoch for proper motion correction
>>> mag_threshold,t_pm = 12.0,2019.5 
>>> # Simplify the reduced star catalog
>>> sc_simplified = sc_reduced.simplify(mag_threshold,t_pm)
>>> print(sc_simplified) # Display basic information about the simplified catalog 
```

### Query Stars over a Specific Sky Area

STARQUERY allows users to perform a conical or rectangular search on the **raw**, **reduced**, or **simplified** star catalog to find stars within a specific area of the sky. Before conducting a catalog query, it is necessary to establish an index database. Once established, there is no need to repeat the process. **Raw** catalogs and **reduced**  catalogs share the same index database, while **simplified** catalogs require a separate index database to be established. If you do not query the **raw or reduced** catalog, there is no need to establish their index database. Due to the large size of the original catalog files, it is not allowed to directly establish an index database for the original catalog. Instead, a database is established for **reduced** catalogs for sharing, which greatly saves time in establishing an index database.

#### Build a catalog index database for the first time

```python
>>> from starcatalogquery import CatalogDB
>>> sc_reduced.build_indices() # Establish catalog index data table for raw/reduced star catalog
>>> sc_database = CatalogDB(sc_reduced._db_path) # Linking databases
>>> sc_database.add_table(sc_reduced._indices_path) # Add the data table to the database
>>> print(sc_database)
```

Add the index data table for simplified star catalog to the catalog index database, and of course, you can also add index data table from other star catalogs to the database.

```python
>>> sc_simplified.build_indices() # Establish star index data table for simplified star catalog
>>> sc_database.add_table(sc_simplified._indices_path) # Add the data table to the database
>>> print(sc_database)
```

Delete some star index data tables from the the database.

```python
>>> sc_database.del_table(sc_reduced._tb_name)
>>> print(sc_database)
```

#### Perform a star catalog query

The preliminary work has been completed, and now we can proceed with the star catalog query operation. Based on the size of the search area, the star catalog query adaptively select the hierarchy level and corresponding tiles with the HEALPix. The combination of these tiles can fully cover the search area.

##### Conical star query on raw/reduced catalog

Extract all stars located in the search area through cutoff magnitude and boundary conditions.

```python
>>> # Set the center pointing in form of [Ra, Dec] in degrees and search radius in degrees
>>> center,radius = [20,30],15 
>>> # Set the cutoff magnitude and the observation epoch
>>> mag_threshold,t_pm = 12.0,2019.5
>>> # Optionally, set the maximum number of brightest stars to output
>>> max_num = 100
>>> # Perform a cone search on the raw star catalog
>>> sc_raw_stars = sc_raw.search_cone(center,radius,mag_threshold,t_pm,max_num=max_num)
>>> print(sc_raw_stars )
```

By filtering out the brightest stars in each tile, extract the stars located in the search area.

```python
>>> max_num_per_tile = 10
>>> sc_raw_stars = sc_raw.search_cone(center,radius,mag_threshold,t_pm,max_num_per_tile=max_num_per_tile)
```

##### Conical star query on simplified catalog

Star query on simplified catalog no longer requires magnitude truncation and prior proper motion correction. 

```python
>>> # Perform a cone search on the simplified star catalog
>>> sc_simplified_stars = sc_simplified.search_cone(center,radius)
>>> print(sc_simplified_stars)
```

However, astronomical corrections can be applied, including 

    - 'proper-motion': Corrects for the motion of stars across the sky due to their velocities.
    - 'aberration': Corrects for the apparent shift in star positions due to the motion of the observer.
    - 'parallax': Corrects for the apparent shift in star positions due to the change in observer's viewpoint as the Earth orbits the Sun.
    - 'deflection': Corrects for the bending of light from stars due to the gravitational field of the Sun, based on general relativity.

```python
>>> astrometry_corrections = {'t':'2019-02-26T20:11:14.347','proper-motion':None,'aberration':(0.55952273, -1.17780654,  7.50324956),'parallax':None}
>>> sc_simplified_stars = sc_simplified.search_cone(center,radius,astrometry_corrections=astrometry_corrections)
```

In *astrometry_corrections* dict, 
    - 't' -> [str] Observation time in UTC, which specifies the time at which corrections are to be applied.
    - 'proper-motion' -> [None] If present, apply proper motion correction.
    - 'aberration' -> [tuple] Observer's velocity relative to Earth's center, (vx, vy, vz) in km/s.
    - 'parallax' -> [None] If present, apply parallax correction.
    - 'deflection' -> [None] If present, apply light deflection correction.

##### Rectangular star query on raw/reduced/Simplified catalog

STARQUERY also provides the capability to perform a rectangle search in the raw, reduced, and simplified star catalogs. Its usage is similar to that of cone search.

```python
>>> # Set a rectangular search area in form of [ra_min, dec_min, ra_max, dec_max] in degrees
>>> radec_box = [5,15,35,45]
>>> # Set the cutoff magnitude and the observation epoch
>>> mag_threshold,t_pm = 12.0,2019.5
>>> # Optionally, set the maximum number of brightest stars to output
>>> max_num = 100
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
>>> # Specify the pixel width in degrees
>>> pixel_width = 0.01
>>> # Calculate the pixel coordinates using the TAN projection
>>> sc_simplified_stars.pixel_xy(pixel_width)
>>> # Retrieve the calculated pixel coordinates
>>> xy = sc_simplified_stars.xy
>>> print(xy)
```

### Calculate the Geometric Invariants and the Asterism Indices of the Filtered Stars

The calculation of triangle invariants draws inspiration from the Astroalign software package developed by Beroiz, M. I. (https://astroalign.quatrope.org/en/latest/). The calculation of quads invariants refers to Astrometry.net(https://astrometry.net).

Steps for the calculation of triangle invariants:

1. **Derive Geometric Invariants:** For each possible triangle formed by groups of three stars, calculates unique geometric invariants (ratios L2/L1 and L1/L0), where L2, L1, and L0 are the sides of each triangle sorted in descending order.
2. **Construct a 2D-Tree Structure:** Utilizing these invariants, a 2D-Tree structure is constructed to facilitate efficient spatial queries.
3. **Associate Invariants with Star Indices:** Each invariant set is linked to the indices of the stars forming the corresponding triangle.

```python
>>> sc_simplified_stars.invariantfeatures(mode_invariants='triangles')
>>> # Retrieve the generated invariants, asterisms, and 2D-Tree structure
>>> invariants = sc_simplified_stars.invariants
>>> asterisms = sc_simplified_stars.asterisms
>>> kdtree = sc_simplified_stars.kdtree
>>> print(invariants,asterisms,kdtree)
```

Steps for the calculation of quad invariants:

1. **Derive Geometric Invariants:** For each possible quad formed by groups of four stars, use the most widely separated pair of stars in the quad to define a local coordinate system, labeling them "A" and "B." The remaining two stars, "C" and "D," have their positions (xC, yC) and (xD, yD) respectively in this coordinate system. The geometric hash code is the 4-vector (xC, yC, xD, yD).
2. **Break symmetries** By requiring that xC ≤ xD and xC + xD ≤ 1, only consider the permutation that meets these conditions.
3. **Construct a 2D-Tree Structure:** Utilizing these invariants, a 4D-Tree structure is constructed to facilitate efficient spatial queries.
4. **Associate Invariants with Star Indices:** Each invariant set is linked to the indices of the stars forming the corresponding quad.

```python
>>> sc_simplified_stars.invariantfeatures(mode_invariants='quads') 
>>> # Retrieve the generated invariants, asterisms, and 4D-Tree structure
>>> invariants = sc_simplified_stars.invariants
>>> asterisms = sc_simplified_stars.asterisms
>>> kdtree = sc_simplified_stars.kdtree
>>> print(invariants,asterisms,kdtree)
```

### Visualization

STARQUERY provides a visualization feature to display the scope of the search area and the coverage of the corresponding catalog tiles. 

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

To facilitate blind matching of star maps, STARQUERY divides the celestial sphere into multi-level equal-area sky regions using the HEALPix algorithm in advance. The default strategy for dividing the celestial sphere ensures that each tile's size is greater than one-quarter of the FOV and less than or equal to half of the FOV. For each level of tiles in the star catalog (from K1 to K11), extract the brightest stars (such as top 5) and integrate them into the K1 level. Then, calculate the geometric invariants and store them into a h5-formatted hash file.

<p align="middle">
  <img src="readme_figs/healpix_list.png" width="400" />
</p>

<p align="middle">
  <img src="readme_figs/healpix_table.png" width="400" />
</p>

```python
>>> k_min,k_max = 1,6
>>> sc_simplified.h5_hashes(k_min,k_max,mode_invariants='triangles') # mode_invariants='quads'
>>> print(sc_simplified.hashed_h5)
```

A h5-formatted star catalog hashed file is generated, recording the center pointing, pixel coordinates of the stars, geometric invariants, and asterism indices for multi-level sky region.

### Read the h5-formatted Star Catalog Hashed File

```python
>>> infile_h5 = 'starcatalogs/indices/at-hyg24_mag12.0_epoch2019.5_triangles_K1_K6.h5'
>>> simplified_catalog.read_h5_hashes(infile_h5)
>>> print(sc_simplified.hashed_data)
```

### Load the Local Offline Star Catalog Database

STARQUERY provides a straightforward method to load local offline star catalog databases, ensuring quick access to the star catalog data.

```python
>>> from starcatalogquery import StarCatalog
>>> dir_from_raw = 'starcatalogs/raw/at-hyg24/' # Path of the raw star catalog
>>> sc_raw = StarCatalog.load(dir_from_raw)
>>> dir_from_reduced = 'starcatalogs/reduced/at-hyg24/' # Path of the reduced starcatalog
>>> sc_reduced = StarCatalog.load(dir_from_reduced)
>>> dir_from_simplified = 'starcatalogs/simplified/at-hyg24/mag12.0/epoch2019.5/' # Path of the simplified starcatalog
>>> sc_simplified = StarCatalog.load(dir_from_simplified)
```

## Change log

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

## Next Release

- Enhance the efficiency of nearest neighbor search to provide support for improving the speed of blind matching of star maps.

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
