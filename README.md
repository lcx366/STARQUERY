# Welcome to the STARQUERY package

[![PyPI version shields.io](https://img.shields.io/pypi/v/starcatalogquery.svg)](https://pypi.python.org/pypi/starcatalogquery/) [![PyPI pyversions](https://img.shields.io/pypi/pyversions/starcatalogquery.svg)](https://pypi.python.org/pypi/starcatalogquery/) [![PyPI status](https://img.shields.io/pypi/status/starcatalogquery.svg)](https://pypi.python.org/pypi/starcatalogquery/) [![GitHub contributors](https://img.shields.io/github/contributors/lcx366/STARQUERY.svg)](https://GitHub.com/lcx366/STARQUERY/graphs/contributors/) [![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/lcx366/STARQUERY/graphs/commit-activity) [![GitHub license](https://img.shields.io/github/license/lcx366/STARQUERY.svg)](https://github.com/lcx366/STARQUERY/blob/master/LICENSE) [![Documentation Status](https://readthedocs.org/projects/starcatalogquery/badge/?version=latest)](http://starcatalogquery.readthedocs.io/?badge=latest) [![Build Status](https://travis-ci.org/lcx366/starcatalogquery.svg?branch=master)](https://travis-ci.org/lcx366/starcatalogquery)

**STARQUERY** is a powerful Python package designed for astronomers and space researchers who need efficient and offline access to star catalog data. By leveraging offline databases, STARQUERY allows users to perform high-speed queries and complex filtering on large astronomical datasets. The package is ideal for tasks such as star map matching, celestial navigation, and space object tracking.

## 🚀 Key Features

1. **Offline Star Catalog Database Creation**
   
   - Easily import star catalog data from popular astronomical sources such as the [STScI Outerspace](https://outerspace.stsci.edu/display/MASTDATA/Catalog+Access) and [The Astronomy Nexus](https://www.astronexus.com/hyg).
   
   - Build a local star catalog database for offline usage, eliminating the need for repeated online queries and improving performance for data-intensive tasks.

2. **Data Simplification**
   
   - Extract essential star attributes (like position, magnitude, etc.) from massive star catalog files to streamline data storage and access, speeding up the query process.

3. **Star Query**
   
   - **Rectangle Query**: Search for stars within a defined rectangular area on the sky.
   
   - **Spherical Cap (Cone) Query**: Search for stars within a specific angular distance from a given point.

4. **Visualization**
   
   - Visualize the search area and catalog tiles, helping users verify that the search area is properly covered by the catalog tiles.

5. **Pixel Coordinate Calculation**
   
   - Convert celestial coordinates (RA, Dec) to pixel coordinates with a specified pixel width, useful for aligning star charts with sensor data, such as from telescopes or cameras.

6. **Invariant Features Computation**
   
   - Calculate geometrically invariant features based on the spatial configuration of stars, enabling robust star matching across different views.

7. **HEALPix-based Sky Area Division**
   
   - Utilize the HEALPix algorithm to divide the celestial sphere into equal-area tiles at multiple levels of granularity, enabling efficient data partitioning and indexing, speeding up spatial queries.

8. **Astronomical Corrections**
   
   - Enhances the accuracy of star positions by applying
     
     - **Proper Motion**: Adjust star positions based on their velocity across the sky.
     
     - **Aberration**: Correct for the apparent shift in star positions due to Observer’s motion.
     
     - **Parallax**: Adjust for the apparent positional shift due to Earth’s orbit around the Sun.
     
     - **Deflection**: Correct for the bending of light caused by gravitational fields of the Sun.

## 🛠️ How to Install

To install STARQUERY, simply use `pip` in your terminal:

```
pip install starcatalogquery
pip install starcatalogquery --upgrade # to upgrade a pre-existing installation
```

If an error message similar to "`ERROR: Could not build wheels for cartopy, which is required to install pyproject.toml-based projects`" is displayed, a good solution is

```bash
mamba install h5py
mamba install cartopy 
```

## **📚** How to Use

Below are some basic examples to help you get started with STARQUERY.

### Build an Offline Star Catalog Database

To start building your database, for example, just download the AT-HYG v2.4 star catalog.

```python
>>> from starcatalogquery import StarCatalog
>>> sc_raw = StarCatalog.get('at-hyg24') # Get the raw star catalog AT-HYG v2.4
>>> print(sc_raw)
```

The `StarCatalog.get()` method fetches the specified star catalog (**at-hyg24** in this example) from the online source, mainly [STScI Outerspace](https://outerspace.stsci.edu/display/GC/WebServices+for+Catalog+Access) and [The Astronomy Nexus](https://www.astronexus.com/hyg). By default, the downloaded catalog is saved to the **current working directory** under the path: `starcatalogs/raw/at-hyg24/`. The catalog is divided into **K5-level tiles**, following the HEALPix hierarchical structure. This means that the celestial sphere is divided into $4^5 * 12 = 12,288$ tiles, each stored as an individual file. STARQUERY supports a wide range of star catalogs, which are listed below with their corresponding identifiers:

</center>

| Star Catalog Name     | Identifier |
|:---------------------:|:----------:|
| HYG v3.7              | hyg37      |
| AT-HYG v2.4           | at-hyg24   |
| GAIA DR3              | gaiadr3    |
| Guide Star Catalog 30 | gsc30      |
| UCAC5                 | ucac5      |
| USNO-B1.0             | usnob      |
| 2MASS                 | 2mass      |

### Simplify the Raw Star Catalog

The raw star catalog typically contains extensive information about stars, resulting in large file sizes that can slow down query performance. To optimize this, we can extract only the essential information, creating a more compact version of the star catalog.

```python
>>> sc_reduced = sc_raw.reduce()
>>> print(sc_reduced)
```

The `reduce()` method extracts only the essential star attributes: 

- **Celestial Position** (RA, Dec) in degrees

- **Proper Motion** in milliarcseconds per year

- **Apparent Magnitude**

- **Distance** in kiloparsecs

- **Epoch**

The reduced catalog is saved to the current working directory under the path: `starcatalogs/reduced/at-hyg24/`.

### Filter the Reduced Star Catalogs

To further improve query efficiency, the reduced star catalog can be filtered based on the detector’s magnitude limit. Additionally, proper motion corrections are applied to adjust star positions to bring them closer to the specified observation time.

```python
>>> mag_threshold = 13.0 # Set the magnitude threshold
>>> t_pm = 2019.5 # Set the observation time for proper motion correction
>>> sc_simplified = sc_reduced.simplify(mag_threshold,t_pm)
>>> print(sc_simplified) 
```

The `simplify()` method filters stars based on a specified magnitude threshold (e.g., mag_threshold = 13.0) and updates star positions according to their proper motion, adjusted to the specified observation time (t_pm = 2019.5). The simplified catalog is saved to the current working directory under the path: `starcatalogs/simplified/at-hyg24/mag13.0/epoch2019.5/`.

### Load the Local Offline Star Catalog Database

Once you’ve downloaded and processed your star catalog, you can easily load the local offline database for further analysis. Here’s how to load the raw, reduced, and simplified star catalogs:

```python
>>> from starcatalogquery import StarCatalog

>>> # Load the raw star catalog
>>> dir_from_raw = 'starcatalogs/raw/at-hyg24/'  # Path to the raw star catalog
>>> sc_raw = StarCatalog.load(dir_from_raw)
>>> print(sc_raw)

>>> # Load the reduced star catalog
>>> dir_from_reduced = 'starcatalogs/reduced/at-hyg24/'  # Path to the reduced star catalog
>>> sc_reduced = StarCatalog.load(dir_from_reduced)
>>> print(sc_reduced)

>>> # Load the simplified star catalog
>>> dir_from_simplified = 'starcatalogs/simplified/at-hyg24/mag13.0/epoch2019.5/'  # Path to the simplified star catalog
>>> sc_simplified = StarCatalog.load(dir_from_simplified)
>>> print(sc_simplified)
```

### Query Stars over a Specific Sky Area

STARQUERY supports both **conical** and **rectangular** queries on the **raw**, **reduced**, or **simplified** star catalogs. Before performing any queries, an index database needs to be prepared to optimize search efficiency.

```python
>>> from starcatalogquery import CatalogDB
>>> sc_database = CatalogDB(sc_simplified._db_path) # Linking databases
>>> sc_simplified.build_indices() # Establish catalog index data table
>>> sc_database.add_table(sc_simplified._indices_path) # Add the data table to the database
>>> print(sc_database)
```

By default, all index files and databases are stored in the current working directory under the path `starcatalogs/indices/`. The index files are organized as CSV files containing hierarchical levels from **K1** to **K11**, along with the **K5_SUB** row number that identifies the position of stars in the **K5****-level files:

| K1  | ... | K5  | ... | K11 | K5_SUB |
|:---:|:---:|:---:|:---:|:---:|:------:|
| 0   | ... | 0   | ... | 267 | 48     |
| 0   | ... | 0   | ... | 366 | 76     |
| 0   | ... | 0   | ... | 421 | 28     |
| 0   | ... | 0   | ... | 846 | 92     |
| 0   | ... | 0   | ... | 492 | 16     |
| 0   | ... | 0   | ... | 775 | 47     |
| ... | ... | ... | ... | ... | ...    |

💡If you need to remove any index tables from the database:

```python
>>> sc_database.del_table(sc_simplified._tb_name)
>>> print(sc_database)
```

#### 🔍 Perform a Star Catalog Query

STARQUERY allows efficient querying of stars over a specific sky area. Depending on the size of the search region, the query system adaptively selects the appropriate hierarchical level and tiles using the **HEALPix** scheme. The package supports both **conical** (circular) and **rectangular** queries on the **raw**, **reduced**, and **simplified** catalogs.

##### Conical Star Query

You can extract all stars within a specified search area by setting a cutoff magnitude and other search parameters.

```python
from starcatalogquery import StarCatalog

# Set search parameters
center = [20, 30]  # Center point [RA, Dec] in degrees
radius = 5  # Search radius in degrees
mag_threshold = 13.0  # Cutoff magnitude
t_pm = 2019.5  # Observation time
max_num = 100  # Optional: Maximum number of stars to return

# Perform a conical search
sc_raw_stars = sc_raw.search_cone(center, radius, mag_threshold, t_pm, max_num=max_num)
print(sc_raw_stars)
```

When performing a star catalog query, especially over a large area, it’s possible that the brightest stars might cluster in specific regions within the search area. This can lead to an uneven distribution of stars in your results, which may not be ideal for applications like star pattern recognition or celestial navigation. To address this, STARQUERY allows you to **limit the number of stars extracted per HEALPix tile**. This ensures a more even selection of bright stars across the entire search region, preventing a concentration of stars in just one part of the sky.

```python
max_num_per_tile = 5
sc_raw_stars = sc_raw.search_cone(center, radius, mag_threshold, t_pm, max_num_per_tile=max_num_per_tile)
```

Queries on the simplified catalog do not require magnitude truncation or prior proper motion correction.

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
>>> # Define astrometry corrections
>>> astrometry_corrections = {
    't': '2019-02-26T20:11:14.347',  # Observation time (UTC)
    'proper-motion': None,  # Apply proper motion correction if specified
    'aberration': (0.5595, -1.1778, 7.5032),  # Observer's velocity (vx, vy, vz) in km/s
    'parallax': None,  # Apply parallax correction if specified
    'deflection': None  # Apply light deflection correction if specified
>>> }

>>> # Perform a conical query with corrections
>>> sc_simplified_stars = sc_simplified.search_cone(center, radius, astrometry_corrections=astrometry_corrections)
>>> print(sc_simplified_stars)
```

##### Rectangular Star Query

The rectangular query works similarly to the conical search but defines a rectangular region instead.

```python
>>> radec_box = [5, 15, 15, 25]  # [ra_min, dec_min, ra_max, dec_max] in degrees
>>> mag_threshold = 13.0  # Cutoff magnitude
>>> t_pm = 2019.5  # Observation time
>>> max_num = 100  # Optional: Maximum number of stars to return

>>> # Rectangle search on the raw catalog
>>> sc_raw_stars = sc_raw.search_box(radec_box, mag_threshold, t_pm, max_num=max_num)

>>> # Rectangle search on the reduced catalog
>>> sc_reduced_stars = sc_reduced.search_box(radec_box, mag_threshold, t_pm, max_num=max_num)

>>> # Rectangle search on the simplified catalog (no magnitude limit needed)
>>> sc_simplified_stars = sc_simplified.search_box(radec_box, max_num=max_num)
>>> print(sc_simplified_stars)
```

### Calculate the Pixel Coordinates of the Filtered Stars

After filtering the stars, you can convert their celestial coordinates (RA, Dec) into pixel coordinates using the **TANGENT (TAN) projection** in the World Coordinate System (WCS). This transformation is particularly useful for aligning star catalogs with image sensors or star trackers.

```python
>>> pixel_width = 0.001 # Set the pixel width in degrees
>>> # Calculate the pixel coordinates using the TAN projection
>>> sc_simplified_stars.pixel_xy(pixel_width)
>>> # Retrieve the calculated pixel coordinates
>>> xy = sc_simplified_stars.xy
>>> print(xy)
```

### Calculate the Geometric Invariants and the Asterism Indices of the Filtered Stars

To enhance star pattern recognition, STARQUERY can compute geometric invariants for triangles or quads formed by groups of stars. 

**1. Triangle Invariants**

More details refer to the Astroalign package developed by Beroiz, M. I. ([Astroalign Documentation](https://astroalign.quatrope.org/en/latest/)).

🌟 **Steps**:

1. **Derive Geometric Invariants**
   
   - For each possible triangle formed by groups of three stars, compute unique geometric invariants as ratios $\left(\frac{L_2}{L_1},\frac{L_1}{L_0}\right)$, where,  $L_2$, $L_1$, $L_0$ are the sides of the triangle sorted in descending order.

2. **Construct a 2D-Tree Structure**
   
   - Use these invariants to build a 2D-Tree for efficient spatial queries.

3. **Associate Invariants with Star Indices**
   
   - Link each invariant set to the indices of the stars that form the corresponding triangle.

**2. Quad Invariants**

More details refer to Astrometry.net ([Astrometry.net Documentation](https://astrometry.net)).

🌟 **Steps**:

1. **Derive Geometric Invariants**
   
   - For each group of four stars, choose the most widely separated pair to define a local coordinate system, labeling them as “A” and “B”.
   
   - The remaining stars “C” and “D” are positioned relative to this coordinate system with coordinates of ($x_C$, $y_C$) and ($x_D$, $y_D$) .
   
   - The resulting **geometric hash code** is the 4-vector:  ($x_C, y_C, x_D, y_D$) .

2. **Break Symmetries**
   
   - Enforce constraints  $x_C \leq x_D$ and  $x_C + x_D \leq 1$ to reduce redundant permutations.

3. **Construct a 4D-Tree Structure**
   
   - Utilize the invariants to build a 4D-Tree for efficient searching.

4. **Associate Invariants with Star Indices**
   
   - Link each set of invariants to the indices of the stars forming the quad.

```python
>>> # Choose the mode of geometric invariants to calculate: 'triangles' or 'quads'
>>> mode_invariants = 'triangles'

>>> # Calculate the geometric invariants
>>> sc_simplified_stars.invariantfeatures(mode_invariants)

>>> # Retrieve the generated invariants, asterisms, and 2D/4D-Tree structure
>>> invariants = sc_simplified_stars.invariants
>>> asterisms = sc_simplified_stars.asterisms
>>> kdtree = sc_simplified_stars.kdtree
>>> print(invariants,asterisms,kdtree)
```

### Visualization

STARQUERY supports visualizing both **conical** and **rectangular** search areas, highlighting the HEALPix tiles that intersect with the search region.

```python
>>> # Visualize the coverage of the box search area
>>> sc_simplified_stars.tiles_draw()
```

<p align="middle">
  <img src="readme_figs/box.png" width="500" />
</p>

```python
>>> # Visualize the coverage of the cone search area
>>> sc_simplified_stars.tiles_draw()
```

<p align="middle">
  <img src="readme_figs/cone.png" width="500" />
</p>

### Sky Division Using HEALPix

STARQUERY pre-divides the celestial sphere into **multi-level, equal-area sky regions** using the **HEALPix algorithm**. This hierarchical tiling system is used to optimize star catalog queries and facilitate blind star map matching.

**Key Features of HEALPix-based Division:**

- Divides the sky into **hierarchical levels** (from **K1** to **K11**).

- For each level, extracts the brightest stars (e.g., top 5 stars) and integrates them into the **K1** level.

- Calculates geometric invariants for star patterns and stores them in **h5-formatted hash files**.

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

The hash files store precomputed 

- **Center Pointing**: RA/Dec of the sky region center.

- **Pixel Coordinates**: Projected coordinates of stars in the region.

- **Geometric Invariants**: Features calculated for triangles or quads.

- **Asterism Indices**: Indices of stars forming each asterism.

### Read the h5-formatted Hash File

Once the hash file is generated, you can load and access the stored data for efficient blind star pattern recognition.

```python
>>> # Path to the hash file
>>> infile_h5 = 'starcatalogs/indices/at-hyg24_mag13.0_epoch2019.5_triangles_K1_K6.h5'
>>> # Load the hash data
>>> sc_hashed = sc_simplified.read_h5_hashes(infile_h5)
>>> print(sc_hashed.hashed_data)
```

## 🔧 Change log

- **1.1.6 — Apr 29, 2025**
  
  - Remove astropy.time from `astrometry_corrections` to speed up the calculations.
  
  - This may be the last version.


- **1.1.5 — Nov 18, 2024**
  
  - Fixed the problem caused by the rectangle near the celestial pole degenerating into a triangle in the spherical rectangle query.
  
  - Polished the usage documentation.

- **1.1.4 — Oct 30, 2024**
  
  - Set the maximum display records for DataFrames to 200 by default.

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

# 🤝 Contributing

We welcome contributions to the STARQUERY project and are grateful for every bit of help. 

- **Bug Reports**: If you find a bug, please create an issue in our issue tracker. Be sure to include detailed information about the bug and steps to reproduce it.
- **Feature Requests**: If you have ideas for new features or improvements, feel free to open an issue to discuss them.

# 📄 Reference

The development of STARQUERY has been influenced and supported by a variety of external resources and tools. Below are some of the key references:

- [**STScI Outerspace**](https://outerspace.stsci.edu/display/GC/WebServices+for+Catalog+Access): The Space Telescope Science Institute provides web services for catalog access, which have been instrumental in the development of STARQUERY.
- [**The Astronomy Nexus**](https://www.astronexus.com/hyg): A valuable resource providing comprehensive star catalog data, such as  the HYG and AT-HYG database. STARQUERY utilizes this resource for accessing detailed astronomical data.
- [**Astroalign**](https://astroalign.quatrope.org/en/latest/index.html): This Python package by Beroiz, M. I. is a significant reference, especially in the development of star pattern recognition features within STARQUERY.
- [**HEALPix**](https://healpix.sourceforge.io): The Hierarchical Equal Area isoLatitude Pixelization (HEALPix) tool has been a reference for implementing the division of the celestial sphere into equal-area sky regions. Learn more about HEALPix at their website.
