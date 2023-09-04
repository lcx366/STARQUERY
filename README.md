# Welcome to the STARQUERY package

[![PyPI version shields.io](https://img.shields.io/pypi/v/starcatalogquery.svg)](https://pypi.python.org/pypi/starcatalogquery/) [![PyPI pyversions](https://img.shields.io/pypi/pyversions/starcatalogquery.svg)](https://pypi.python.org/pypi/starcatalogquery/) [![PyPI status](https://img.shields.io/pypi/status/starcatalogquery.svg)](https://pypi.python.org/pypi/starcatalogquery/) [![GitHub contributors](https://img.shields.io/github/contributors/lcx366/STARQUERY.svg)](https://GitHub.com/lcx366/STARQUERY/graphs/contributors/) [![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/lcx366/STARQUERY/graphs/commit-activity) [![GitHub license](https://img.shields.io/github/license/lcx366/STARQUERY.svg)](https://github.com/lcx366/STARQUERY/blob/master/LICENSE) [![Documentation Status](https://readthedocs.org/projects/starcatalogquery/badge/?version=latest)](http://starcatalogquery.readthedocs.io/?badge=latest) [![Build Status](https://travis-ci.org/lcx366/starcatalogquery.svg?branch=master)](https://travis-ci.org/lcx366/starcatalogquery)

This package is an archive of scientific routines for establishing an offline star catalog query database. Currently, operations on star catalogue query include:

1. From the remote server ([MIKULSKI ARCHIVE&
   SPACE TELESCOPES](https://archive.stsci.edu)) to obtain star catalog data, and build an offline star catalog database locally;
2. Simplify the original star catalogs to reduce the storage volume of star catalogs;
3. Filter out the required stars according to the selected sky area, including rectangle search and spherical cap(cone) search;
4. Given a pixel width, calculate the pixel coordinates of the filtered stars;
5. Visualize the scope of the search area and the coverage of the corresponding catalog tiles.
6. According to the HEALPix algorithm, the celestial sphere is divided into multiple sky areas, and a feature library is established in each sky area for blind matching of star maps.
7. Load the local offline star catalog database;

## How to Install

On Linux, macOS and Windows architectures, the binary wheels can be installed using pip by executing one of the following commands:

```
pip install starcatalogquery
pip install starcatalogquery --upgrade # to upgrade a pre-existing installation
```

## How to use

### Build an offline star catalog database

Get the star catalog data from the remote server ([MIKULSKI ARCHIVE&
SPACE TELESCOPES](https://archive.stsci.edu)) , and build an offline star catalog database locally. Currently, available star catalogs include 'hygv35', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.

```python
>>> from starcatalogquery import StarCatalog
>>> # Get the star catalog HYGv3.5 with the tile size of 5 deg
>>> hygv35_raw = StarCatalog.get('hygv35',5)
```

We get an instance of class StarCatalogRaw with

    Attributes:
        - tiles_dir: Directory of the star catalog files
        - tiles_num: Total number of the tile files
        - tile_size: Geometric size of the tile in [deg]
        - sc_size: File size of the star catalog
        - validity: Validity of the star catalog
        - sc_name: Name of the star catalog. Available options include 'hygv35', 'gsc12', 'gsc242', 'gaiadr3', '2mass', 'ucac5', 'usnob', etc.
        - _mode: Type of the star catalog. Available options include 'raw', 'reduced', 'simplified', where
            'raw' represents the original star catalog, which covers all information about the stars,
            'reduced' represents the reduced star catalog, which contains the position, proper motion, apparent magnitude, epoch of the stars,
            'simplified' represents the minimalist star catalog, which only includes the position and apparent magnitude of the stars at a specific epoch.
        - stars_num: Total number of stars in the catalog
        - mag: Range of apparent magnitudes for the catalog
        - description: Catalog Summary
    Methods:
        - load: Load the raw star catalog files from the local database.
        - reduce: Reduce the raw star catalog so that it only contains necessary information of stars such as the celestial position, proper motion, apparent magnitude, and epoch, etc.
        - search_box: Perform a rectangle search on the raw star catalog.
        - search_cone: Perform a cone search on the raw star catalog.
        - _search_draw: Visualize the scope of the search area and the coverage of the corresponding tiles.  

### Reduce the raw star catalogs

The original star catalog contains all the catalog information of the stars, which takes up a huge storage space and slows down the speed of star query. Therefore, we need to extract the necessary information of stars from it, and use the method `.reduce()` to build a reduced star catalog database.

```python
>>> hygv35_reduced = hygv35_raw.reduce()
```

The reduced star catalog only contains the celestial position, proper motion, apparent magnitude, and epoch of stars.

### Simplify the reduced star catalogs

In order to further reduce the size of the star catalog and improve its query efficiency, we filter out the reduced star catalog according to the limit magnitude of the detector, and make proper motion corrections to obtain a minimalist star catalog.

```python
>>> mag_threshold,t_pm = 9.0,2022.0 
>>> hygv35_simplified = hygv35_reduced.simplify(mag_threshold,t_pm)
```

The simplified(minimalist) star catalog only includes the celetial position and apparent magnitude of stars at a specific epoch.

### Query information about stars in a specific sky area

Perform a cone search of stars on the raw or reduced star catalog.

```python
>>> # Set center pointing in format of [Ra,Dec] in [deg] and search radius in [deg]
>>> center,radius = [20,30],15 
>>> # Set cutoff magnitude and observation epoch
>>> mag_threshold,t_pm = 9.0,2022.0
>>> # Set the maximum number of brightest stars to output
>>> max_num = 100 # optinal, default = None
>>> hygv35_raw_stars = hygv35_raw.search_cone(center,radius,mag_threshold,t_pm,max_num)
>>> hygv35_reduced_stars = hygv35_reduced.search_cone(center,radius,mag_threshold,t_pm,max_num)
```

Perform a rectangle search of stars on the raw or reduced star catalog.

```python
>>> # Set a rectangular search area in format of [ra_min,dec_min,ra_max,dec_max] in [deg]
>>> radec_box = [5,15,35,45]
>>> # Set cutoff magnitude and observation epoch
>>> mag_threshold,t_pm = 9.0,2022.0
>>> # Set the maximum number of brightest stars to output
>>> max_num = 100 # optinal, default = None
>>> hygv35_raw_stars = hygv35_raw.search_box(radec_box,mag_threshold,t_pm,max_num)
>>> hygv35_reduced_stars = hygv35_reduced.search_box(radec_box,mag_threshold,t_pm,max_num)
```

Perform a cone/rectangle search of stars on the simplified star catalog.

```python
>>> hygv35_simplified_stars = hygv35_simplified.search_cone(center,radius,max_num)
>>> hygv35_simplified_stars = hygv35_simplified.search_box(radec_box,max_num)
```

We get an instance of class Stars, with

    Attributes:
        - sc_name: Name of the starcatalog.
        - center: Center pointing in format of [ra_c,dec_c] in [deg]
        - df: Dataframe of the stars
        - max_num: Number of stars in the dataframe
        - radec: Celestial coordinates of stars
    Methods:
        - pixel_xy: Calculate the pixel coordinates of stars in a sky area.

### Calculate the pixel coordinates of the filtered stars

Given the pixel width of the detector, calculate the pixel coordinates of the filtered stars using the *TAN* projection in WCS transformations.

```python
>>> pixel_width = 0.01 # pixel width in [deg]
>>> hygv35_simplified_stars.pixel_xy(0.01)
>>> print(hygv35_simplified_stars.xy)
```

### Calculate the triangle invariants and the asterism indices of the filtered stars

```python
>>> hygv35_simplified_stars.invariantfeatures()
>>> print(hygv35_simplified_stars.invariants,hygv35_simplified_stars.asterisms,hygv35_simplified_stars.kdtree)
```

### Visualization

Visualize the scope of the search area and the coverage of the corresponding catalog tiles.

```python
>>> box_area = {'box':[5,15,35,45]} # {'box':[ra_min,dec_min,ra_max,dec_max]}
>>> cone_area = {'cone':[20,30,15]} # {'cone':[ra_c,dec_c,radius]}
>>> hygv35_simplified._search_draw(box_area)
>>> hygv35_simplified._search_draw(cone_area)
>>> # ._search_draw is also available for hygv35_raw and hygv35_reduced
```

<p align="middle">
  <img src="readme_figs/box.png" width="500" />
</p>

<p align="middle">
  <img src="readme_figs/cone.png" width="500" />
</p>

### Divide the sky into multiple equal-area sky regions

For blind matching of star maps, the celestial sphere is pre-divided into multiple equal-area sky regions using the HEALPix algorithm.

<p align="middle">
  <img src="readme_figs/healpix_list.png" width="400" />
</p>

<p align="middle">
  <img src="readme_figs/healpix_table.png" width="400" />
</p>

By default, the following strategy is adopted to divide the sky area:

- For 58.6 > FOV >= 29.3, k = 1, radius of cone search = 37.2; 

- For 29.3 > FOV >= 14.7, k = 2, radius of cone search = 17.0;

- For 14.7 > FOV >= 7.33, k = 3, radius of cone search = 8.3;

- For 7.33 > FOV >= 3.66, k = 4, radius of cone search = 4.1;

- For 3.66 > FOV >= 1.83, k = 5, radius of cone search = 2.1;

```python
>>> fov,pixel_width = 8,0.01 # in [deg]
>>> # Set the maximum number of brightest stars in each sky area
>>> max_num = 30 
>>> outh5,ratio = hygv35_simplified.h5_indices(fov,pixel_width,max_num) # Ratio of the solid angles spanned by the square and the cone
```

A h5-formatted star catalog indices file is generated, which records the center pointing, pixel coordinates of the stars, triangle invariants and asterism indices of each sky area.

### Read and parse the h5-formatted star catalog indices file

```python
>>> from starcatalogquery import StarCatalog
>>> infile_h5 = 'starcatalogs/indices/hygv35/k2_mag9.0_mcp30_2022.0.h5'
>>> fp_radecs,stars_xy,stars_invariants,stars_asterisms = StarCatalog.read_h5_indices(infile_h5)
```

### Load the local offline star catalog database

#### Load the raw or reduced star catalog

```python
>>> from starcatalogquery import StarCatalog
>>> dir_from_raw = '/Volumes/TOSHIBA/starcatalogs/raw/hygv35/res5/' # Path of the raw starcatalog
>>> hygv35_raw = StarCatalog.load(dir_from_raw)
>>> # dir_from_reduced = '/Volumes/TOSHIBA/starcatalogs/reduced/hygv35/res5/' # Path of the reduced starcatalog
>>> # hygv35_reduced = StarCatalog.load(dir_from_reduced)
```

#### Load the simplified star catalog

```python
>>> from starcatalogquery import StarCatalog
>>> dir_from_simplified = '/Volumes/TOSHIBA/starcatalogs/simplified/hygv35/res5/mag9.0/epoch2022.0/' # Path of the starcatalog
>>> hygv35_simplified = StarCatalog.load(dir_from_simplified)
```

## Change log

- **0.1.11 — Sep 04, 2023**
  
  - Adjusted the parameters involved in dividing the celestial sphere into multiple equal-area sky regions using the HEALPix algorithm, as well as the corresponding radius of cone search used for blind matching of star maps.
  - Adjusted the arrangement of the equal-area sky regions so that it gradually moves from the celestial equator to the poles.
  - minor bugs fixed

- **0.1.10 — Jul 23, 2023**
  
  - Added the the HEALPix algorithm for dividing the celestial sphere into multiple equal-area sky regions.

- **0.1.8 — Jul 03, 2023**
  
  - Get rid of dependency on **pyshtools**
    - Added function `from_cap` to replace `pysh.SHGrid.from_cap` in *catalog_query.py*
    - Using `SphericalCircle` in *astropy* to replace `MakeCircleCoord` in *pyshtools.utils* 

- **0.1.7 — Jun 16, 2023**
  
  - Simplified parameter input of `StarCatalog.load` for star catalog loading.

- **0.1.5 — May 13, 2023**
  
  - Add method `.invariantfeatures()` to class `Stars`, which calculates the triangle invariants and constructs a 2D Tree; and records the asterism indices for each triangle.

- **0.1.0 — May 10,  2023**
  
  - The ***starcatalogquery*** package was released.

## Reference

- [MIKULSKI ARCHIVE&SPACE TELESCOPES](https://archive.stsci.edu)
- [Astroalign](https://astroalign.quatrope.org/en/latest/index.html)
- [HEALPix](https://healpix.sourceforge.io)
