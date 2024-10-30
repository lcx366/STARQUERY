import pandas as pd

from .classes import StarCatalog, CatalogDB
from .utils import data_prepare

# Load Earth Orientation Parameters (EOP) and Leap Seconds (LS) files
data_prepare.iers_load()

# Load the JPL Solar System Planetary Ephemeris (SSPE) file
data_prepare.sspe_load('DE440S')

# Configure pandas display options
pd.set_option('display.max_columns', None) # Display all columns in DataFrame outputs
pd.set_option("display.precision", 8) # Set the default floating-point precision to 8 digits
pd.set_option('display.width', 1000) # Set the maximum display width for DataFrames
pd.set_option('display.max_rows', 200) # Set the maximum display records for DataFrames

# The __init__.py file initializes the starcatalogquery package by loading essential data and configuring pandas settings.
# The iers_load function loads Earth Orientation Parameters and Leap Seconds, which are crucial for accurate astronomical calculations.
# The sspe_load function loads the JPL Solar System Planetary Ephemeris, providing precise positions of solar system bodies.
# Pandas options are configured to ensure a comprehensive and precise display of DataFrame outputs, facilitating data analysis and debugging.
