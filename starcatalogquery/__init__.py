import pandas as pd

from .classes import StarCatalog, CatalogDB
from .utils import data_prepare

# Load the Earth Orientation Parameters (EOP) and Leap Seconds (LS) file
data_prepare.iers_load()

# Load the JPL Solar System Planetary Ephemeris (SSPE)
data_prepare.sspe_load('DE440S')

# Configure pandas options
pd.set_option('display.max_columns', None)  # Display all columns in DataFrame outputs
pd.set_option("display.precision", 8)       # Set the default floating point precision to 8

# The __init__.py file initializes the starcatalogquery package by loading essential data and configuring pandas settings.
# The iers_load function loads Earth Orientation Parameters and Leap Seconds, which are crucial for accurate astronomical calculations.
# The sspe_load function loads the JPL Solar System Planetary Ephemeris, used for precise solar system body positions.
# Pandas options are set to ensure comprehensive and precise display of DataFrame outputs, which is useful for data analysis and debugging.
