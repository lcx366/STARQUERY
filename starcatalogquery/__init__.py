import pandas as pd

from .classes import StarCatalog
from .utils import data_prepare

# Load Earth Orientation Parameters (EOP) and leap seconds.
data_prepare.iers_load()

# Load the JPL Solar System Planetary Ephemeris (DE440S),
# providing accurate positions of solar system bodies for celestial computations.
data_prepare.sspe_load('DE440S')

# Configure pandas display settings for improved DataFrame readability.
pd.set_option('display.max_columns', None)       # Show all columns in DataFrame outputs
pd.set_option('display.precision', 8)            # Use 8-digit precision for floating-point values
pd.set_option('display.width', 1000)             # Set the maximum width of DataFrame outputs
pd.set_option('display.max_rows', 200)           # Limit maximum rows displayed in DataFrame output

# This initialization file prepares the starcatalogquery package by:
# - Loading Earth orientation and leap second data for accurate time transformations,
# - Initializing planetary ephemerides for solar system object positioning,
# - Configuring pandas for consistent and detailed tabular output during analysis and debugging.