import os
import numpy as np
import pymeteo.skewt as skewt
import pymeteo.thermo as met
import metpy.calc as mpcalc
from metpy.units import units
import xarray as xr
import sys

file = sys.argv[0]

# Read data
ds = xr.open_dataset(file)
p = ds['pressure'].isel({'sounding': 0}).values * units.hPa
z = ds['altitude'].isel({'sounding': 0}).values * units.meters
T = ds['temperature'].isel({'sounding': 0}).values * units.degC
Td = ds['dewPoint'].isel({'sounding': 0}).values * units.degC
wind_speed = ds['windSpeed'].isel({'sounding': 0}).values * (units.meter/units.second)
wind_dir = ds['windDirection'].isel({'sounding': 0}).values * units.degrees

# Calculate input for skewt
u, v = mpcalc.wind_components(wind_speed, wind_dir)
th = met.theta(np.array(T)+met.T00, np.array(p)*100)
w = met.es(np.array(Td)+met.T00) / met.es(np.array(T)+met.T00)
pp = met.es(np.array(T)+met.T00) / (np.array(p)*100)
qv = 0.622 * pp * w

# Create skewt

skewt.plot(None, np.array(z), th, np.array(p)*100,
           qv, np.array(u), np.array(v), 'sounding.pdf',
           title=os.path.basename(file))
