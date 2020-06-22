"""
Script to rename variable names in level1 and level2
files
"""

import glob
import tqdm
import numpy as np
import xarray as xr

files_l2 = sorted(glob.glob('/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EURE*.nc'))
files_l1 = sorted(glob.glob('/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_1/*/*/*.nc'))

rename_dict = {
    'altitude': 'alt',
    'ascentRate': 'dz',
    'dewPoint': 'dp',
    'flight_time': 'flight_time',
    'humidity': 'rh',
    'latitude': 'lat',
    'longitude': 'lon',
    'mixingRatio': 'mr',
    'mixing_ratio': 'mr',
    'relative_humidity': 'rh',
    'specific_humidity': 'q',
    'pressure': 'p',
    'temperature': 'ta',
    'windDirection': 'wdir',
    'windSpeed':'wspd',
    'wind_direction': 'wdir',
    'wind_speed': 'wspd',
    'dew_point': 'dp',
    'wind_u': 'u',
    'wind_v': 'v'
}

attrs_to_delete = ['surface_altitude', 'latitude_of_launch_location', 'longitude_of_launch_location', 'python_version', 'converted_by', 'contact_person', 'institution', 'location'] 

files = np.hstack([files_l1, files_l2])

for file in tqdm.tqdm(files):
    #if file == '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export/level_1/ATL/Vaisala/EUREC4A_ATL_sounding_ascent_20200127_1059.nc':
    #    import pdb; pdb.set_trace()
    ds = xr.open_dataset(file)
    ds.load()
    for var in list(ds.variables):
        if var in rename_dict.keys():
            ds = ds.rename({var: rename_dict[var]})

    for var in list(ds.variables):
        ds[var].encoding['zlib'] = True
        if '_FillValue' in ds[var].encoding.keys():
            if np.isnan(ds[var].encoding['_FillValue']):
                ds[var].encoding['_FillValue'] = 9.96921e36
            if ds[var].encoding['_FillValue'] == 9.969209968386869e36:
                ds[var].encoding['_FillValue'] = 9.96921e36
        if 'coordinates' in ds[var].attrs.keys():
            del ds[var].attrs['coordinates']
            ds[var].encoding['coordinates'] = "flight_time lat lon" 
    attrs = list(ds.attrs.keys())
    for attr in attrs:
        if attr in attrs_to_delete:
            del ds.attrs[attr]
    for var in list(ds.variables):
        if 'coordinates' in ds[var].encoding.keys():
            ds[var].encoding['coordinates'] = "flight_time lat lon"
    ds.flight_time.encoding['units'] = "seconds since 1970-01-01 00:00:00 UTC"
    ds.launch_time.encoding['units'] = "seconds since 1970-01-01 00:00:00 UTC"

    ds.to_netcdf(file)

