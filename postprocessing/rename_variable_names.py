"""
Script to rename variable names in level1 and level2
files
"""

import argparse
import glob
import tqdm
import numpy as np
import xarray as xr

# create parser
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--inputfilefmt", metavar="INPUT_FILE_FMT",
                        help="Level1 and 2 files or fileformat \n"
                             "including wildcards e.g. $EXPORT_PATH/level_2/EUREC4A*.nc",
                        default=None,
                        required=True,
                        nargs='+')

# parse the arguments
args = vars(parser.parse_args())

if len(args['inputfilefmt']) == 0:
    files = sorted(glob.glob(args['inputfilefmt'][0]))
else:
    files = sorted(args['inputfilefmt'])

rename_dict = {
    'altitude': 'alt',
    'ascentRate': 'dz',
    'ascent_rate': 'dz',
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

rename_attrs_dict = {
    'git-version': 'version',
    'git_version': 'version'
    }

attrs_to_delete = ['platform_location', 'surface_altitude', 'latitude_of_launch_location', 'longitude_of_launch_location', 'python_version', 'converted_by', 'contact_person', 'institution', 'location'] 

vars_to_delete = ['altitude_WGS84']

for file in tqdm.tqdm(files):
    #if file == '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export/level_1/ATL/Vaisala/EUREC4A_ATL_sounding_ascent_20200127_1059.nc':
    #    import pdb; pdb.set_trace()
    ds = xr.open_dataset(file)
    ds.load()
    for var in list(ds.variables):
        if var in rename_dict.keys():
            ds = ds.rename({var: rename_dict[var]})

    for attr in list(ds.attrs.keys()):
        if attr in rename_attrs_dict.keys():
            ds.attrs[rename_attrs_dict[attr]] = ds.attrs[attr]
            del ds.attrs[attr]

    for var in list(ds.variables):
        ds[var].encoding['zlib'] = True
        if '_FillValue' in ds[var].encoding.keys():
            if np.isnan(ds[var].encoding['_FillValue']):
                ds[var].encoding['_FillValue'] = 9.96921e36
            if ds[var].encoding['_FillValue'] == 9.969209968386869e36:
                ds[var].encoding['_FillValue'] = 9.96921e36
        if 'coordinates' in ds[var].attrs.keys():
            del ds[var].attrs['coordinates']
            ds[var].encoding['coordinates'] = "sounding flight_time lat lon"
    attrs = list(ds.attrs.keys())
    for attr in attrs:
        if attr in attrs_to_delete:
            del ds.attrs[attr]
    for var in vars_to_delete:
        if var in ds.data_vars:
            del ds[var]
    for var in list(ds.variables):
        if 'coordinates' in ds[var].encoding.keys():
            ds[var].encoding['coordinates'] = "sounding flight_time lat lon"
    ds.flight_time.encoding['units'] = "seconds since 2020-01-01 00:00:00"
    ds.launch_time.encoding['units'] = "seconds since 2020-01-01 00:00:00"

    ds.to_netcdf(file, unlimited_dims=['sounding'])

