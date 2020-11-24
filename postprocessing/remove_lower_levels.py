"""
Script to delete lower measurements in the level 2
data as the soundings from the ships cannot be
trusted below 40 m
"""
import argparse
import glob
import numpy as np
import xarray as xr
import tqdm

# create parser
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--inputfilefmt", metavar="INPUT_FILE_FMT",
                        help="Level2 files or fileformat \n"
                             "including wildcards e.g. $EXPORT_PATH/level_2/EUREC4A*.nc",
                        default=None,
                        required=True,
                        nargs='+')

# parse the arguments
args = vars(parser.parse_args())

if len(args['inputfilefmt']) == 0:
    level2_files = sorted(glob.glob(args['inputfilefmt'][0]))
else:
    level2_files = sorted(args['inputfilefmt'])

for file in tqdm.tqdm(level2_files):
    print(file)
    ds_in = xr.load_dataset(file)

    if 'BCO' in file:
        #  No remove of lower levels in case of BCO soundings
        ds_out = ds_in
    else:
        ds_ = ds_in
        lowest_mask = np.where(ds_in.altitude<40, True, False)
        for var in ds_in.data_vars:
            if (len(ds_in[var].dims) == 2) and ('_FillValue' in ds_in[var].encoding.keys()):
                if ds_in[var].dims[0] == 'sounding':
                    ds_[var].values[:,lowest_mask] = np.nan
                elif ds_in[var].dims[1] == 'soudning':
                    ds_[var].values[lowest_mask,:] = np.nan
        ds_out = ds_

    for var in ds_out.data_vars:
        ds_out[var].encoding = {'zlib': True}

    # Delete some attributes
    del ds_out.attrs['surface_altitude']
    try:
        del ds_out.attrs['python_version']
    except:
        pass
    ds_out.to_netcdf(file+'2', unlimited_dims=['sounding'])
