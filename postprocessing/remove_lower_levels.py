"""
Script to delete lower measurements in the level 2
data as the soundings from the ships cannot be
trusted below 40 m
"""

level2_files = [
                '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_Meteor_soundings.nc',
                '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_BCO_soundings.nc',
                '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_RonBrown_soundings.nc',
                '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_MS-Merian_soundings.nc',
                '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_Atalante_soundings_Vaisala.nc',
                '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_Atalante_soundings_Meteomodem.nc'
                ]
import os
import numpy as np
import xarray as xr
import tqdm

for file in tqdm.tqdm(level2_files):
    print(file)
    ds_in = xr.load_dataset(file)
    #os.remove(file)
    if file == '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_BCO_soundings.nc':
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
