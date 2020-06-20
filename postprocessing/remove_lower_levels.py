"""
Script to delete lower measurements in the level 2
data as the soundings from the ships cannot be
trusted below 40 m
"""

level2_files = [
                '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/EUREC4A_Meteor_soundings.nc',
                '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/EUREC4A_BCO_soundings.nc',
                '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/EUREC4A_RH-Brown_soundings.nc',
                '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/EUREC4A_MS-Merian_soundings.nc',
                '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/EUREC4A_Atalante_soundings_Vaisala.nc',
                '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export/level_2/EUREC4A_Atalante_soundings_MeteoModem.nc'
                ]
import numpy as np
import xarray as xr

for file in level2_files:
    print(file)
    ds_in = xr.open_dataset(file)
    if file == '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/EUREC4A_BCO_soundings.nc':
        ds_out = ds_in
    else:
        ds_ = ds_in
        lowest_mask = np.where(ds_in.altitude<40, True, False)
        for var in ds_in.data_vars:
            if len(ds_in[var].dims) == 2:
                if ds_in[var].dims[0] == 'sounding':
                    ds_[var][0,lowest_mask] = np.nan
                else:
                    ds_[var][lowest_mask,0] = np.nan
        ds_out = ds_

    for var in ds_out.data_vars:
        ds_out[var].encoding = {'zlib': True}

    # Delete some attributes
    del ds_out.attrs['surface_altitude']
    try:
        del ds_out.attrs['python_version']
    except:
        pass
    ds_out.to_netcdf(file+'2')
