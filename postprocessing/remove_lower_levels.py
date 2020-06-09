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
import xarray as xr

for file in level2_files:
    ds_in = xr.open_dataset(file)
    if file == '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/EUREC4A_BCO_soundings.nc':
        ds_out = ds_in
    else:
        ds_out = ds_in.sel(altitude=slice(40,30990))

    for var in ds_out.data_vars:
        ds_out[var].encoding = {'zlib': True}

    # Delete some attributes
    del ds_out.attrs['surface_altitude']
    del ds_out.attrs['python_version']

    ds_out.to_netcdf(file+'2')
