"""
Script to change units
"""
import os
import glob
import tqdm
import argparse
import numpy as np
import xarray as xr
import logging

logging.basicConfig(filename='change_units.log', level=logging.INFO)


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inputfilefmt', required=True,
                    help='Input filename format (level 1 files), e.g. EUREC4A_*soundings_*.nc')
args = vars(parser.parse_args())

files = sorted(glob.glob(args['inputfilefmt']))

logging.info('Total number of soundings found: {}'.format(len(files)))


def hPa2Pa(hPa):
    """
    Convert hPa to Pa
    """
    assert np.all(hPa[~np.isnan(hPa)] < 1100), 'Is this really hPa?'
    return hPa * 100.


def gkg2kgkg(gkg):
    return gkg / 1000.


def degC2K(degC):
    assert np.all(degC[~np.isnan(degC)] < 150), 'Is this really deg C?'
    return degC + 273.15


def percent2frac(percent):
    return percent/100.


units_convert = {'p': {'units_old': 'hPa',
                       'units_new': 'Pa',
                       'converter': hPa2Pa},
                'mr': {'units_old': 'g/kg',
                       'units_new': 'kg/kg',
                       'converter': gkg2kgkg},
                'q': {'units_old': 'g/kg',
                       'units_new': 'kg/kg',
                       'converter': gkg2kgkg},
                'dp': {'units_old': 'degrees_Celsius',
                       'units_new': 'K',
                       'converter': degC2K},
                'ta': {'units_old': 'degrees_Celsius',
                       'units_new': 'K',
                       'converter': degC2K},
                'rh': {'units_old': '%',
                       'units_new': '1',
                       'converter': percent2frac},
                'dp': {'units_old': 'degrees_Celsius',
                       'units_new': 'K',
                       'converter': degC2K}
                }


for f, file in enumerate(tqdm.tqdm(files)):
    ds = xr.load_dataset(file)
    for var in ds.data_vars:
        if var in units_convert.keys():
            var_dict = units_convert[var]
            if ds[var].attrs['units'] == var_dict['units_old']:
                ds[var].values = var_dict['converter'](ds[var].values)
                ds[var].attrs['units'] = var_dict['units_new']
    ds['flight_time'].encoding['_FillValue'] = False
    ds.to_netcdf(file, unlimited_dims=['sounding'])

