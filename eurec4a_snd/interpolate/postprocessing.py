"""
Functions for post-processing data
"""

import numpy as np
import metpy.calc
from metpy.units import units

def calc_saturation_pressure(temperature_K, method='hardy1988'):
    """
    Calculate saturation water vapor pressure
    
    Input
    -----
    temperature_K : array
        array of temperature in Kevlin or dew point temperature for actual vapor pressure
    method : str
        Formula used for calculating the saturation pressure
            'hardy1988' : ITS-90 Formulations for Vapor Pressure, Frostpoint Temperature,
                Dewpoint Temperature, and Enhancement Factors in the Range –100 to +100 C,
                Bob Hardy, Proceedings of the Third International Symposium on Humidity and Moisture,
                1998 (same as used in Aspen software after May 2018)
    
    Return
    ------
    e_sw : array
        saturation pressure (Pa)
    
    Examples
    --------
    >>> calc_saturation_pressure([273.15])
    array([ 611.2129107])
    
    >>> calc_saturation_pressure([273.15, 293.15, 253.15])
    array([  611.2129107 ,  2339.26239586,   125.58350529])
    """

    if method == 'hardy1988':
        g = np.empty(8)
        g[0] = -2.8365744*10**3
        g[1] = -6.028076559*10**3
        g[2] = 1.954263612*10**1
        g[3] = -2.737830188*10**(-2)
        g[4] = 1.6261698*10**(-5)
        g[5] = 7.0229056*10**(-10)
        g[6] = -1.8680009*10**(-13)
        g[7] = 2.7150305

        e_sw = np.zeros_like(temperature_K)

        for t, temp in enumerate(temperature_K):
            ln_e_sw = np.sum([g[i]*temp**(i-2) for i in range(0,7)]) + g[7]*np.log(temp)
            e_sw[t] = np.exp(ln_e_sw)
        return e_sw


def calc_mixing_ratio_hardy(dew_point_K, pressure_Pa):
    e_s = calc_saturation_pressure(dew_point_K)
    mixing_ratio = metpy.calc.mixing_ratio(e_s * units.pascal, pressure_Pa * units.pascal)
    return mixing_ratio


def set_global_attributes(ds, global_attrs_dict):
    """
    Set global attributes
    """
    for attribute, value in global_attrs_dict.items():
        ds.attrs[attribute] = value
    return ds


def compress_dataset(ds):
    """
    Apply internal netCDF4 compression
    """
    for var in ds.data_vars:
        ds[var].encoding['zlib'] = True
    return ds


def set_additional_var_attributes(ds, meta_data_dict):
    """
    Set further descriptive variable
    attributes for interpolated variables
    """
    for var in ds.variables:
        try:
            meta_data_var = meta_data_dict[var]
        except KeyError:
            continue
        for key, value in meta_data_var.items():
            ds[var].attrs[key] = value
    return ds


def write_dataset(ds, filename):
    ds = compress_dataset(ds)
    ds.to_netcdf(filename, unlimited_dims=['sounding'])


def get_direction(ds_interp, ds):
    direction = None
    # First source of direction
    direction_msg = None
    try:
        if '309057' in ds.source or '309052' in ds.source:
            direction_msg = 'ascending'
        elif '309056' in ds.source or '309053' in ds.source:
            direction_msg = 'descending'
    except AttributeError:
        print('No attribute source found')

    # Second source of direction
    median_ascent = np.nanmedian(ds.ascentRate.values)
    if median_ascent > 0:
        direction_data = 'ascending'
    elif median_ascent < 0:
        direction_data = 'descending'
    else:
        print('Direction not retrievable')
    
    if (direction_msg == direction_data) or (direction_msg is None):
        return direction_data
    else:
        print('Direction mismatch')

def calc_ascentrate(height, time):
    """
    Calculate the ascent rate

    Input
    -----
    sounding : obj
        sounding class containing gpm
        and flight time

    Return
    ------
    soundning : obj
        sounding including the ascent rate
    """
    ascent_rate = np.diff(height)/(np.diff(time))
    ascent_rate = np.ma.concatenate(([0], ascent_rate))  # 0 at first measurement
    return ascent_rate


