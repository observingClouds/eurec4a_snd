"""
Helper functions
"""
import tempfile
import os
import re
import inspect
import platform
import json
import datetime as dt
import numpy as np
import logging
from pathlib import Path, PureWindowsPath
import configparser
from configparser import ExtendedInterpolation
from netCDF4 import date2num
import xarray as xr


class UnitChangedError(Exception):
    pass


class UnexpectedUnit(Exception):
    pass


class RegexDict(dict):
    """
    Dictionary with capability of taking regular expressions
    """
    def get_matching(self, event):
        return (self[key] for key in self if re.match(key, event))

    def get_matching_combined(self, event):
        """
        Find matching keys and return combined dictionary

        >>> d = {'EUREC4A_*':{'a':0}, 'EUREC4A_BCO':{'p':1}}
        >>> rd = RegexDict(d)
        >>> rd.get_matching_combined('EUREC4A_BCO')
        {'a': 0, 'p': 1}
        """
        matching = self.get_matching(event)
        dall = {}
        for d in matching:
            dall.update(d)
        return dall


def get_global_attrs(cfg_file, key):
    """
    Get global attributes from configuration file
    """
    with open(cfg_file, 'r') as f:
        j = json.load(f)
        rd = RegexDict(j['global_meta_data'])
        return rd.get_matching_combined(key)


def unixpath(path_in):
    """
    Convert windows path to unix path syntax
    depending on the used OS
    """
    if platform.system() == 'Windows':
        path_out = Path(PureWindowsPath(path_in))
    else:
        path_out = Path(path_in)
    return path_out


def setup_logging(verbose):
    assert verbose in ["DEBUG", "INFO", "WARNING", "ERROR"]
    # Get filename of calling script
    frame = inspect.stack()[1]
    module = inspect.getmodule(frame[0])
    filename = module.__file__

    logging.basicConfig(
        level=logging.getLevelName(verbose),
        format="%(levelname)s - %(name)s - %(funcName)s - %(message)s",
        handlers=[
            logging.FileHandler("{}.log".format(filename)),
            logging.StreamHandler()
        ])


def load_configuration(configuration_file=None):
    """
    Loads the configuration file PATH.ini.
    1. If provided load configuration_file
    2. Attempt to load from home directory
    3. Attempt to load from relative path inside BCO-git structure

    Args:
        configuration_file: optional: complete path to the configuration file.

    Returns:
        instance of ConfigParser class with extended interpolation.
    """
    # Get filename of calling script
    frame = inspect.stack()[1]
    module = inspect.getmodule(frame[0])
    filename = module.__file__

    dir_path = os.path.dirname(os.path.realpath(filename))
    ini_path = "/".join(dir_path.split("/")[:-1]) + "/eurec4a_snd/config/meta_information.ini"
    if not isinstance(configuration_file, str):
        possible_file_in_userdir = Path("~/meta_information.ini").expanduser()
        if os.path.isfile(possible_file_in_userdir):
            configuration_file = possible_file_in_userdir
        elif os.path.isfile(ini_path):
            configuration_file = ini_path
        if configuration_file is None or not os.path.isfile(configuration_file):
            raise FileNotFoundError(
                "No Configuration File 'meta_information.ini' found. Please create one"
                " in your home directory "
                "or provide the path via the argument parsing -c.")
        else:
            logging.info("Using configuration file: %s" % configuration_file)

    conf = configparser.ConfigParser(interpolation=ExtendedInterpolation())
    conf.read(configuration_file)
    return conf


def convert_bufr_to_json(bufr_fn, logger=None):
    """
    Convert bufr file to json with ecCodes
    software
    """
    tmp_folder = tempfile.mkdtemp()
    tmp_output_json_fn = os.path.join(tmp_folder, 'tmp_bufr.json')
    r = os.system("bufr_dump -j s {} > {}".format(bufr_fn, tmp_output_json_fn))
    if logger is None:
        print("Converted {} to {}".format(bufr_fn, tmp_output_json_fn))
    else:
        logger.debug("Converted {} to {}".format(bufr_fn, tmp_output_json_fn))
    return tmp_output_json_fn


def flatten_json(y):
    """
    Flatten structured json array
    
    Input
    -----
    y : list or dict

    Return
    ------
    list : list
        list containing dicts for each key value pair
    """
    global r
    out = {}
    r = 0

    def flatten(x, name='', number=0):
        """helper function"""
        global r
        if type(x) is dict:
            for a in x:
                flatten(x[a], a)
        elif type(x) is list:
            i = 0
            for a in x:
                flatten(a)
                i += 1
                r += 1
        else:
            out['{:07g}_{}'.format(r, name)] = x

    flatten(y)
    return out


def read_json(json_fn):
    """
    Read and flatten json
    """
    with open(json_fn) as file:
        struct_json = json.load(file)

    bfr_json_flat = flatten_json(struct_json)
    keys = bfr_json_flat.keys()
    key_keys = []
    for key in keys:
        if 'key' in key:
            key_keys.append(key[:-4])

    return bfr_json_flat, key_keys


def convert_json_to_arrays(json_flat, key_keys):
    """
    Convert json data to array
    """
    class Sounding:
        """
        Class containing sounding data
        """
        def __init__(self):
            self.station_lat = None
            self.station_lon = None
            self.sounding_start_time = None
            self.time = []
            self.time_unit = None
            self.pressure = []
            self.pressure_unit = None
            self.temperature = []
            self.temperature_unit = None
            self.dewpoint = []
            self.dewpoint_unit = None
            self.windspeed = []
            self.windspeed_unit = None
            self.winddirection = []
            self.winddirection_unit = None
            self.gpm = []
            self.gpm_unit = None
            self.displacement_lat = []
            self.displacement_lat_unit = None
            self.displacement_lon = []
            self.displacement_lon_unit = None
            self.extendedVerticalSoundingSignificance = []
            self.meta_data = {}

        def to_dataset(self):
            """
            Convert sounding to dataset
            """
            date_unit = "seconds since 1970-01-01 00:00:00 UTC"
            var_dict = {
                'launch_time': date2num(self.sounding_start_time, date_unit),
                'flight_time': self.time,
                'altitude': self.gpm,
                'ascentRate': self.ascentrate,
                'dewPoint': self.dewpoint,
                'humidity': self.relativehumidity,
                'latitude': self.latitude,
                'longitude': self.longitude,
                'mixingRatio': self.mixingratio,
                'pressure': self.pressure,
                'temperature': self.temperature,
                'windDirection': self.winddirection,
                'windSpeed': self.windspeed
            }

            ds = xr.Dataset()
            for var, value in var_dict.items():
                if isinstance(value, np.ndarray):
                    ds[var] = xr.DataArray([value], dims=['sounding', 'level'])
                else:
                    ds[var] = xr.DataArray([value], dims=['sounding'])
            return ds


    def _ensure_measurement_integrity(self):
        """
        Test integrity of each measurement unit

        Measurements of the sonde in the bufr file
        contain usually:
            - time since launch
            - pressure
            - gpm
            - location displacement
            - temperature
            - dewpoint
            - wind direction
            - wind speed

        This is a complete unit. However, there are,
        dependining on the bufr format additional
        measurements, which might not consist of
        a complete measurement set. This is for ex.
        the case for the entry which contains the
        "absoluteWindShearIn1KmLayerBelow"

        Here the completeness of the measurement
        is checked and corrected otherwise by
        adding nan values to the not measured
        values.
        """
        if len(self.time) > len(self.temperature):
            self.temperature.append(np.nan)
        if len(self.time) > len(self.gpm):
            self.gpm.append(np.nan)
        if len(self.time) > len(self.dewpoint):
            self.dewpoint.append(np.nan)
        if len(self.time) > len(self.pressure):
            self.pressure.append(np.nan)
        if len(self.time) > len(self.windspeed):
            self.windspeed.append(np.nan)
        if len(self.time) > len(self.winddirection):
            self.winddirection.append(np.nan)
        if len(self.time) > len(self.displacement_lat):
            self.displacement_lat.append(np.nan)
        if len(self.time) > len(self.extendedVerticalSoundingSignificance):
            self.extendedVerticalSoundingSignificance.append(0)
        return

    s = Sounding()

    for key_key in key_keys:
        if json_flat[key_key+'_key'] == 'latitude':
            s.station_lat = json_flat[key_key+'_value']
        elif json_flat[key_key+'_key'] == 'longitude':
            s.station_lon = json_flat[key_key+'_value']
        elif json_flat[key_key+'_key'] == 'pressure':
            s.pressure.append(json_flat[key_key+'_value'])
            if s.pressure_unit is None:
                s.pressure_unit = json_flat[key_key+'_units']
            # Unit consistency test
            elif s.pressure_unit != json_flat[key_key+'_units']:
                raise UnitChangedError('{} and {} are not same unit'.format(s.pressure_unit,
                                                                            json_flat[key_key+'_units']))
        elif json_flat[key_key+'_key'] == 'windSpeed':
            s.windspeed.append(json_flat[key_key+'_value'])
            if s.windspeed_unit is None:
                s.windspeed_unit = json_flat[key_key+'_units']
            # Unit consistency test
            elif s.windspeed_unit != json_flat[key_key+'_units']:
                raise UnitChangedError('{} and {} are not same unit'.format(s.windspeed_unit,
                                                                            json_flat[key_key+'_units']))
        elif json_flat[key_key+'_key'] == 'windDirection':
            s.winddirection.append(json_flat[key_key+'_value'])
            if s.winddirection_unit is None:
                s.winddirection_unit = json_flat[key_key+'_units']
            # Unit consistency test
            elif s.winddirection_unit != json_flat[key_key+'_units']:
                raise UnitChangedError('{} and {} are not same unit'.format(s.winddirection_unit,
                                                                            json_flat[key_key+'_units']))
        elif json_flat[key_key+'_key'] == 'extendedVerticalSoundingSignificance':
            s.extendedVerticalSoundingSignificance.append(json_flat[key_key+'_value'])
        elif json_flat[key_key+'_key'] == 'nonCoordinateGeopotentialHeight':
            s.gpm.append(json_flat[key_key+'_value'])
            if s.gpm_unit is None:
                s.gpm_unit = json_flat[key_key+'_units']
            # Unit consistency test
            elif s.gpm_unit != json_flat[key_key+'_units']:
                raise UnitChangedError('{} and {} are not same unit'.format(s.gpm_unit,
                                                                            json_flat[key_key+'_units']))
        elif json_flat[key_key+'_key'] == 'airTemperature':
            s.temperature.append(json_flat[key_key+'_value'])
            if s.temperature_unit is None:
                s.temperature_unit = json_flat[key_key+'_units']
            # Unit consistency test
            elif s.temperature_unit != json_flat[key_key+'_units']:
                raise UnitChangedError('{} and {} are not same unit'.format(s.temperature_unit,
                                                                            json_flat[key_key+'_units']))
        elif json_flat[key_key+'_key'] == 'dewpointTemperature':
            s.dewpoint.append(json_flat[key_key+'_value'])
            if s.dewpoint_unit is None:
                s.dewpoint_unit = json_flat[key_key+'_units']
            # Unit consistency test
            elif s.dewpoint_unit != json_flat[key_key+'_units']:
                raise UnitChangedError('{} and {} are not same unit'.format(s.dewpoint_unit,
                                                                            json_flat[key_key+'_units']))
        elif json_flat[key_key+'_key'] == 'latitudeDisplacement':
            s.displacement_lat.append(json_flat[key_key+'_value'])
            if s.displacement_lat_unit is None:
                s.displacement_lat_unit = json_flat[key_key+'_units']
            # Unit consistency test
            elif s.displacement_lat_unit != json_flat[key_key+'_units']:
                raise UnitChangedError('{} and {} are not same unit'.format(s.displacement_lat_unit,
                                                                            json_flat[key_key+'_units']))
        elif json_flat[key_key+'_key'] == 'longitudeDisplacement':
            s.displacement_lon.append(json_flat[key_key+'_value'])
            if s.displacement_lon_unit is None:
                s.displacement_lon_unit = json_flat[key_key+'_units']
            # Unit consistency test
            elif s.displacement_lon_unit != json_flat[key_key+'_units']:
                raise UnitChangedError('{} and {} are not same unit'.format(s.displacement_lon_unit,
                                                                            json_flat[key_key+'_units']))
        elif json_flat[key_key+'_key'] == 'timePeriod':
            # Check conistency of data
            _ensure_measurement_integrity(s)
            s.time.append(json_flat[key_key+'_value'])

            if s.time_unit is None:
                s.time_unit = json_flat[key_key+'_units']
            # Unit consistency test
            elif s.time_unit != json_flat[key_key+'_units']:
                raise UnitChangedError('{} and {} are not same unit'.format(s.time_unit,
                                                                            json_flat[key_key+'_units']))
        elif json_flat[key_key+'_key'] == 'year':
            year = json_flat[key_key+'_value']
        elif json_flat[key_key+'_key'] == 'month':
            month = json_flat[key_key+'_value']
        elif json_flat[key_key+'_key'] == 'day':
            day = json_flat[key_key+'_value']
        elif json_flat[key_key+'_key'] == 'hour':
            hour = json_flat[key_key+'_value']
        elif json_flat[key_key+'_key'] == 'minute':
            minute = json_flat[key_key+'_value']
        elif json_flat[key_key+'_key'] == 'second':
            second = json_flat[key_key+'_value']

        # Meta data
        elif json_flat[key_key+'_key'] == 'radiosondeSerialNumber':
            s.meta_data['sonde_serial_number'] = json_flat[key_key+'_value']
        elif json_flat[key_key+'_key'] == 'softwareVersionNumber':
            s.meta_data['softwareVersionNumber'] = json_flat[key_key+'_value']
        elif json_flat[key_key+'_key'] == 'radiosondeType':
            s.meta_data['radiosondeType'] = json_flat[key_key+'_value']
        elif json_flat[key_key+'_key'] == 'unexpandedDescriptors':
            try:
                s.meta_data['bufr_msg'] = json_flat[key_key+'_value']
            except KeyError:
                descriptors = read_unexpandedDescriptors(key_key, json_flat)
                msg_fmt = search_bufr_msg_format(descriptors)
                s.meta_data['bufr_msg'] = msg_fmt
        elif json_flat[key_key+'_key'] == 'radiosondeOperatingFrequency':
            s.meta_data['sonde_frequency'] = str(json_flat[key_key+'_value']) +' '+ json_flat[key_key+'_units']

    _ensure_measurement_integrity(s)

    s.sounding_start_time = dt.datetime(year,
                                        month,
                                        day,
                                        hour,
                                        minute,
                                        second)

    return s


def read_unexpandedDescriptors(current_key, json_flat):
    """
    Reading unexpandedDescriptors in BUFR file.

    unecpandedDescriptor contains the information which
    bufr messages are included in the message
    """

    descriptors = []
    i = int(current_key)
    while True:
        try:
            _ = json_flat["{0:07g}_value".format(i)]
            break
        except KeyError:
            try:
                descriptors.append(json_flat['{0:07g}_'.format(i)])
                i += 1
            except KeyError:
                break
    return descriptors


def search_bufr_msg_format(descriptors):
    """
    Search for common bufr message format in descriptors
    """

    common_descriptors = [309052, 309053, 309055, 309056, 309057]
    for descriptor in descriptors:
        if descriptor in common_descriptors:
            return descriptor
    return None


def replace_missing_data(sounding):
    """
    Removing Nones from sounding measurements
    """
    def replace_none(entry):
        """
        Replace None with NaN
        """
        if entry is None:
            return np.nan
        else:
            return entry

    variables = ['displacement_lat', 'displacement_lon', 'pressure', 'windspeed',
                 'winddirection', 'temperature', 'dewpoint', 'gpm', 'time']

    for var in variables:
        sounding.__dict__[var] = list(map(replace_none, sounding.__dict__[var]))

    return sounding


def decode_extendedVerticalSoundingSignificance(decimal):
    """
    Decode the extendedVerticalSoudingSignificance

    The extededVerticalSoundingSignificance is indicating
    specific height levels and values in an atmospheric
    sounding.
    These height levels might be original levels or slighlty
    interpolated once to describe for example the atmospheric
    conditions at the standard pressure levels.

    This indicator might be a combination of the following
    for a single level:

    1  Surface
    2  Standard level
    3  Tropopause level
    4  Maximum wind level
    5  Significant temperature level
    6  Significant humidity level
    7  Significant wind level
    8  Beginning of missing temperature data
    9  End of missing temperature data
    10 Beginning of missing humidity data
    11 End of missing humidity data
    12 Beginning of missing wind data
    13 End of missing wind data
    14 Top of wind sounding
    15 Level determined by regional decision
    16 Reserved
    17 Pressure level originally indicated by height as the vertical coordinate
    18 Missing value

    Input
    -----
    decimal : integer
        extendedVerticalSoundingSignificance as decimal

    Result
    ------
    bits : array
        list of Significances that are set

    Example
    -------
    >>> decode_extendedVerticalSoundingSignificance(16384)
    array([4])

    >>> decode_extendedVerticalSoundingSignificance(20)
    array([14, 16])
    """
    mask = np.array("{0:018b}".format(decimal), dtype='c').astype(bool)
    bits = np.where(mask)[0]+1
    return bits


def convert_list_to_array(sounding):
    """
    Convert datatype of sounding
    """
    variables = ['displacement_lat', 'displacement_lon', 'pressure', 'windspeed',
                 'winddirection', 'temperature', 'dewpoint', 'gpm', 'time',
                 'extendedVerticalSoundingSignificance']

    for var in variables:
        sounding.__dict__[var] = np.array(sounding.__dict__[var])

    return sounding


def read_bufr_sounding(bufr_file, logging):
    json_file = convert_bufr_to_json(bufr_file, logging)
    json_flat, keys = read_json(json_file)

    sounding = convert_json_to_arrays(json_flat, keys)
    sounding = replace_missing_data(sounding)
    sounding = convert_list_to_array(sounding)

    return sounding


def calculate_coordinates(origin, offset):
    """
    Calculate positon of measurement

    Input
    -----
    origin : float
        latitude or longitude of launch position
    offset : float
        latitudinal or longitudinal displacement
        since launch from origin

    Return
    ------
    float : position of measurement
    """
    return origin + offset


def calc_ascentrate(sounding):
    """
    Calculate the ascent rate

    Input
    -----
    sounding : obj
        sounding class containing gpm
        and flight time

    Return
    ------
    sounding : obj
        sounding including the ascent rate
    """
    ascent_rate = np.diff(sounding.gpm)/(np.diff(sounding.time))
    ascent_rate = np.ma.concatenate(([0], ascent_rate))  # 0 at first measurement
    sounding.ascentrate = ascent_rate

    return sounding


def calc_temporal_resolution(sounding):
    """
    Calculate temporal resolution of sounding

    Returns the most common temporal resolution
    by calculating the temporal differences
    and returning the most common difference.

    Input
    -----
    sounding : obj
        sounding class containing flight time
        information

    Return
    ------
    temporal_resolution : float
        temporal resolution
    """
    time_differences = np.abs(np.diff(np.ma.compressed(sounding.time)))
    time_differences_counts = np.bincount(time_differences.astype(int))
    most_common_diff = np.argmax(time_differences_counts)
    temporal_resolution = most_common_diff
    return temporal_resolution


def bufr_specific_handling(sounding):
    """
    Apply bufr message specific functions

    Depending on the BUFR format, that data
    has to be prepared differently

    BUFR309053 (dropsonde)
    BUFR309056 (radiosonde descent)
    BUFR309057 (radiosonde ascent)
    - Remove last entries of time, latitude, longitude because those
        belong to the 'absoluteWindShearIn1KmLayerAbove'/
        'absoluteWindShearIn1KmLayerBelow' entries

    """
    variables = ['latitude', 'longitude', 'pressure', 'windspeed',
                 'winddirection', 'temperature', 'dewpoint', 'gpm', 'time']

    if sounding.meta_data['bufr_msg'] == 309053:
        # Nothing to do so far
        pass
    elif sounding.meta_data['bufr_msg'] == 309056:
        # Nothing to do so far
        pass
    elif sounding.meta_data['bufr_msg'] == 309057:
        # Nothing to do so far
        pass
    return sounding


def get_sounding_direction(bufr_msg):
    """
    Get direction of sounding

    1: upward
    -1: downward
    """
    if str(bufr_msg) == '309052':
        return 1
    elif str(bufr_msg) == '309053':
        return -1
    elif str(bufr_msg) == '309056':
        return -1
    elif str(bufr_msg) == '309057':
        return 1
    else:
        raise NotImplementedError('The bufr message format {} is not implemented'.format(bufr_msg))


def kelvin_to_celsius(kelvin):
    """
    Convert Kelvin to Celsius
    """
    return kelvin - 273.15


def celsius_to_kelvin(celsius):
    """
    Convert Celsius to Kelvin
    """
    return 273.15 + celsius


def pascal_to_hectoPascal(pascal):
    """
    Convert Pa to hPa
    """
    return pascal/100.


converter_dict = {'K-->C': kelvin_to_celsius,
                  'K-->degC': kelvin_to_celsius,
                  'Pa-->hPa': pascal_to_hectoPascal
                  }


def calc_relative_humidity(sounding):
    """
    Calculate relative humidity
    """
    relative_humidity = 100*(np.exp((17.625*sounding.dewpoint)/
        (243.04+sounding.dewpoint))/np.exp((17.625*sounding.temperature)/
        (243.04+sounding.temperature)))
    return relative_humidity


def convert_Tdew_to_measuredRH(sounding, manufacturer='Vaisala'):
    """
    Convert dewpoint temperatures to relative
    humidity

    This function uses the exact invers formula,
    of which is used by Vaisala MW41 and MW31 to
    convert the measured RH to dewpoint temperature
    for the BUFR output.

    Input
    -----
    sounding : sounding obj with temperature and dewpoint in Celsius

    Output
    ------
    rh_measured : array
      measured relative humidities

    NOTE: Due to numerical uncertainties (floating point), the conversion from
    dewpoint temperature back to the measured relative humidity
    might not be exact.
    """
    if manufacturer == 'Vaisala':
        dewpoint_depr = sounding.temperature - sounding.dewpoint
        temperature_K = celsius_to_kelvin(sounding.temperature)
        dewpoint_K = celsius_to_kelvin(sounding.dewpoint)
        a = 4 * (temperature_K - 273.15) * dewpoint_depr - 2 * 2711.5 * dewpoint_depr
        b = temperature_K * 30 - dewpoint_K*temperature_K - dewpoint_K * 30
        rh_measured = 100 * np.exp(-a / b)
    elif manufacturer == 'MeteoModem':
        # Following Tetens, 1930
        T = sounding.temperature
        Td = sounding.dewpoint
        e = vapor_pressure(Td)
        es = vapor_pressure(T)
        rh_measured = e/es * 100
    return rh_measured


def vapor_pressure(T, formula='Tetens1930'):
    """
    Calculate vapor pressure according to
    given formula
    """
    if formula == 'Tetens1930':
        e = 6.11 * 10**((7.5*T)/(237.3+T))
    return e


def calc_vapor_pressure(sounding):
    """
    Calculate water vapor pressure
    """
    vapor_pressure = (sounding.relativehumidity/100.) * (611.2 * np.exp((17.62 * (sounding.temperature))/(243.12 + sounding.temperature)))
    return vapor_pressure


def calc_saturation_pressure(temperature_K, method='hardy1998'):
    """
    Calculate saturation water vapor pressure
    
    Input
    -----
    temperature_K : array
        array of temperature in Kevlin or dew point temperature for actual vapor pressure
    method : str
        Formula used for calculating the saturation pressure
            'hardy1998' : ITS-90 Formulations for Vapor Pressure, Frostpoint Temperature,
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

    if method == 'hardy1998':
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
            ln_e_sw = np.sum([g[i]*temp**(i-2) for i in range(0, 7)]) + g[7]*np.log(temp)
            e_sw[t] = np.exp(ln_e_sw)
        return e_sw


def calc_wv_mixing_ratio(sounding, vapor_pressure):
    """
    Calculate water vapor mixing ratio
    """
    wv_mix_ratio = 1000.*((0.622*vapor_pressure)/(100.*sounding.pressure - vapor_pressure))
    return wv_mix_ratio


def expected_unit_check(sounding):
    """
    Check if units are as expected
    and try to convert accordingly
    """

    variables = ['displacement_lat', 'displacement_lon', 'pressure', 'windspeed',
                 'winddirection', 'temperature', 'dewpoint', 'gpm', 'time']

    expected_bufr_units = ['deg', 'deg', 'Pa', 'm/s', 'deg', 'K', 'K', 'gpm', 's']
    expected_output_units = ['deg', 'deg', 'hPa', 'm/s', 'deg', 'degC', 'degC', 'gpm', 's']

    for v, var in enumerate(variables):
        if (sounding.__dict__[var+'_unit'] != expected_output_units[v]):
            # Convert data to expected unit
            ## Find converter function
            try:
                func = converter_dict['-->'.join([sounding.__dict__[var+'_unit'], expected_output_units[v]])]
            except KeyError:
                raise UnexpectedUnit('Unit {} was expected, but got {}. Conversion was not successful'.format(expected_bufr_units[v],
                                                                           sounding.__dict__[var+'_unit']))
            else:
                sounding.__dict__[var] = func(sounding.__dict__[var])
                sounding.__dict__[var+'_unit'] = expected_output_units[v]
        else:
            pass

    return sounding


def nan_argsort(array, direction=1):
    """
    Sorting with handling nan values
    depending on the sounding direction

    Input
    -----
    array : array-like
        data to sort
    direction : int
        integer (-1, 1) to indicate on which end
        of the sorted array nan values should be
        saved.

    Result
    ------
    indices : array
        Indices that would sort the input array
    """
    tmp = array.copy().astype('float')
    tmp[np.isnan(array)] = -np.inf*direction
    return np.argsort(tmp)


def sort_sounding_by_time(sounding):
    """
    Sort sounding by altitude
    """
    sorter = nan_argsort(sounding.time, sounding.direction)
    sounding.time = sounding.time[sorter]
    sounding.ascentrate = sounding.ascentrate[sorter]
    sounding.gpm = sounding.gpm[sorter]
    sounding.pressure = sounding.pressure[sorter]
    sounding.temperature = sounding.temperature[sorter]
    sounding.relativehumidity = sounding.relativehumidity[sorter]
    sounding.dewpoint = sounding.dewpoint[sorter]
    sounding.mixingratio = sounding.mixingratio[sorter]
    sounding.windspeed = sounding.windspeed[sorter]
    sounding.winddirection = sounding.winddirection[sorter]
    sounding.latitude = sounding.latitude[sorter]
    sounding.longitude = sounding.longitude[sorter]
    sounding.extendedVerticalSoundingSignificance = sounding.extendedVerticalSoundingSignificance[sorter]

    return sounding


def exclude_1000hPa_gpm(sounding):
    """
    BUFR files include values calculated for 1000 hPa
    even when the sounding starts at an higher
    elevation.

    These values are those, where time contain
    a missing value.

    This function returns the sounding
    without the missing data in the time
    dimension.
    """
    nan_mask = ~np.isnan(sounding.time)
    sounding = exclude_sounding_level(sounding, nan_mask)

    return sounding


def exclude_sounding_level(sounding, nan_mask):
    """
    Function to exclude sounding
    """
    sounding.time = sounding.time[nan_mask]
    sounding.gpm = sounding.gpm[nan_mask]
    sounding.pressure = sounding.pressure[nan_mask]
    sounding.temperature = sounding.temperature[nan_mask]
    sounding.relativehumidity = sounding.relativehumidity[nan_mask]
    sounding.dewpoint = sounding.dewpoint[nan_mask]
    sounding.mixingratio = sounding.mixingratio[nan_mask]
    sounding.windspeed = sounding.windspeed[nan_mask]
    sounding.winddirection = sounding.winddirection[nan_mask]
    sounding.latitude = sounding.latitude[nan_mask]
    sounding.longitude = sounding.longitude[nan_mask]
    sounding.extendedVerticalSoundingSignificance = sounding.extendedVerticalSoundingSignificance[nan_mask]

    return sounding


def exclude_specific_extendedVerticalSoundingSignificance_levels(sounding, significance_bits):
    """
    Exclude levels with specific extendedVerticalSoundingSignificance

    Exclude sounding levels that contain one or more significance bits
    and no additional one.

    Input
    -----
    sounding : sounding object
        sounding

    significance_bits : array like
        significance bits that should trigger removal of sounding level

    Note: Only those levels will be excluded, where all significance bits
          that are set are also included in significance_bits.

    Example:
    exclude_specific_extendedVerticalSoundingSignificance_levels(sounding, [1,3])
    would exclude the level with the bits [1], [1,3] and [3],
    but does not exclude e.g. the levels [], [1,4], [3,5], [1,4,8,..], ....
    """
    # Get levels where extendedVerticalSoundingSignificance is not 0
    significance_levels = np.where(sounding.extendedVerticalSoundingSignificance != 0)[0]
    to_delete_mask = np.zeros(len(sounding.time), dtype=bool)
    for significance_level in significance_levels:
        current_level = sounding.extendedVerticalSoundingSignificance[significance_level]
        current_level = set(decode_extendedVerticalSoundingSignificance(current_level))
        to_delete_mask[significance_level] = ((current_level == set(significance_bits)) or (current_level == {4}))
    sounding = exclude_sounding_level(sounding, ~to_delete_mask)

    return sounding


def correct_meteomodem_surface(sounding, bufr_file):
    """
    Correct specific meteo modem soundings
    """
    d = {'ATALANTE_2020012512_1.309057.BFR': 21.676335267277707,
         'ATALANTE_2020012612_1.309057.BFR': 20.299818275809564,
         'ATALANTE_2020012612_2.309057.BFR': 20.490235079280602,
         'ATALANTE_2020012618_1.309057.BFR': 20.395027551835543,
         'ATALANTE_2020012618_2.309057.BFR': 20.72719422108181,
         'ATALANTE_2020012618_3.309057.BFR': 20.72719422108181,
         'ATALANTE_2020012618_4.309057.BFR': 20.490235079280602,
         'ATALANTE_2020012618_5.309057.BFR': 20.961170412362186,
         'ATALANTE_2020012618_6.309057.BFR': 21.096516728956811,
         'ATALANTE_2020012700_1.309057.BFR': 21.228696671077643,
         'ATALANTE_2020012700_2.309057.BFR': 21.645982812790148,
         'ATALANTE_2020012700_3.309057.BFR': 21.676335267277707,
         'ATALANTE_2020012800_1.309057.BFR': 22.734527932378281,
         'ATALANTE_2020012800_2.309057.BFR': np.nan,
         'ATALANTE_2020020218_1.309057.BFR': 20.631811617269207,
         'ATALANTE_2020020218_2.309057.BFR': 21.192243399966955,
         'ATALANTE_2020020300_1.309057.BFR': 23.250129401644031,
         'ATALANTE_2020020300_2.309057.BFR': 20.770057518368084,
         'ATALANTE_2020020300_3.309057.BFR': 21.228696671077643,
         'ATALANTE_2020020300_4.309057.BFR': 21.132797893968846,
         'ATALANTE_2020020300_5.309057.BFR': 21.594661125167445,
         'ATALANTE_2020020300_6.309057.BFR': 22.347032338755625,
         'ATALANTE_2020020300_7.309057.BFR': 21.12715294391278,
         'ATALANTE_2020020306_1.309057.BFR': 21.16565351341211,
         'ATALANTE_2020020306_2.309057.BFR': 21.350899212666299,
         'ATALANTE_2020020306_3.309057.BFR': 21.114158685504709,
         'ATALANTE_2020020306_4.309057.BFR': 21.047255593956379,
         'ATALANTE_2020020306_5.309057.BFR': 20.653281422508336,
         'ATALANTE_2020020306_6.309057.BFR': 20.927958315808677,
         'ATALANTE_2020020306_7.309057.BFR': 21.357437545602338,
         'ATALANTE_2020020312_1.309057.BFR': 21.510617824089007,
         'ATALANTE_2020020312_2.309057.BFR': 21.195165907041719,
         'ATALANTE_2020020312_3.309057.BFR': np.nan,
         'ATALANTE_2020020312_4.309057.BFR': 21.000788492386096,
         'ATALANTE_2020020312_5.309057.BFR': 21.000788492386096,
         'ATALANTE_2020020312_6.309057.BFR': 21.13786737261098,
         'ATALANTE_2020020312_7.309057.BFR': 21.877213391027734,
         'ATALANTE_2020020312_8.309057.BFR': 21.469277211341332,
         'ATALANTE_2020020318_1.309057.BFR': 21.93987011060981,
         'ATALANTE_2020020318_3.309057.BFR': 21.543164982699349,
         'ATALANTE_2020020318_4.309057.BFR': 21.431461248468246,
         'ATALANTE_2020020318_5.309057.BFR': 21.31908609493534,
         'ATALANTE_2020021100_1.309057.BFR': 21.534170247830687,
         'ATALANTE_2020021800_1.309057.BFR': 21.216973566186809,
         'ATALANTE_2020021800_2.309057.BFR': 21.483871685146806,
         'ATALANTE_2020021800_4.309057.BFR': 22.123277349708609,
         'ATALANTE_2020021800_5.309057.BFR': 21.931339923645719}
    if os.path.basename(bufr_file) in d.keys():
        Td_corrected = d[os.path.basename(bufr_file)]
        print('Correct Td from {} to {}'.format(sounding.dewpoint[0], Td_corrected))
        sounding.dewpoint[0] = Td_corrected
        if np.isnan(Td_corrected):
            sounding.temperature[0] = np.nan
    else:
        print('Unknown BFR file. No correction will be applied.')
    return sounding
