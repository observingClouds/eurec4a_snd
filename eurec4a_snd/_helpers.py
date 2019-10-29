"""
Helper functions for bufr import
"""
import tempfile
import os
import numpy as np
import json
import datetime as dt

class UnitChangedError(Exception):
    pass

class UnexpectedUnit(Exception):
    pass

def convert_bufr_to_json(bufr_fn):
    tmp_folder = tempfile.mkdtemp()
    tmp_output_json_fn = tmp_folder+'/tmp_bufr.json'
    r = os.system("bufr_dump -j s {} > {}".format(bufr_fn, tmp_output_json_fn))
    print("Converted {} to {}".format(bufr_fn, tmp_output_json_fn))
    return tmp_output_json_fn

def flatten_json(y):
    global r
    out = {}
    r=0
    def flatten(x, name='', number=0):
        global r
        if type(x) is dict:
            for a in x:
                flatten(x[a], a)
        elif type(x) is list:
            i = 0
            for a in x:
                number=flatten(a)
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
    with open(json_fn) as f:
        p=json.load(f)

    bfr_json_flat = flatten_json(p)
    keys = bfr_json_flat.keys()
    key_keys = []
    for key in keys:
        if 'key' in key:
            key_keys.append(key[:-4])
    
    return bfr_json_flat, key_keys

def calculate_sounding_start(sounding):
    import datetime as dt
    sounding.sounding_start_time = dt.datetime(year,
                                               month,
                                               day,
                                               hour,
                                               minute,
                                               second)
    return sounding

def convert_json_to_arrays(json_flat, key_keys):
    """
    Convert json data to array
    """
    
    class Sounding:
        def __init__(self):
            self.station_lat = None
            self.station_lon = None
            self.sounding_start_time = None
            self.time  = []
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
            self.meta_data = {}
    
    s = Sounding()
            
    year_ = None
    month_ = None
    day_ = None
    hour_ = None
    minute_ = None
    seconds_ = None

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
            elif (s.pressure_unit != json_flat[key_key+'_units']):
                raise UnitChangedError('{} and {} are not same unit'.format(s.pressure_unit,
                                                                            json_flat[key_key+'_units']))
        elif json_flat[key_key+'_key'] == 'windSpeed':
            s.windspeed.append(json_flat[key_key+'_value'])
            if s.windspeed_unit is None:
                s.windspeed_unit = json_flat[key_key+'_units']
            # Unit consistency test
            elif (s.windspeed_unit != json_flat[key_key+'_units']):
                raise UnitChangedError('{} and {} are not same unit'.format(s.windspeed_unit,
                                                                            json_flat[key_key+'_units']))
        elif json_flat[key_key+'_key'] == 'windDirection':
            s.winddirection.append(json_flat[key_key+'_value'])
            if s.winddirection_unit is None:
                s.winddirection_unit = json_flat[key_key+'_units']
            # Unit consistency test
            elif (s.winddirection_unit != json_flat[key_key+'_units']):
                raise UnitChangedError('{} and {} are not same unit'.format(s.winddirection_unit,
                                                                            json_flat[key_key+'_units']))
        elif json_flat[key_key+'_key'] == 'nonCoordinateGeopotentialHeight':
            s.gpm.append(json_flat[key_key+'_value'])
            if s.gpm_unit is None:
                s.gpm_unit = json_flat[key_key+'_units']
            # Unit consistency test
            elif (s.gpm_unit != json_flat[key_key+'_units']):
                raise UnitChangedError('{} and {} are not same unit'.format(s.gpm_unit,
                                                                            json_flat[key_key+'_units']))
        elif json_flat[key_key+'_key'] == 'airTemperature':
            s.temperature.append(json_flat[key_key+'_value'])
            if s.temperature_unit is None:
                s.temperature_unit = json_flat[key_key+'_units']
            # Unit consistency test
            elif (s.temperature_unit != json_flat[key_key+'_units']):
                raise UnitChangedError('{} and {} are not same unit'.format(s.temperature_unit,
                                                                            json_flat[key_key+'_units']))
        elif json_flat[key_key+'_key'] == 'dewpointTemperature':
            s.dewpoint.append(json_flat[key_key+'_value'])
            if s.dewpoint_unit is None:
                s.dewpoint_unit = json_flat[key_key+'_units']
            # Unit consistency test
            elif (s.dewpoint_unit != json_flat[key_key+'_units']):
                raise UnitChangedError('{} and {} are not same unit'.format(s.dewpoint_unit,
                                                                            json_flat[key_key+'_units']))
        elif json_flat[key_key+'_key'] == 'latitudeDisplacement':
            s.displacement_lat.append(json_flat[key_key+'_value'])
            if s.displacement_lat_unit is None:
                s.displacement_lat_unit = json_flat[key_key+'_units']
            # Unit consistency test
            elif (s.displacement_lat_unit != json_flat[key_key+'_units']):
                raise UnitChangedError('{} and {} are not same unit'.format(s.displacement_lat_unit,
                                                                            json_flat[key_key+'_units']))
        elif json_flat[key_key+'_key'] == 'longitudeDisplacement':
            s.displacement_lon.append(json_flat[key_key+'_value'])
            if s.displacement_lon_unit is None:
                s.displacement_lon_unit = json_flat[key_key+'_units']
            # Unit consistency test
            elif (s.displacement_lon_unit != json_flat[key_key+'_units']):
                raise UnitChangedError('{} and {} are not same unit'.format(s.displacement_lon_unit,
                                                                            json_flat[key_key+'_units']))
        elif json_flat[key_key+'_key'] == 'timePeriod':
            s.time.append(json_flat[key_key+'_value'])
            if s.time_unit is None:
                s.time_unit = json_flat[key_key+'_units']
            # Unit consistency test
            elif (s.time_unit != json_flat[key_key+'_units']):
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
                # Error probably caused, because there are several unexpandedDescriptors
                # which seems to be only the case for dropsondes?!
                s.meta_data['bufr_msg'] = 309053
                
    
    s.sounding_start_time = dt.datetime(year,
                                        month,
                                        day,
                                        hour,
                                        minute,
                                        second)
    
    return s

def replace_missing_data(sounding):
    # Remove Nones from lists and
    def replace_none(entry):
        if entry == None:
            return np.nan
        else:
            return entry

    variables = ['displacement_lat', 'displacement_lon', 'pressure', 'windspeed',
                 'winddirection', 'temperature', 'dewpoint', 'gpm', 'time']
    
    for var in variables:
        sounding.__dict__[var] =list(map(replace_none, sounding.__dict__[var]))
    
    return sounding

def convert_list_to_array(sounding):
    """
    Convert datatype of sounding
    """
    variables = ['displacement_lat', 'displacement_lon', 'pressure', 'windspeed',
                 'winddirection', 'temperature', 'dewpoint', 'gpm', 'time']
    
    for var in variables:
        sounding.__dict__[var] = np.array(sounding.__dict__[var])
    
    return sounding

def calculate_coordinates(origin, offset):
    return origin + offset

def bufr_specific_handling(sounding):
    """
    Apply bufr message specific functions
    
    Depending on the BUFR format, that data
    has to be prepared differently
    
    BUFR309053 (dropsonde)
    BUFR309056 (radiosonde descent)
    BUFR309057 (radiosonde ascent)
    - Remove last entries of time, latitude, longitude because those
        belong to the 'absoluteWindShearIn1KmLayerAbove'/ 'absoluteWindShearIn1KmLayerBelow' entries
    
    """
    variables = ['latitude', 'longitude', 'pressure', 'windspeed',
                 'winddirection', 'temperature', 'dewpoint', 'gpm', 'time']
    
    if sounding.meta_data['bufr_msg'] == 309053:
        # Nothing to do so far
        pass
    elif sounding.meta_data['bufr_msg'] == 309056:
        if np.isnan(sounding.time[0]):
            for var in variables:
                sounding.__dict__[var] = sounding.__dict__[var][1:] 
        sounding.latitude = sounding.latitude[:-1]
        sounding.longitude = sounding.longitude[:-1]
        sounding.pressure = sounding.pressure[:-1]
        sounding.time = sounding.time[:-1]
    elif sounding.meta_data['bufr_msg'] == 309057:
        sounding.latitude = sounding.latitude[:-1]
        sounding.longitude = sounding.longitude[:-1]
        sounding.pressure = sounding.pressure[:-1]
        sounding.time = sounding.time[:-1]
    return sounding

def get_sounding_direction(bufr_msg):
    """
    Get direction of sounding
    
    1: upward
    -1: downward
    """
    
    if str(bufr_msg) == '309053':
        return -1
    elif str(bufr_msg) == '309056':
        return -1
    elif str(bufr_msg) == '309057':
        return 1
    else:
        raise NotImplementedError('The bufr message format {} is not implemented'.format(bufr_msg))

def kelvin_to_celsius(kelvin):
    return kelvin - 273.15

def pascal_to_hectoPascal(pascal):
    return pascal/100.

converter_dict = {'K-->C': kelvin_to_celsius,
                  'K-->degC': kelvin_to_celsius,
                  'Pa-->hPa': pascal_to_hectoPascal
                 }

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
    
