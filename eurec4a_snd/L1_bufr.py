#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script to convert ASCII files of radiosonde type RS41

Original version by: Johannes Kiliani/Lukas Frank
"""

# insert some subroutines if possible
import time
import platform
from pathlib import Path, PureWindowsPath
import shutil
import os.path
import sys
import subprocess as sp
import configparser
from configparser import ExtendedInterpolation
import argparse
import logging
import numpy as np
import netCDF4
from netCDF4 import Dataset, default_fillvals, num2date, date2num

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from config import cfg_creator as configupdater
from _helpers import *
import _thermo as thermo

# ====================================================
# General MPI-BCO settings:
# ====================================================


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
    dir_path = os.path.dirname(os.path.realpath(__file__))
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


def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-c', '--configfile', metavar="meta_information.ini", help='Provide a meta_information.ini configuration file. \n'
                                                                       'If not provided it will be searched for at:\n'
                                                                       '1. ~/meta_information.ini\n'
                                                                       '2. ../../../meta_information.ini', required=False, default=None)
    parser.add_argument("-i", "--inputfile", metavar="INPUT_FILE",
                        help="Single sonde file (bufr) or file format\n"
                             "including wildcards",
                        default=None,
                        required=False,
                        type=unixpath)

    parser.add_argument("-p", "--inputpath", metavar='/some/example/path',
                        help="Path to the folder containing sonde bufr files.",
                        default=None,
                        required=False,
                        type=unixpath)

    parser.add_argument("-o", "--outputfolder", metavar="/some/example/path/",
                        help="Output folder for converted files (netCDF). You can\n"
                             " although not recommended also define an output file\n"
                             "(format). However, please share only those with the\n"
                             " the default filename.\n"
                             " The following formats can be used:\n"
                             "\t {platform}\t platform name\n"
                             "\t {location}\t platform location\n"
                             "\t {direction}\t sounding direction\n"
                             "\t {date}\t\t date of sounding release\n"
                             "\t\t or self defined date format with\n"
                             "\t\t %%Y%%m%%d %%H%%M and so on\n"
                             "\t\t and others to format the output folder dynamically.",
                        default=None,
                        required=False)

    parser.add_argument("-a", "--additional-variables",
                        help="Provide further additional variables from the BUFR file that\n"
                             "should be included in the output file.\n"
                             "currently supported are:\n"
                             "\t extendedVerticalSoundingSignificance\n",
                        default=[],
                        nargs='+',
                        required=False)

    parser.add_argument('-s', '--significant_levels', metavar=[1,...,18],
                        help='Define extendedVerticalSoundingSignificance types that should\n'
                             'not be included in the converted files.',
                        required=False, default=[], nargs='+', type=int)

    parser.add_argument('-v', '--verbose', metavar="DEBUG",
                        help='Set the level of verbosity [DEBUG, INFO,'
                        ' WARNING, ERROR]',
                        required=False, default="INFO")

    parsed_args = vars(parser.parse_args())

    if (parsed_args["inputpath"] is not None) and (parsed_args["inputfile"] is not None):
        parser.error(
            "either --inputpath or --inputfile should be used. Use --inputpath"
            "for several files and --inputfile for single ones")

    if (parsed_args["inputpath"] is None) and (parsed_args["inputfile"] is None):
        parser.error(
            "either --inputpath or --inputfile must be defined. Use --inputpath"
            "for several files and --inputfile for single ones.")

    return parsed_args


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
    logging.basicConfig(
        level=logging.getLevelName(verbose),
        format="%(levelname)s - %(name)s - %(funcName)s - %(message)s",
        handlers=[
            logging.FileHandler("{}.log".format(__file__)),
            logging.StreamHandler()
        ])


def main(args={}):
    # Set up global configuration of BCO-MPI-GIT:
    try:
        args = get_args()
    except:
        assert args!={}, "Arguments are missing"

    setup_logging(args['verbose'])

    logging.debug('Gathering version information')
    try:
        import eurec4a_snd
        __version__ = eurec4a_snd.__version__
        package_version_set = True
    except (ModuleNotFoundError, AttributeError):
        logging.debug('No eurec4a_snd package version found')
        __version__ = 'see git_version'
        package_version_set = False

    try:
        git_module_version = sp.check_output(
            ["git", "describe", "--always", "--dirty"], stderr=sp.STDOUT).strip()
        git_version_set = True
    except (sp.CalledProcessError, FileNotFoundError):
        logging.debug('No git-version could be found.')
        git_module_version = "--"
        git_version_set = False

    if (not git_version_set and not package_version_set):
        logging.warning('No version of the converter could be found!'
                        ' Please consider the installation via conda'
                        ' or if this is not working clone the git re'
                        'pository')

    logging.info('Version of script: {} (conda package), {} (git version)'.format(__version__, git_module_version))

    logging.debug('Search for configuration file and load in case present')
    try:
        config = load_configuration(args["configfile"])
    except FileNotFoundError:
        logging.info('No configuration file could be found and will now'
                     ' be created with your help')
        configupdater.update_config(os.path.abspath(os.path.dirname(__file__)) +
            '/config/meta_information_template.ini',
            Path('~/meta_information.ini').expanduser())
        sys.exit("Config file has been created at {0}. Please restart script with the option -c {0}".format(Path('~/meta_information.ini').expanduser()))
        if args["outputfolder"] is None and (args["inputpath"] is None and args["inputpath"] is None):
            sys.exit("No config file found! Outputfolder and Inputpath"
                     " or Inputfile need to be provided!")
        else:
            logging.warning("The file meta_information.ini could not be found and"
                            " no metainformation will be added to the output! "
                            " This is not recommended!")
            pass
    else:
        if args["inputpath"] is None:
            args["inputpath"] = config["FILES"]["INPUT_DAT"]
        if args["outputfolder"] is None:
            args["outputfolder"] = config["FILES"]["OUTPUT_DAT2NC"]

    time_in = time.time()
    date_unit = "seconds since 1970-01-01 00:00:00 UTC"

    logging.debug('Create filelist')
    if args['inputfile'] is None:
        filelist = args['inputpath'].glob('*.bfr')
    else:
        filelist = [args['inputfile']]
    filelist = sorted(filelist)

    logging.info('Files to process {}'.format([file.name for file in filelist]))
    for ifile, bufr_file in enumerate(filelist):
        logging.info('Reading file number {}'.format(ifile))

        json_file = convert_bufr_to_json(bufr_file, logging)
        json_flat, keys = read_json(json_file)
        if args['verbose'] != 'DEBUG':
            shutil.rmtree(os.path.dirname(json_file))

        sounding = convert_json_to_arrays(json_flat, keys)
        sounding = replace_missing_data(sounding)
        sounding = convert_list_to_array(sounding)

        sounding.latitude = calculate_coordinates(
            sounding.station_lat,
            sounding.displacement_lat)
        sounding.longitude = calculate_coordinates(
            sounding.station_lon,
            sounding.displacement_lon)

        sounding = bufr_specific_handling(sounding)

        try:
            serial = sounding.meta_data['sonde_serial_number']
        except KeyError:
            serial = '--'
        try:
            sondetype = sounding.meta_data['radiosondeType']
        except KeyError:
            sondetype = '--'
        try:
            sondefreq = sounding.meta_data['sonde_frequency']
        except KeyError:
            sondefreq = '--'

        logging.debug('Get sounding direction (ascending/descending)')
        sounding.direction = get_sounding_direction(sounding.meta_data['bufr_msg'])
        if sounding.direction == 1:
            # Upward
            direction_str = 'Ascent'
        elif sounding.direction == -1:
            # Downward
            direction_str = 'Descent'

        sounding = expected_unit_check(sounding)

        # after all needed header information is read, the reduced data field
        # is masked for NaN values and an output file produced afterward:
        logging.debug('Mask invalid sounding data')
        sounding.time = np.ma.masked_invalid(sounding.time)
        sounding.gpm = np.ma.masked_invalid(sounding.gpm)
        sounding.pressure = np.ma.masked_invalid(sounding.pressure)
        sounding.temperature = np.ma.masked_invalid(sounding.temperature)
        sounding.dewpoint = np.ma.masked_invalid(sounding.dewpoint)
        sounding.windspeed = np.ma.masked_invalid(sounding.windspeed)
        sounding.winddirection = np.ma.masked_invalid(sounding.winddirection)
        sounding.latitude = np.ma.masked_invalid(sounding.latitude)
        sounding.longitude = np.ma.masked_invalid(sounding.longitude)

        # Calculate additional variables
        e = thermo.es(sounding.dewpoint, sounding.pressure*100)
        if sondetype == 123:
            sounding.relativehumidity = convert_Tdew_to_measuredRH(sounding)
        elif sondetype == 177:
            sounding.relativehumidity = convert_Tdew_to_measuredRH(sounding, 'MeteoModem')
        else:
            logging.error('Sonde type {} is not known. Unknown how to retrieve measured RH'.format(sondetype))
            sys.exit()
        sounding.mixingratio = (thermo.Rd/thermo.Rv)*e/(sounding.pressure*100-e)*1000

        # Remove unwanted expandedVerticalSoundingSignificance levels
        sounding = exclude_specific_extendedVerticalSoundingSignificance_levels(sounding, args['significant_levels'])

        # Remove 1000hPa reduced gpm
        sounding = exclude_1000hPa_gpm(sounding)

        # Ascent rate
        sounding = calc_ascentrate(sounding)

        # Sort sounding by flight time
        sounding = sort_sounding_by_time(sounding)

        # Find temporal resolution
        time_resolution = calc_temporal_resolution(sounding)

        # Create outputfile with time information from file
        sounding_date = sounding.sounding_start_time
        YYYYMMDDHHMM = sounding.sounding_start_time.strftime('%Y%m%d%H%M')

        outpath_fmt = unixpath(args['outputfolder'])
        outpath = Path(sounding.sounding_start_time.strftime(outpath_fmt.as_posix()))

        if outpath.suffix == '.nc':
            outfile = Path(outpath.as_posix().format(platform=config['PLATFORM']['platform_name_short'],
                                                     location=config['PLATFORM']['platform_location'].
                                                               replace(' ', '').
                                                               replace(',', '').
                                                               replace(';', ''),
                                                     direction='{}Profile'.format(direction_str),
                                                     date=sounding_date.strftime('%Y%m%d_%H%M')))
        else:
            outfile = Path(os.path.join(outpath, \
                "{platform}_Sounding{direction}_{location}_{date}.nc".\
                format(platform=config['PLATFORM']['platform_name_short'],
                       location=config['PLATFORM']['platform_location'].
                                 replace(' ', '').
                                 replace(',', '').
                                 replace(';', ''),
                       direction='{}Profile'.format(direction_str),
                       date=sounding_date.strftime('%Y%m%d_%H%M'))))

        if not outfile.parent.exists():
            os.makedirs(outfile.parent)

        # Creation of output NetCDF file
        fo = Dataset(outfile, 'w', format='NETCDF4')

        # assign NetCDF file attributes from meta data
        fo.title = 'Sounding data containing temperature, pressure, humidity,'\
                   ' latitude, longitude, wind direction, wind speed, and time'
        # Platform information
        fo.platform_name = '{long} ({short})'.format(
            long=config['PLATFORM']['platform_name_long'],
            short=config['PLATFORM']['platform_name_short'])
        fo.location = config['PLATFORM']['platform_location']
        fo.surface_altitude = config['PLATFORM']['platform_altitude']

        # Instrument metadata
        fo.instrument = config['INSTRUMENT']['instrument_description']
        fo.number_of_Probe = serial
        fo.sonde_type = str(sondetype)
        fo.sonde_frequency = sondefreq

        # Information about launch
        fo.date_YYYYMMDD = sounding_date.strftime('%Y%m%d')
        fo.time_of_launch_HHmmss = sounding_date.strftime('%H%M%S')
        fo.launch_unixtime = date2num(sounding.sounding_start_time, date_unit)
        fo.latitude_of_launch_location = '{0:5.2f} deg N'.\
            format(sounding.station_lat)
        fo.longitude_of_launch_location = '{0:6.2f} deg E'.\
            format(sounding.station_lon)

        # Information about output
        fo.resolution = "{:g} sec".format(time_resolution)
        fo.source = Path(PureWindowsPath(bufr_file)).absolute().as_posix()
        fo.git_version = git_module_version
        fo.created_with = '{file} with its last modifications on {time}'.\
            format(time=time.ctime(os.path.getmtime(os.path.realpath(__file__))),
                   file=os.path.basename(__file__))
        fo.created_on = str(time.ctime(time.time()))
        fo.contact_person = '{name} ({mail})'.format(
            name=config['OUTPUT']['contact_person_name'],
            mail=config['OUTPUT']['contact_person_email'])
        fo.institution = config['OUTPUT']['institution']
        fo.converted_by = '{name} ({mail})'.format(
            name=config['OUTPUT']['executive_person_name'],
            mail=config['OUTPUT']['executive_person_email'])
        fo.python_version = "{} (with numpy:{}, netCDF4:{}, eurec4a_snd:{})".\
            format(sys.version, np.__version__, netCDF4.__version__,
                   __version__)
        fo.Conventions = 'CF-1.7'
        fo.featureType = "trajectory"

        # Define Dimension (record length) from ASCII record counter
        fo.createDimension('levels', len(sounding.pressure))
        fo.createDimension('sounding', None)
        fo.createDimension('str_dim', 1000)
        fillval = default_fillvals['f4']

        # Creation of NetCDF Variables, including description and unit
        nc_prof = fo.createVariable(
            'sounding', 'S1', ('sounding', 'str_dim'),
            fill_value='',
            zlib=True)
        nc_prof.cf_role = "sounding_id"
        nc_prof.long_name = 'sounding identifier'
        nc_prof.description = 'unique string describing the soundings origin'

        nc_launchtime = fo.createVariable('launch_time', 'f8', ('sounding'),
            zlib=True)
        nc_launchtime.long_name = "time at which the sonde has been launched"
        nc_launchtime.units = 'seconds since 1970-01-01 00:00:00 UTC'
        nc_launchtime.calendar = 'gregorian'
        nc_launchtime.standard_name = 'time'

        nc_tindex = fo.createVariable(
            'flight_time', 'f4', ('sounding', 'levels'),
            fill_value=fillval,
            zlib=True)
        nc_tindex.long_name = 'time passed since launch'
        nc_tindex.standard_name = 'time'
        nc_tindex.units = 'seconds since {launch}'.format(
            launch=sounding_date.strftime('%Y-%m-%d %H:%M:%S UTC'))
        nc_tindex.axis = 'T'
        nc_tindex.calendar = "gregorian"
        nc_vvert = fo.createVariable(
            'ascentRate', 'f4', ('sounding', 'levels'),
            fill_value=fillval,
            zlib=True)
        nc_vvert.long_name = 'ascent/descent rate of balloon or other measuring device'
        nc_vvert.description = 'ascent rate is positive/ descent rate is negative'
        nc_vvert.units = 'm/s'
        nc_vvert.coordinates = "flight_time longitude latitude pressure"
        nc_alti = fo.createVariable(
            'altitude', 'f4', ('sounding', 'levels'),
            fill_value=fillval,
            zlib=True)
        nc_alti.standard_name = 'geopotential_height'
        nc_alti.units = 'gpm'
        nc_alti.coordinates = "flight_time longitude latitude pressure"
        nc_pres = fo.createVariable(
            'pressure', 'f4', ('sounding', 'levels'),
            fill_value=fillval,
            zlib=True)
        nc_pres.standard_name = 'air_pressure'
        nc_pres.units = 'hPa'
        nc_pres.axis = 'Z'
        nc_pres.positive = 'down'
        nc_temp = fo.createVariable(
            'temperature', 'f4', ('sounding', 'levels'),
            fill_value=fillval,
            zlib=True)
        nc_temp.standard_name = 'air_temperature'
        nc_temp.units = 'degrees_Celsius'
        nc_temp.coordinates = "flight_time longitude latitude pressure"
        nc_rh = fo.createVariable(
            'humidity', 'f4', ('sounding', 'levels'),
            fill_value=fillval,
            zlib=True)
        nc_rh.standard_name = 'relative_humidity'
        nc_rh.units = '%'
        nc_rh.coordinates = "flight_time longitude latitude pressure"
        nc_dewp = fo.createVariable(
            'dewPoint', 'f4', ('sounding', 'levels'),
            fill_value=fillval,
            zlib=True)
        nc_dewp.standard_name = 'dew_point_temperature'
        nc_dewp.units = 'degrees_Celsius'
        nc_dewp.coordinates = "flight_time longitude latitude pressure"
        nc_mix = fo.createVariable(
            'mixingRatio', 'f4', ('sounding', 'levels'),
            fill_value=fillval,
            zlib=True)
        nc_mix.long_name = 'water vapor mixing ratio'
        nc_mix.standard_name = 'humidity_mixing_ratio'
        nc_mix.units = 'g/kg'
        nc_mix.coordinates = "flight_time longitude latitude pressure"
        nc_vhori = fo.createVariable(
            'windSpeed', 'f4', ('sounding', 'levels'),
            fill_value=fillval,
            zlib=True)
        nc_vhori.standard_name = 'wind_speed'
        nc_vhori.units = 'm/s'
        nc_vhori.coordinates = "flight_time longitude latitude pressure"
        nc_vdir = fo.createVariable(
            'windDirection', 'f4', ('sounding', 'levels'),
            fill_value=fillval,
            zlib=True)
        nc_vdir.standard_name = 'wind_from_direction'
        nc_vdir.units = 'degrees'
        nc_vdir.coordinates = "flight_time longitude latitude pressure"
        nc_lat = fo.createVariable(
            'latitude', 'f4', ('sounding', 'levels'),
            fill_value=fillval,
            zlib=True)
        nc_lat.long_name = 'latitude'
        nc_lat.standard_name = 'latitude'
        nc_lat.units = 'degrees_north'
        nc_lat.axis = 'Y'
        nc_long = fo.createVariable(
            'longitude', 'f4', ('sounding', 'levels'),
            fill_value=fillval,
            zlib=True)
        nc_long.long_name = 'longitude'
        nc_long.standard_name = 'longitude'
        nc_long.units = 'degrees_east'
        nc_long.axis = 'X'
        if 'extendedVerticalSoundingSignificance' in args['additional_variables']:
            nc_evss = fo.createVariable(
                'extendedVerticalSoundingSignificance', 'i4', ('sounding', 'levels'),
                fill_value=0,
                zlib=True)
            nc_evss.long_name = 'extended vertical soudning significance'
            nc_evss.description = 'see BUFR code flag table to decode'
            nc_evss.units = '-'

        sounding_name = '{platform}__{lat:5.2f}_{lon:5.2f}__{launchtime}'.\
                        format(platform=config['PLATFORM']['platform_name_short'],
                               lat=sounding.station_lat,
                               lon=sounding.station_lon,
                               launchtime=str(YYYYMMDDHHMM))
        sounding_name_parts = []
        for char in sounding_name:
            sounding_name_parts.extend(char)

        nc_prof[0, 0:len(sounding_name_parts)] = sounding_name_parts
        nc_launchtime[0] = date2num(sounding.sounding_start_time, date_unit)

        nc_tindex[0, :] = sounding.time
        nc_vvert[0, :] = sounding.ascentrate
        nc_alti[0, :] = sounding.gpm
        nc_pres[0, :] = sounding.pressure
        nc_temp[0, :] = sounding.temperature
        nc_rh[0, :] = sounding.relativehumidity
        nc_dewp[0, :] = sounding.dewpoint
        nc_mix[0, :] = sounding.mixingratio
        nc_vhori[0, :] = sounding.windspeed
        nc_vdir[0, :] = sounding.winddirection
        nc_lat[0, :] = sounding.latitude
        nc_long[0, :] = sounding.longitude
        if 'extendedVerticalSoundingSignificance' in args['additional_variables']:
            nc_evss[0, :] = sounding.extendedVerticalSoundingSignificance
        fo.close()
        logging.info('DONE: {input} converted to {output}'.format(
            input=filelist[ifile],
            output=outfile))

    time_out = time.time()
    logging.debug('System time: {}s'.format(time_out - time_in))


if __name__ == "__main__":
    main()
