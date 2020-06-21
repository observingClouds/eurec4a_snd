#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script to convert archive files (.mwx) from the Vaisala MW41
sounding system to netCDF
"""

import os
from pathlib import Path, PureWindowsPath
import subprocess as sp
import argparse
import logging
import tempfile
import glob
from xml.dom import minidom
from datetime import timedelta
import datetime as dt
import netCDF4
from netCDF4 import num2date, default_fillvals
import time #as time_module
import tqdm
import json
import numpy as np
import pandas as pd
import xarray as xr

import sys
sys.path.append('.')
from _mwx_helpers import *
from _helpers import calc_saturation_pressure, unixpath, setup_logging, load_configuration, get_global_attrs

# User-configuration
level = 'L1'
output_filename_format = '{campaign}_{platform_short}_{instrument_id}_{direction}_{level}_%Y%m%d_%H%M_{version}.nc'  # including path

json_config_fn = 'mwx_config.json'

def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-c', '--configfile', metavar="meta_information.ini", help='Provide a meta_information.ini configuration file. \n'
                                                                       'If not provided it will be searched for at:\n'
                                                                       '1. ~/meta_information.ini\n'
                                                                       '2. ../../../meta_information.ini', required=False, default=None)
    parser.add_argument("-i", "--inputfile", metavar="INPUT_FILE",
                        help="Single sonde file (mwx) or file format\n"
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

    parser.add_argument('-v', '--verbose', metavar="DEBUG",
                        help='Set the level of verbosity [DEBUG, INFO,'
                        ' WARNING, ERROR]',
                        required=False, default="INFO")
    
    parser.add_argument("--platform_name_long", metavar='Platform XYZ',
                        help="Complete name of the measurement platform e.g. Barbados Cloud Observatory",
                        default=None,
                        required=False,
                        type=str)
    parser.add_argument("--platform_name_short", metavar='XYZ',
                        help="Short name of the measurement platform e.g. BCO",
                        default=None,
                        required=False,
                        type=str)
    parser.add_argument("--platform_location", metavar='Platform location',
                        help="Location of the measurement platform e.g. Deebles Point, Barbados, West Indies",
                        default=None,
                        required=False,
                        type=str)
    parser.add_argument("--instrument_id", metavar='Instrument identifier',
                        help="Instrument identifier e.g. Vaisala-RS or Meteomodem-RS",
                        default='RS',
                        required=False,
                        type=str)
    parser.add_argument("--platform_id", metavar='Platform identifier',
                        help="Platform identifier as used in config e.g. Atalante or BCO",
                        default=None,
                        required=False,
                        type=str)
    parser.add_argument("--campaign", metavar='Campaign',
                        help="Campaign name as used in config e.g. EUREC4A",
                        default='EUREC4A',
                        required=False,
                        type=str)
    parser.add_argument("--round-like-BUFR", metavar='BOOLEAN',
                        help="Switch on/off rounding of output values as output in BUFR files",
                        default=False,
                        required=False,
                        type=bool)

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


def main(args={}):
    # Set up global configuration of BCO-MPI-GIT:
    try:
        args = get_args()
    except:
        sys.exit() 

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
            ["git", "describe", "--always", "--dirty"], stderr=sp.STDOUT).strip().decode()
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
        if (args["inputpath"] is None) and (args["inputfile"] is None):
            args["inputpath"] = config["FILES"]["INPUT_MWX"]
        if args["outputfolder"] is None:
            args["outputfolder"] = config["FILES"]["OUTPUT_MWX2NC"]

    output_folder = args['outputfolder']
    basename = os.path.basename(output_folder)
    dirname = os.path.dirname(output_folder)
    # Only filename given --> append './'
    if basename != '' and dirname == '':
        output_format = os.path.join('./', output_folder)
    # Only path given --> append standard filename format
    elif basename == '' and dirname != '':
        output_format = os.path.join(output_folder, output_filename_format)
    # Path with filename given --> nothing to change
    else:
        output_format = output_folder
    
    # Create output folder if it does not exist
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    # Loading standard config file
    with open(json_config_fn, 'r') as f:
        j=json.load(f)
    sync_sounding_values = j['sync_sounding_items']
    std_sounding_values =  j['std_sounding_items']
    radiosondes_values = j['radiosondes_sounding_items']
    meta_data_dict = j['meta_data']

    campaign = args['campaign']
    platform_id = args['platform_id']
    instrument_id = args['instrument_id']

    # Function definitions
    f_flighttime = lambda radiorx: begin_time_dt + timedelta(seconds=radiorx-float(launch_time))

    logging.debug('Create filelist')
    if args['inputfile'] is None:
        mwx_files = sorted(glob.glob(os.fspath(os.path.join(args['inputpath'],'*.mwx'))))
    else:
        mwx_files = sorted(glob.glob(os.fspath(args['inputfile'])))
    mwx_files = sorted(mwx_files)

    # Load station meta data
    if args['platform_name_long'] is None:
        platform_name_long = config['PLATFORM']['platform_name_long']
    else:
        platform_name_long = args['platform_name_long']
    if args['platform_name_short'] is None:
        platform_name_short = config['PLATFORM']['platform_name_short']
    else:
        platform_name_short = args['platform_name_short']
    if args['platform_location'] is None:
        platform_location = config['PLATFORM']['platform_location']
    else:
        platform_location = args['platform_location']

    for mwx_file in tqdm.tqdm(mwx_files):

        # Decompress/open mwx file
        with MWX(mwx_file) as mwx:
            decompressed_files = mwx.decompressed_files

            # Get the files SynchronizedSoundingData.xml, Soundings.xml, ...
            a1, sync_filename = check_availability(decompressed_files, 'SynchronizedSoundingData.xml', True)
            a2, snd_filename = check_availability(decompressed_files, 'Soundings.xml', True)
            a3, radio_filename = check_availability(decompressed_files, 'Radiosondes.xml', True)
            if np.any([not a1, not a2, not a3]):
                logging.warning('No sounding data found in {}. Skipped'.format(mwx_file))
                continue

            # Read Soundings.xml to get base time
            itemlist = read_xml(snd_filename)
            for i, item in enumerate(itemlist):
                begin_time = item.attributes['BeginTime'].value
                launch_time = item.attributes['LaunchTime'].value
                altitude = item.attributes['Altitude'].value
            begin_time_dt = dt.datetime.strptime(begin_time,'%Y-%m-%dT%H:%M:%S.%f')

            # Read sounding data
            pd_snd = get_sounding_profile(sync_filename, sync_sounding_values)

            # Read Radiosounding.xml to get sounding metadata
            sounding_meta_dict = get_sounding_metadata(radio_filename, radiosondes_values)

        # Create flight time
        pd_snd['flight_time'] = pd_snd.RadioRxTimePk.apply(f_flighttime)

        # Round sounding measurements similar to BUFR message
        if args['round_like_BUFR']:
            logging.debug('Data is rounded similar to BUFR message output')
            pd_snd_rnd = pd_snd.round(decimals={'Humidity':2, 'WindDir':0, 'Height':1, 'Temperature':2,
                                   'Latitude':5, 'Longitude':5, 'WindSpeed': 1, 'Pressure': 2, 'Altitude': 1})
        else:
            logging.debug('Data is not rounded')
            pd_snd_rnd = pd_snd

        # Split ascending and descending sounding
        direction_dict = {0: 'ascent', 1: 'descent'}
        pd_snd_asc = pd_snd_rnd.loc[pd_snd_rnd.Dropping == 0]
        pd_snd_dsc = pd_snd_rnd.loc[pd_snd_rnd.Dropping == 1]

        # Write output
        for s,sounding in enumerate([pd_snd_asc, pd_snd_dsc]):
            if len(sounding) < 2:
                logging.warning('Sounding ({}) does not contain data. Skip sounding-direction of {}'.format(direction_dict[s], mwx_file)) #, direction_dict[sounding.Dropping.values[0]]))
                continue
            xr_snd = xr.Dataset.from_dataframe(sounding)

            # Calc extra variables
            ## Ascent rate
            ascent_rate = calc_ascent_rate(xr_snd)
            ## Dew point temperature
            dewpoint = convert_RH_to_dewpoint(xr_snd.Temperature.values, xr_snd.Humidity.values)
            ## Mixing ratio
            e_s = calc_saturation_pressure(xr_snd.Temperature.values)  # calc_vapor_pressure(xr_snd, method='hardy')
            mixing_ratio = calc_wv_mixing_ratio(xr_snd, e_s)*xr_snd.Humidity.values/100.
            ## Launch time as type(datetime)
            flight_time_unix = sounding.flight_time.values.astype(float)/1e9
            launch_time_unix = flight_time_unix[0]
            launch_time_dt = num2date(launch_time_unix, "seconds since 1970-01-01")
            ## Resolution
            resolution = calc_temporal_resolution(flight_time_unix)
            ## Sounding ID
            sounding_id = '{platform}__{direction}__{lat:3.2f}_{lon:4.2f}__{time}'.format(platform = platform_id,
                                                                             direction = direction_dict[s],
                                                                             lat = xr_snd.Latitude.values[0],
                                                                             lon = xr_snd.Longitude.values[0],
                                                                             time = launch_time_dt.strftime('%Y%m%d%H%M'))
            ## Sounding type
            if sounding_meta_dict["SondeTypeName"] == 'RS41-SGP':
                sonde_type = '123'
            else:
                raise SondeTypeNotImplemented('SondeTypeName {} is not implemented'.format(sounding_meta_dict["SondeTypeName"]))

            xr_output = xr.Dataset()
            xr_output['launch_time'] = xr.DataArray([launch_time_unix], dims = ['sounding'])
            xr_output['flight_time'] = xr.DataArray([flight_time_unix], dims = ['sounding', 'levels'])
            xr_output['sounding_id'] = xr.DataArray([sounding_id], dims = ['sounding'])
            xr_output['ascentRate'] = xr.DataArray([ascent_rate], dims = ['sounding', 'levels'])
            #xr_output['altitude'] = xr.DataArray([xr_snd.GeometricHeight.values], dims = ['sounding', 'levels']) # This is different to BUFR
            xr_output['pressure'] = xr.DataArray([xr_snd.Pressure.values], dims = ['sounding', 'levels'])
            xr_output['altitude'] = xr.DataArray([xr_snd.Height.values], dims = ['sounding', 'levels'])
            xr_output['temperature'] = xr.DataArray([xr_snd.Temperature.values - 273.15], dims = ['sounding', 'levels'])
            xr_output['humidity'] = xr.DataArray([xr_snd.Humidity.values], dims = ['sounding', 'levels'])
            xr_output['dewPoint'] = xr.DataArray([dewpoint - 273.15], dims = ['sounding', 'levels'])
            xr_output['mixingRatio'] = xr.DataArray([mixing_ratio], dims = ['sounding', 'levels'])
            xr_output['windSpeed'] = xr.DataArray([xr_snd.WindSpeed.values], dims = ['sounding', 'levels'])
            xr_output['windDirection'] = xr.DataArray([xr_snd.WindDir.values], dims = ['sounding', 'levels'])
            xr_output['latitude'] = xr.DataArray([xr_snd.Latitude.values], dims = ['sounding', 'levels'])
            xr_output['longitude'] = xr.DataArray([xr_snd.Longitude.values], dims = ['sounding', 'levels'])
        #     xr_output['standard_level_flag'] = xr.DataArray()

            # Write attributes
            ## Variable
            for variable in xr_output.data_vars:
                if variable in meta_data_dict.keys():
                    variable_meta_data = meta_data_dict[variable]
                    for attr, value in variable_meta_data.items():
                        xr_output[variable].attrs[attr] = value

            ## Global
            xr_output.attrs['title'] = "Sounding data during the {} campaign (level 1)".format(campaign)
            xr_output.attrs['campaign_id'] = campaign
            xr_output.attrs['platform_id'] = f'{platform_id}'
            xr_output.attrs['instrument_id'] = f'{instrument_id}'
            xr_output.attrs['platform_location'] = platform_location
            xr_output.attrs['surface_altitude'] = '{:3.1f} m'.format(float(altitude))
            xr_output.attrs['instrument'] = f'Radiosonde {sounding_meta_dict["SondeTypeName"]} by Vaisala'
            xr_output.attrs['number_of_probe'] = sounding_meta_dict['SerialNbr']
            xr_output.attrs['sonde_type'] = sonde_type
            xr_output.attrs['sonde_frequency'] = f'{sounding_meta_dict["Frequency"]} MHz'
            xr_output.attrs['date_YYYYMMDD'] = launch_time_dt.strftime('%Y%m%d')
            xr_output.attrs['time_of_launch_HHmmss'] = launch_time_dt.strftime('%H%M%S')
            xr_output.attrs['launch_unixtime'] = launch_time_unix
            xr_output.attrs['latitude_of_launch_location'] = f'{np.round(xr_snd.Latitude.values[0],2)} deg N' # better from Soundings.xml?
            xr_output.attrs['longitude_of_launch_location'] = f'{np.round(xr_snd.Longitude.values[0],2)} deg E' # better from Soundings.xml?
            xr_output.attrs['source'] = mwx_file
            xr_output.attrs['git_version'] = git_module_version
            xr_output.attrs['created_with'] = '{file} with its last modifications on {time}'.\
                   format(time=time.ctime(os.path.getmtime(os.path.realpath(__file__))),
                   file=os.path.basename(__file__))
            xr_output.attrs['created_on'] = str(time.ctime(time.time()))
            xr_output.attrs['python_version'] = "{} (with numpy:{}, netCDF4:{}, eurec4a_snd:{})".\
                   format(sys.version, np.__version__, netCDF4.__version__, __version__)
            xr_output.attrs['resolution'] = f'{resolution} sec'
            xr_output.attrs['Conventions'] = "CF-1.7"
            xr_output.attrs['featureType'] = "trajectory"

            # Overwrite standard attrs with those defined in config file
            # Get global meta data from mwx_config.json
            glob_attrs_dict = get_global_attrs(json_config_fn, f'{campaign}_{platform_id}_{instrument_id}_{level}')
            for attrs, value in glob_attrs_dict.items():
                xr_output.attrs[attrs] = value

            # Reduce dtype to float instead of double
            xr_output.sounding_id.encoding = {'dtype': 'S1000', 'char_dim_name': 'str_dim'}
            for variable in ['altitude', 'ascentRate', 'dewPoint', 'humidity', 'latitude', 'longitude',
                             'mixingRatio', 'pressure', 'temperature', 'windDirection', 'windSpeed']:
                xr_output[variable].encoding['dtype'] = 'f4'

            for variable in xr_output.data_vars:
                xr_output[variable].encoding['zlib'] = True

            if package_version_set:
                version = __version__
            else:
                version = git_module_version
            filename = output_format.format(campaign=campaign,
                                            platform_short=platform_name_short,
                                            direction=direction_dict[sounding.Dropping.values[0]],
                                            instrument_id=args["instrument_id"],
                                            version=version,
                                            )
            filename = launch_time_dt.strftime(filename)
            xr_output.to_netcdf(filename, unlimited_dims=['sounding'])
            logging.info('File converted to {}'.format(filename))


if __name__ == "__main__":
    main()

