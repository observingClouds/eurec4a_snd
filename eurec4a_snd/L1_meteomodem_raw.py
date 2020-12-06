#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script to convert sounding files from Meteomodem soundings
 saved as .cor ASCII files to netCDF
"""

import os
from pathlib import Path
import subprocess as sp
import argparse
import logging
import glob
import datetime as dt
import netCDF4
from netCDF4 import num2date
import time
import tqdm
import json
import numpy as np
import xarray as xr

import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from _mwx_helpers import *
from _helpers import calc_saturation_pressure, unixpath, setup_logging,\
    load_configuration, get_global_attrs, read_bufr_sounding
import quality_check as qc

# User-configuration
level = 'L1'
output_filename_format = '{campaign}_{platform_short}_{instrument_id}_{level}-{direction}_%Y%m%dT%H%M_{version}.nc'  # including path

json_config_fn = '/'.join([os.path.dirname(os.path.abspath(__file__)),'config/mwx_config.json'])

def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-c', '--configfile', metavar="meta_information.ini", help='Provide a meta_information.ini configuration file. \n'
                                                                       'If not provided it will be searched for at:\n'
                                                                       '1. ~/meta_information.ini\n'
                                                                       '2. ../../../meta_information.ini', required=False, default=None)
    parser.add_argument("-i", "--inputfile", metavar="INPUT_FILE",
                        help="Single sonde file (cor) or file format\n"
                             "including wildcards",
                        default=None,
                        required=False,
                        type=unixpath)

    parser.add_argument("-b", "--BUFRfile", metavar="INPUT_FILE",
                        help="Single sonde file (BUFR) or file format\n"
                             "including wildcards for metadata only",
                        default=None,
                        required=False,
                        type=unixpath)

    parser.add_argument("-p", "--inputpath", metavar='/some/example/path',
                        help="Path to the folder containing sonde cor files.",
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
    parser.add_argument("--platform", metavar='Platform identifier',
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


def main(args=None):
    # Set up global configuration of BCO-MPI-GIT:
    if args is None:
        args = {}
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

    meta_data_dict = j['meta_data']

    campaign = args['campaign']
    platform = args['platform']
    instrument_id = args['instrument_id']

    logging.debug('Create filelist')
    if args['inputfile'] is None:
        mwx_files = sorted(glob.glob(os.fspath(os.path.join(args['inputpath'],'*.cor'))))
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
        # Read file
        pd_snd = pd.read_csv(mwx_file, delimiter='\t')

        # Get date information from filename
        date_str = mwx_file.split('_')[-2]
        date_dt = dt.datetime.strptime(date_str, '%Y%m%d%H').date()
        first_time_hour =  np.round(pd_snd.Time[0]/(60*60))
        if (first_time_hour > 12) and (dt.datetime.strptime(date_str, '%Y%m%d%H').hour == 0):
            date_dt = date_dt - dt.timedelta(days=1)  # hour is the forecast hour and therefor at midnight -1 need to be subtracted
        elif (first_time_hour < 12) and (dt.datetime.strptime(date_str, '%Y%m%d%H').hour == 0):
            date_dt = date_dt
        else:
            date_dt = date_dt


        # Create flight time
        t = lambda s: dt.datetime.combine(date_dt, dt.time(hour=0,minute=0) )+ dt.timedelta(seconds=s)
        pd_snd['flight_time'] = pd_snd.Time.apply(t)

        # Rename variables
        rename_dict = {'Time': 'Time', 'Altitude': 'Height', 'Latitude': 'Latitude', 'Longitude': 'Longitude',
                       'WindF': 'WindSpeed', 'WindD': 'WindDir', 'U': 'Humidity', 'T': 'Temperature',
                       'DP': 'Dewpoint', 'Flag': 'Flag', 'Ascent': 'AscentRate', 'Press': 'Pressure',
                       'VE': 'WindU', 'VN': 'WindV'
                       }

        pd_snd.rename(columns=rename_dict, inplace=True)

        # Convert temperature to Kelvin
        pd_snd['Temperature'] += 273.15

        # Convert radians to degree
        pd_snd['Latitude'] = np.rad2deg(pd_snd['Latitude'])
        pd_snd['Longitude'] = np.rad2deg(pd_snd['Longitude'])

        # Round sounding measurements similar to BUFR message
        if args['round_like_BUFR']:
            logging.debug('Data is rounded similar to BUFR message output')
            pd_snd_rnd = pd_snd.round(decimals={'Humidity':2, 'WindDir':0, 'Height':1, 'Temperature':2,
                                   'Latitude':5, 'Longitude':5, 'WindSpeed': 1, 'Pressure': 2, 'Altitude': 1})
        else:
            logging.debug('Data is not rounded')
            pd_snd_rnd = pd_snd

        # Quality check
        pd_snd_rnd = qc.TU_sensor(pd_snd_rnd, logging)

        # Split ascending and descending sounding
        direction_dict = {0: 'ascent', 1: 'descent'}
        idx_max_hgt = np.argmax(pd_snd_rnd.Height)
        pd_snd_asc = pd_snd_rnd.iloc[0:idx_max_hgt + 1].copy()
        pd_snd_asc['Dropping'] = 0
        pd_snd_dsc = pd_snd_rnd.iloc[idx_max_hgt + 1:].copy()
        pd_snd_dsc['Dropping'] = 1


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
            e_s = calc_saturation_pressure(xr_snd.Temperature.values)
            mixing_ratio = calc_wv_mixing_ratio(xr_snd, e_s)*xr_snd.Humidity.values/100.
            ## Launch time as type(datetime)
            flight_time_unix = sounding.flight_time.values.astype(float)/1e9
            launch_time_unix = flight_time_unix[0]
            launch_time_dt = num2date(launch_time_unix, "seconds since 1970-01-01")
            ## Resolution
            resolution = calc_temporal_resolution(flight_time_unix)
            ## Sounding ID
            sounding_id = '{platform}__{direction}__{lat:3.2f}_{lon:4.2f}__{time}'.format(platform = platform,
                                                                             direction = direction_dict[s],
                                                                             lat = xr_snd.Latitude.values[0],
                                                                             lon = xr_snd.Longitude.values[0],
                                                                             time = launch_time_dt.strftime('%Y%m%d%H%M'))

            # Read metadata from BUFR file
            ## COR files do not contain any metadata
            bufr_file = args['BUFRfile']
            if bufr_file is not None:
                meta_data_avail = True
            else:
                # Guess filename
                bufr_file = mwx_file.replace('COR','BFR').replace('.cor', '.309057.BFR')
                if os.path.isfile(bufr_file):
                    meta_data_avail = True
                else:
                    logging.warning("No sonde specific metadata is available, please consider adding paths to the"
                                    " BUFR files as well.")
                    meta_data_avail = False

            if meta_data_avail:
                sounding_bfr = read_bufr_sounding(bufr_file, logging)


            xr_output = xr.Dataset()
            xr_output['launch_time'] = xr.DataArray([sounding.flight_time.values[0]], dims = ['sounding'])
            xr_output['flight_time'] = xr.DataArray([sounding.flight_time], dims = ['sounding', 'level'])
            xr_output['sounding'] = xr.DataArray([sounding_id], dims = ['sounding'])
            xr_output['ascentRate'] = xr.DataArray([ascent_rate], dims = ['sounding', 'level'])
            xr_output['pressure'] = xr.DataArray([xr_snd.Pressure.values], dims = ['sounding', 'level'])
            xr_output['altitude'] = xr.DataArray([xr_snd.Height.values], dims = ['sounding', 'level'])
            xr_output['temperature'] = xr.DataArray([xr_snd.Temperature.values - 273.15], dims = ['sounding', 'level'])
            xr_output['humidity'] = xr.DataArray([xr_snd.Humidity.values], dims = ['sounding', 'level'])
            xr_output['dewPoint'] = xr.DataArray([dewpoint - 273.15], dims = ['sounding', 'level'])
            xr_output['mixingRatio'] = xr.DataArray([mixing_ratio], dims = ['sounding', 'level'])
            xr_output['windSpeed'] = xr.DataArray([xr_snd.WindSpeed.values], dims = ['sounding', 'level'])
            xr_output['windDirection'] = xr.DataArray([xr_snd.WindDir.values], dims = ['sounding', 'level'])
            xr_output['latitude'] = xr.DataArray([xr_snd.Latitude.values], dims = ['sounding', 'level'])
            xr_output['longitude'] = xr.DataArray([xr_snd.Longitude.values], dims = ['sounding', 'level'])
            xr_output['method'] = xr.DataArray([xr_snd.Flag.values], dims = ['sounding', 'level'])

            # Write attributes
            ## Variable
            for variable in xr_output.data_vars:
                if variable in meta_data_dict.keys():
                    variable_meta_data = meta_data_dict[variable]
                    for attr, value in variable_meta_data.items():
                        xr_output[variable].attrs[attr] = value

            ## Coordinates
            for coordinate in xr_output.coords:
                if coordinate in meta_data_dict.keys():
                    coordinate_meta_data = meta_data_dict[coordinate]
                    for attr, value in coordinate_meta_data.items():
                        xr_output[coordinate].attrs[attr] = value

            ## Global
            xr_output.attrs['title'] = "Sounding data during the {} campaign (level 1)".format(campaign)
            xr_output.attrs['campaign_id'] = campaign
            xr_output.attrs['platform'] = f'{platform}'
            xr_output.attrs['instrument_id'] = f'{instrument_id}'
            xr_output.attrs['platform_location'] = platform_location
            altitude = config['PLATFORM']['platform_altitude']
            xr_output.attrs['surface_altitude'] = altitude

            # Metadata from BUFR file
            if meta_data_avail:
                xr_output.attrs['instrument'] = f'Radiosonde GPSonde M10 by MeteoModem'
                xr_output.attrs['number_of_probe'] = sounding_bfr.meta_data['sonde_serial_number']
                xr_output.attrs['sonde_type'] = str(sounding_bfr.meta_data['radiosondeType'])
                xr_output.attrs['sonde_frequency'] = f"{sounding_bfr.meta_data['sonde_frequency']}"

            xr_output.attrs['date_YYYYMMDD'] = launch_time_dt.strftime('%Y%m%d')
            xr_output.attrs['time_of_launch_HHmmss'] = launch_time_dt.strftime('%H%M%S')
            xr_output.attrs['launch_unixtime'] = launch_time_unix
            xr_output.attrs['latitude_of_launch_location'] = f'{np.round(xr_snd.Latitude.values[0],2)} deg N' # better from Soundings.xml?
            xr_output.attrs['longitude_of_launch_location'] = f'{np.round(xr_snd.Longitude.values[0],2)} deg E' # better from Soundings.xml?
            if meta_data_avail:
                xr_output.attrs['source'] = 'data: {}, metadata: {}'.format(mwx_file, bufr_file)
            else:
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
            glob_attrs_dict = get_global_attrs(json_config_fn, f'{campaign}_{platform}_{instrument_id}_{level}')
            for attrs, value in glob_attrs_dict.items():
                xr_output.attrs[attrs] = value

            # Reduce dtype to float instead of double
            xr_output.sounding.encoding = {'dtype': 'str'}
            xr_output.method.encoding['dtype'] = 'int16'
            xr_output.method.attrs['flag_values'] = np.array([np.short(v) for v in xr_output.method.attrs['flag_values'].split(' ')])
            for variable in ['altitude', 'ascentRate', 'dewPoint', 'humidity', 'latitude', 'longitude',
                             'mixingRatio', 'pressure', 'temperature', 'windDirection', 'windSpeed']:
                xr_output[variable].encoding['dtype'] = 'f4'

            for variable in xr_output.data_vars:
                xr_output[variable].encoding['zlib'] = True

            if package_version_set:
                version = 'v{}'.format(__version__)
            else:
                version = git_module_version

            filename = output_format.format(campaign=campaign,
                                            platform_short=platform,
                                            direction=direction_dict[sounding.Dropping.values[0]],
                                            instrument_id=args["instrument_id"],
                                            version=version,
                                            level=level
                                            )
            filename = launch_time_dt.strftime(filename)

            xr_output.to_netcdf(filename, unlimited_dims=['sounding'],
                                encoding={'flight_time':{'units':"seconds since 2020-01-01", 'dtype':'float'},
                                          'launch_time':{'units':"seconds since 2020-01-01", 'dtype':'float'}
                                          })
            logging.info('File converted to {}'.format(filename))


if __name__ == "__main__":
    main()

