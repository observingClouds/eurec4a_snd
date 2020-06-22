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
import xarray as xr

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from config import cfg_creator as configupdater
from _helpers import *
import _thermo as thermo

json_config_fn = 'mwx_config.json'


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
        filelist = list(args['inputpath'].glob('*.bfr'))
        filelist.extend(list(args['inputpath'].glob('*.BFR')))
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
            direction_str = 'ascent'
        elif sounding.direction == -1:
            # Downward
            direction_str = 'descent'

        sounding = expected_unit_check(sounding)

        # Correct surface values of MeteoModem soundings
        if sondetype == 177:
            sounding = correct_meteomodem_surface(sounding, bufr_file)


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
        #sounding = exclude_1000hPa_gpm(sounding)

        # Ascent rate
        sounding = calc_ascentrate(sounding)

        # Sort sounding by flight time
        sounding = sort_sounding_by_time(sounding)

        # Find temporal resolution
        resolution = calc_temporal_resolution(sounding)

        # Get global attributes
        campaign = args['campaign']
        platform_id = args['platform_id']
        instrument_id = args['instrument_id']
        level = 'L1'
        if package_version_set:
            version = __version__
        else:
            version = git_module_version
        glob_attrs_dict = get_global_attrs(json_config_fn, f'{campaign}_{platform_id}_{instrument_id}_{level}')
        platform_location = glob_attrs_dict['platform_location']

        # Create outputfile with time information from file
        sounding_date = sounding.sounding_start_time
        YYYYMMDDHHMM = sounding.sounding_start_time.strftime('%Y%m%d%H%M')

        outpath_fmt = unixpath(args['outputfolder'])
        outpath = Path(sounding.sounding_start_time.strftime(outpath_fmt.as_posix()))

        if outpath.suffix == '.nc':
            outfile = Path(outpath.as_posix().format(campaign=campaign,
                                                     platform=platform_id,
                                                     level=level,
                                                     direction='{}'.format(direction_str),
                                                     date=sounding_date.strftime('%Y%m%dT%H%M'),
                                                     version=version))
        else:
            outfile = Path(os.path.join(outpath, \
                "{campaign}_{platform}_{level}-{direction}_{date}_{version}.nc".\
                format(campaign=campaign,
                       platform=platform_id,
                       level=level,
                       direction='{}'.format(direction_str),
                       date=sounding_date.strftime('%Y%m%dT%H%M'),
                       version=version)))

        if not outfile.parent.exists():
            os.makedirs(outfile.parent)

        xr_output = sounding.to_dataset()
        ## Global
        xr_output.attrs['title'] = "EUREC4A level 1 sounding data".format(campaign)
        xr_output.attrs['campaign_id'] = campaign
        xr_output.attrs['platform_id'] = f'{platform_id}'
        xr_output.attrs['instrument_id'] = f'{instrument_id}'
        xr_output.attrs['platform_location'] = platform_location
        xr_output.attrs['contact_person'] = '{name} ({mail})'.format(
            name=config['OUTPUT']['contact_person_name'],
            mail=config['OUTPUT']['contact_person_email'])
        xr_output.attrs['institution'] = config['OUTPUT']['institution']
        xr_output.attrs['converted_by'] = '{name} ({mail})'.format(
            name=config['OUTPUT']['executive_person_name'],
            mail=config['OUTPUT']['executive_person_email'])
        if 'surface_altitude' in glob_attrs_dict.keys():
            altitude = glob_attrs_dict['surface_altitude']
        else:
            altitude = config['PLATFORM']['platform_altitude']
        xr_output.attrs['surface_altitude'] = altitude
        if sondetype == 177:
            xr_output.attrs['instrument'] = f'Radiosonde GPSonde M10 by MeteoModem'
        elif sondetype == 123:
            xr_output.attrs['instrument'] = f'Radiosonde RS41-SGP by Vaisala'
        xr_output.attrs['number_of_probe'] = serial
        xr_output.attrs['sonde_type'] = str(sondetype)
        xr_output.attrs['sonde_frequency'] = sondefreq
        xr_output.attrs['date_YYYYMMDD'] = sounding_date.strftime('%Y%m%d')
        xr_output.attrs['time_of_launch_HHmmss'] = sounding_date.strftime('%H%M%S')
        xr_output.attrs['launch_unixtime'] = float(xr_output['launch_time'].values[0])
        xr_output.attrs[
            'latitude_of_launch_location'] = '{0:5.2f} deg N'.\
            format(sounding.station_lat)
        xr_output.attrs[
            'longitude_of_launch_location'] = '{0:6.2f} deg E'.\
            format(sounding.station_lon)
        xr_output.attrs['source'] = str(Path(PureWindowsPath(bufr_file)))
        xr_output.attrs['git_version'] = git_module_version
        xr_output.attrs['created_with'] = '{file} with its last modifications on {time}'. \
            format(time=time.ctime(os.path.getmtime(os.path.realpath(__file__))),
                   file=os.path.basename(__file__))
        xr_output.attrs['created_on'] = str(time.ctime(time.time()))
        xr_output.attrs['python_version'] = "{} (with numpy:{}, netCDF4:{}, eurec4a_snd:{})". \
            format(sys.version, np.__version__, netCDF4.__version__, __version__)
        xr_output.attrs['resolution'] = f'{resolution} sec'
        xr_output.attrs['Conventions'] = "CF-1.7"
        xr_output.attrs['featureType'] = "trajectory"

        # Overwrite standard attrs with those defined in config file
        # Get global meta data from mwx_config.json
        glob_attrs_dict = get_global_attrs(json_config_fn, f'{campaign}_{platform_id}_{instrument_id}_{level}')
        for attrs, value in glob_attrs_dict.items():
            xr_output.attrs[attrs] = value

        if 'extendedVerticalSoundingSignificance' in args['additional_variables']:
            xr_output['extendedVerticalSoundingSignificance'] = xr.DataArray([sounding.extendedVerticalSoundingSignificance], dims=['sounding','levels'])

        sounding_name = '{platform}__{lat:5.2f}_{lon:5.2f}__{launchtime}'.\
                        format(platform=config['PLATFORM']['platform_name_short'],
                               lat=sounding.station_lat,
                               lon=sounding.station_lon,
                               launchtime=str(YYYYMMDDHHMM))

        xr_output['sounding_id'] = xr.DataArray([sounding_name], dims = ['sounding'])

        with open(json_config_fn, 'r') as f:
            j = json.load(f)
        meta_data_dict = j['meta_data']
        for variable in xr_output.data_vars:
            if variable in meta_data_dict.keys():
                variable_meta_data = meta_data_dict[variable]
                for attr, value in variable_meta_data.items():
                    xr_output[variable].attrs[attr] = value

        # Reduce dtype to float instead of double
        xr_output.sounding_id.encoding = {'dtype': 'S1000', 'char_dim_name': 'str_dim'}
        for variable in ['altitude', 'ascentRate', 'dewPoint', 'humidity', 'latitude', 'longitude',
                         'mixingRatio', 'pressure', 'temperature', 'windDirection', 'windSpeed']:
            xr_output[variable].encoding['dtype'] = 'f4'
        xr_output['flight_time'].encoding['dtype'] = 'f8'

        for variable in xr_output.data_vars:
            xr_output[variable].encoding['zlib'] = True

        xr_output.to_netcdf(outfile, unlimited_dims=['sounding'])

        # sounding_name_parts = []
        # for char in sounding_name:
        #     sounding_name_parts.extend(char)

        # nc_prof[0, 0:len(sounding_name_parts)] = sounding_name_parts
        # nc_launchtime[0] = date2num(sounding.sounding_start_time, date_unit)
        #
        # nc_tindex[0, :] = sounding.time
        # nc_vvert[0, :] = sounding.ascentrate
        # nc_alti[0, :] = sounding.gpm
        # nc_pres[0, :] = sounding.pressure
        # nc_temp[0, :] = sounding.temperature
        # nc_rh[0, :] = sounding.relativehumidity
        # nc_dewp[0, :] = sounding.dewpoint
        # nc_mix[0, :] = sounding.mixingratio
        # nc_vhori[0, :] = sounding.windspeed
        # nc_vdir[0, :] = sounding.winddirection
        # nc_lat[0, :] = sounding.latitude
        # nc_long[0, :] = sounding.longitude
        # if 'extendedVerticalSoundingSignificance' in args['additional_variables']:
        #     nc_evss[0, :] = sounding.extendedVerticalSoundingSignificance
        # fo.close()

        logging.info('DONE: {input} converted to {output}'.format(
            input=filelist[ifile],
            output=outfile))

    time_out = time.time()
    logging.debug('System time: {}s'.format(time_out - time_in))


if __name__ == "__main__":
    main()
