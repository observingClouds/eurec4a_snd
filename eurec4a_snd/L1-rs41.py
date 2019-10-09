"""
Script to convert ASCII files of radiosonde type RS41

Created by: Johannes Kiliani/Lukas Frank
"""
__author__ = "Johannes Kiliani/Lukas Frank"
__date__ = "$Nov 19, 2018 05:27:06 PM$"

# insert some subroutines if possible
import time
import datetime
import calendar
import os.path
import sys
import glob
import subprocess as sp
import glob
import configparser
from configparser import ExtendedInterpolation
import argparse
import logging
import numpy as np
import netCDF4
from netCDF4 import Dataset, default_fillvals, num2date

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
    ini_path = dir_path.split("BCO")[0] + "BCO/PATH.ini"
    if not isinstance(configuration_file, str):
        if os.path.isfile("~/PATH.ini"):
            configuration_file = "~/PATH.ini"
        elif os.path.isfile(ini_path):
            configuration_file = ini_path
        if configuration_file is None or not os.path.isfile(configuration_file):
            raise FileNotFoundError(
                "No Configuration File 'PATH.ini' found. Please create one"
                " in your home directory "
                "or provide the path via the argument parsing -c.")
        else:
            logging.info("Using configuration file: %s" % configuration_file)

    conf = configparser.ConfigParser(interpolation=ExtendedInterpolation())
    conf.read(configuration_file)
    return conf


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--configfile', metavar="PATH.ini", help='Provide a PATH.ini configuration file. \n'
                                                                       'If not provided it will be searched for at:\n'
                                                                       '1. ~/PATH.ini\n'
                                                                       '2. ../../../PATH.ini', required=False)
    parser.add_argument("-i", "--inputfile", metavar="INPUT_FILE",
                        help="Single rs41 data file (txt) or file format\n"
                             "including wildcards", default=None, required=False)

    parser.add_argument("-p", "--inputpath", metavar='/some/example/path/',
                        help="Path to the folder containing radiosonde txt files",
                        default=None,
                        required=False)

    parser.add_argument("-o", "--outputfolder", metavar="/some/example/path/",
                        help="Output folder for converted files (netCDF)",
                        default=None,
                        required=False)

    parser.add_argument('-z', '--compression', metavar="COMPRESSION_LEVEL",
                        help="Set the Level of compression for the output (1-9)",
                        required=False, default=4)

    parser.add_argument('-v', '--verbose', metavar="DEBUG",
                        help='Set the level of verbosity [DEBUG, INFO,'
                        ' WARNING, ERROR]',
                        required=False, default="INFO")

    parser.add_argument('-d', '--date', metavar="YYYYMMDD", help='Provide the desired date to be processed. '
                                                                 'Fomat: YYYYMMDD', required=False, default=None)
    parsed_args = vars(parser.parse_args())

    if (parsed_args["date"] is None) and (parsed_args["inputfile"] is None):
        parser.error(
            "either --date or --inputfile needs to be set. (--date not yet implemented)")

    return parsed_args


def setup_logging(verbose):
    assert verbose in ["DEBUG", "INFO", "WARNING", "ERROR"]
    logging.basicConfig(
        level=logging.getLevelName(verbose),
        format="%(levelname)s - %(name)s - %(funcName)s - %(message)s",
        handlers=[
            logging.FileHandler(f"{__file__}.log"),
            logging.StreamHandler()
        ])


# Set up global configuration of BCO-MPI-GIT:
args = get_args()
setup_logging(args['verbose'])

try:
    config = load_configuration(args["configfile"])
except FileNotFoundError:
    if args["outputfolder"] is None and (args["inputpath"] is None and args["inputpath"] is None):
        sys.exit("No config file found! Outputfolder and Inputpath"
                 " or Inputfile need to be provided!")
    else:
        pass
else:
    if args["inputpath"] is None:
        args["inputpath"] = config["BCO_RADIOSONDES"]["PATH"]
    if args["outputfolder"] is None:
        args["outputfolder"] = config["CONVERT_BCO_RADIOSONDES"]["OUTPUT_PATH"]

try:
    git_module_version = sp.check_output(
        ["git", "describe", "--always"]).strip()
except:
    logging.info('No git-version could be found. Please consider'
                 'pulling the git repository.')
    git_module_version = "--"

compression_level = int(args["compression"])

time_in = time.time()

# Creating file list according to given arguments
if args['inputfile'] is None:
    filelist = glob.glob(args['inputpath'] + '*.txt')
else:
    filelist = glob.glob(args['inputfile'])
filelist = sorted(filelist)

logging.info('Files to process {}'.format(filelist))
for ifile in range(0, len(filelist)):
    logging.info('Reading file number {}'.format(ifile))
    linecount = 0
    with open(filelist[ifile], 'rb') as fa:
        for line in fa:  # loop over all lines within ASCII file
            logging.info(line.decode('ISO-8859-1'))
            linecount += 1  # Counts all lines within header
            string = line.split()
            if linecount == 1:
                soundid = string[-1]
            if linecount == 2:
                datest = string[-1][0:4] + string[-1][5:7] + string[-1][8:10]
                timest = string[-1][11:13] + \
                    string[-1][14:16] + string[-1][17:19]
            if linecount == 3:
                serial = string[-1]
            if linecount == 4:
                sondetype = string[-1]
            if linecount == 10:
                break

    utctime = calendar.timegm(datetime.datetime(int(datest[0:4]), int(datest[4:6]), int(datest[6:8]),
                                                int(timest[0:2]), int(timest[2:4]), int(timest[4:6])).timetuple())
    data = np.genfromtxt(
        filelist[ifile], skip_header=9, invalid_raise=False)

    num_rows = data.shape[0]  # length of raw array
    logging.info('Number of rows: {}'.format(num_rows))

    # after all needed header information is read, the reduced data field is masked for
    # NaN values and an output file produced afterward:

    tindex = np.ma.masked_invalid(data[:, 1])
    vvert_m = np.ma.masked_invalid(data[:, 9])
    alti_m = np.ma.masked_invalid(data[:, 2])
    pres_m = np.ma.masked_invalid(data[:, 3])
    temp_m = np.ma.masked_invalid(data[:, 4])
    rh_m = np.ma.masked_invalid(data[:, 5])
    # mixing ratio is not including in input file
    # but remains for compatibility with airport
    # sounding files
    mix_m = np.ma.masked_invalid(data[:, 5])
    dewp_m = np.ma.masked_invalid(data[:, 6])
    vhori_m = np.ma.masked_invalid(data[:, 8])
    vdir_m = np.ma.masked_invalid(data[:, 7])
    lat_m = np.ma.masked_invalid(data[:, 11])
    long_m = np.ma.masked_invalid(data[:, 12])

    rh_m = np.ma.masked_outside(rh_m, 0., 100.)
    pres_m = np.ma.masked_less(pres_m, 5.)
    vhori_m = np.ma.masked_greater(vhori_m, 100.)
    mix_m = np.ma.masked_greater(mix_m, -1., )

    # Find temporal resolution
    # using most common time difference
    _, indices = np.unique(np.diff(tindex), return_inverse=True)
    time_resolution = (np.diff(tindex)[np.argmax(np.bincount(indices))])

    # Create outputfile with time information from file
    YYYYMM = num2date(
        np.int(utctime), 'seconds since 1970-01-01 00:00:00 UTC').strftime('%Y%m')
    YYYYMMDDHHMM = num2date(
        np.int(utctime), 'seconds since 1970-01-01 00:00:00 UTC').strftime('%Y%m%d%H%M')
    outpath = args['outputfolder'] + YYYYMM + '/'
    if not os.path.exists(outpath):
        success = sp.call(["mkdir", "-p", outpath])
    outfile = outpath + \
        "rs__Vaisala__{}__{}check6.nc".format(time_resolution, YYYYMMDDHHMM)

    # Creation of output NetCDF file
    fo = Dataset(outfile, 'w', format='NETCDF4')

    # assign NetCDF file attributes from strings read from ASCII file header
    fo.title = 'Sounding data containing temperature, pressure, humidity,' \
               ' latitude, longitude, wind direction, wind speed, and time'
    # Platform information
    fo.platform_name = 'Barbados Cloud Observatory'
    fo.surface_altitude = '20.'
    fo.location = 'Deebles Point, Barbados, West Indies'

    # Instrument metadata
    fo.instrument = "radiosonde by Vaisala"
    fo.number_of_Probe = serial
    fo.radiosonde_type = sondetype

    # Information about launch
    fo.date_YYYYMMDD = datest
    fo.time_of_launch_HHmmss = timest
    fo.launch_unixtime = utctime
    fo.latitude_of_launch_location = '{0:5.2f} deg N'.format(np.round(lat_m[0]), 2)
    fo.longitude_of_launch_location = '{0:6.2f} deg E'.format(np.round(long_m[0]), 2)

    # Information about output
    fo.resolution = "{:g} sec".format(time_resolution)
    fo.source = filelist[ifile]
    fo.git_version = git_module_version
    fo.created_with = '{file} with its last modifications on {time}'.\
        format(time=time.ctime(os.path.getmtime(os.path.realpath(__file__))),
               file=os.path.basename(__file__))
    fo.created_on = str(time.ctime(time.time()))
    fo.python_version = "{} (with numpy:{}, netCDF4:{})".\
        format(sys.version, np.__version__, netCDF4.__version__)
    fo.Conventions = 'CF-1.7'
    fo.featureType = "trajectory"

    # Define Dimension (record length) from ASCII record counter
    fo.createDimension('levels', num_rows)
    prof_dim = fo.createDimension('trajectory', 1)
    str_dim = fo.createDimension('str_dim', 1000)
    fillval = default_fillvals['f4']

    # Creation of NetCDF Variables, including description and unit
    nc_prof = fo.createVariable(
        'trajectory', 'S1', ('trajectory', 'str_dim'), fill_value='')
    nc_prof.cf_role = "trajectory_id"
    nc_prof.long_name = 'trajectory identifier'
    nc_prof.description = 'unique string describing the trajectories origin'

    nc_launchtime = fo.createVariable('launch_time', 'f8', ('trajectory'))
    nc_launchtime.long_name = "time at which the sonde has been launched"
    nc_launchtime.units = 'seconds since 1970-01-01 00:00:00 UTC'
    nc_launchtime.calendar = 'gregorian'

    nc_tindex = fo.createVariable(
        'flight_time', 'f4', ('trajectory', 'levels'), fill_value=fillval)
    nc_tindex.long_name = 'time passed since launch'
    nc_tindex.standard_name = 'time'
    nc_tindex.units = 'seconds since launchtime'
    nc_tindex.axis = 'T'
    nc_tindex.calendar = "gregorian"
    nc_vvert = fo.createVariable(
        'ascentRate', 'f4', ('trajectory', 'levels'), fill_value=fillval)
    nc_vvert.long_name = 'ascent/descent rate of balloon or other measuring device'
    nc_vvert.description = 'ascent rate is positive/ descent rate is negative'
    nc_vvert.units = 'm/s'
    nc_vvert.coordinates = "flight_time longitude latitude pressure"
    nc_alti = fo.createVariable(
        'altitude', 'f4', ('trajectory', 'levels'), fill_value=fillval)
    nc_alti.standard_name = 'altitude'
    nc_alti.units = 'm'
    nc_alti.coordinates = "flight_time longitude latitude pressure"
    nc_pres = fo.createVariable(
        'pressure', 'f4', ('trajectory', 'levels'), fill_value=fillval)
    nc_pres.standard_name = 'air_pressure'
    nc_pres.units = 'hPa'
    nc_pres.axis = 'Z'
    nc_pres.positive = 'down'
    nc_temp = fo.createVariable(
        'temperature', 'f4', ('trajectory', 'levels'), fill_value=fillval)
    nc_temp.standard_name = 'air_temperature'
    nc_temp.units = 'degrees_Celsius'
    nc_temp.coordinates = "flight_time longitude latitude pressure"
    nc_rh = fo.createVariable(
        'humidity', 'f4', ('trajectory', 'levels'), fill_value=fillval)
    nc_rh.standard_name = 'relative_humidity'
    nc_rh.units = '%'
    nc_rh.coordinates = "flight_time longitude latitude pressure"
    nc_dewp = fo.createVariable(
        'dewPoint', 'f4', ('trajectory', 'levels'), fill_value=fillval)
    nc_dewp.standard_name = 'dew_point_temperature'
    nc_dewp.units = 'degrees_Celsius'
    nc_dewp.coordinates = "flight_time longitude latitude pressure"
    nc_mix = fo.createVariable(
        'mixingRatio', 'f4', ('trajectory', 'levels'), fill_value=fillval)
    nc_mix.long_name = 'water vapor mixing ratio'
    nc_mix.standard_name = 'humidity_mixing_ratio'
    nc_mix.units = 'g/kg'
    nc_mix.coordinates = "flight_time longitude latitude pressure"
    nc_vhori = fo.createVariable(
        'windSpeed', 'f4', ('trajectory', 'levels'), fill_value=fillval)
    nc_vhori.standard_name = 'wind_speed'
    nc_vhori.units = 'm/s'
    nc_vhori.coordinates = "flight_time longitude latitude pressure"
    nc_vdir = fo.createVariable(
        'windDirection', 'f4', ('trajectory', 'levels'), fill_value=fillval)
    nc_vdir.standard_name = 'wind_from_direction'
    nc_vdir.units = 'degrees'
    nc_vdir.coordinates = "flight_time longitude latitude pressure"
    nc_lat = fo.createVariable(
        'latitude', 'f4', ('trajectory', 'levels'), fill_value=fillval)
    nc_lat.long_name = 'latitude'
    nc_lat.standard_name = 'latitude'
    nc_lat.units = 'degrees_north'
    nc_lat.axis = 'Y'
    nc_long = fo.createVariable(
        'longitude', 'f4', ('trajectory', 'levels'), fill_value=fillval)
    nc_long.long_name = 'longitude'
    nc_long.standard_name = 'longitude'
    nc_long.units = 'degrees_east'
    nc_long.axis = 'X'

    trajectory_name = '{plattform}__{lat}_{lon}__{launchtime}'.\
                      format(plattform='BCO',
                             lat=lat_m[0],
                             lon=long_m[0],
                             launchtime=str(utctime))
    trajectory_name_parts = []
    for char in trajectory_name:
        trajectory_name_parts.extend(char)

    nc_prof[0, 0:len(trajectory_name_parts)] = trajectory_name_parts

    nc_tindex[0, :] = tindex[:]
    nc_vvert[0, :] = vvert_m[:]
    nc_alti[0, :] = alti_m[:]
    nc_pres[0, :] = pres_m[:]
    nc_temp[0, :] = temp_m[:]
    nc_rh[0, :] = rh_m[:]
    nc_dewp[0, :] = dewp_m[:]
    nc_mix[0, :] = mix_m[:]
    nc_vhori[0, :] = vhori_m[:]
    nc_vdir[0, :] = vdir_m[:]
    nc_lat[0, :] = lat_m[:]
    nc_long[0, :] = long_m[:]

    fo.close()
    logging.info('File Done')

time_out = time.time()
logging.debug('System time: ', time_out - time_in, ' s')
