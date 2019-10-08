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
import numpy as np
from netCDF4 import Dataset, default_fillvals, num2date
import subprocess as sp
import glob
import configparser
from configparser import ExtendedInterpolation
import argparse


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

        if not os.path.isfile(configuration_file):
            raise FileNotFoundError(
                "No Configuration File 'PATH.ini' found. Please create one in your home directory "
                "or provide the path via the argument parsing -c.")
        else:
            print("Using configuration file: %s" % configuration_file)

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

    parser.add_argument('-d', '--date', metavar="YYYYMMDD", help='Provide the desired date to be processed. '
                                                                 'Fomat: YYYYMMDD', required=False, default=None)
    parsed_args = vars(parser.parse_args())

    if (parsed_args["date"] is None) and (parsed_args["inputfile"] is None):
        parser.error(
            "either --date or --inputfile needs to be set. (--date not yet implemented)")

    return parsed_args


# Set up global configuration of BCO-MPI-GIT:
args = get_args()

try:
    config = load_configuration(args["configfile"])
except FileNotFoundError:
    print("No config file found! Outputfolder and Inputpath or Inputfile need to be provided!")

if args["inputpath"] is None:
    args["inputpath"] = config["BCO_RADIOSONDES"]["PATH"]
if args["outputfolder"] is None:
    args["outputfolder"] = config["CONVERT_BCO_RADIOSONDES"]["OUTPUT_PATH"]

try:
    git_module_version = sp.check_output(
        ["git", "describe", "--always"]).strip()
except:
    git_module_version = "--"

compression_level = int(args["compression"])

time_in = time.time()

# <<<<<<< L1-rs41.py
# Creating file list according to given arguments
if args['inputfile'] is None:
    filelist = glob.glob(args['inputpath'] + '*.txt')
else:
    filelist = glob.glob(args['inputfile'])
filelist = sorted(filelist)
# =======
inroot = '/media/lukas/Lidar_LF/Radiosondes/ASCII/201902/'
outroot = '/media/lukas/Lidar_LF/Radiosondes/NetCDF/'

filelist0 = os.listdir(inroot)
filelist0.sort()
# >>>>>>> L1-rs41.py

print('Files to process {}'.format(filelist))
for ifile in range(0, len(filelist)):
    print('Reading file number {}'.format(ifile))
    linecount = 0
    with open(filelist[ifile], 'rb') as fa:
        for line in fa:  # loop over all lines within ASCII file
            print(line)
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
        filelist[ifile], skip_header=9, invalid_raise=False, encoding="ISO-8859-1")

    num_rows = data.shape[0]  # length of raw array
    print('Number of rows: {}'.format(num_rows))

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
        "rs__Vaisala__{}check6.nc".format(time_resolution, YYYYMMDDHHMM)

    # Creation of output NetCDF file
    fo = Dataset(outfile, 'w', format='NETCDF4')

    # assign NetCDF file attributes from strings read from ASCII file header
    fo.title = 'Sounding data containing temperature, pressure, humidity,' \
               ' latitude, longitude, wind direction, wind speed, and time'
    fo.featureType = "profile"
    fo.date_YYYYMMDD = datest
    fo.time_of_launch_HHmmss = timest
    fo.Launch_unixtime = utctime
    fo.location = 'Deebles Point, Barbados, West Indies'
    fo.Latitude_of_launch_location = '13.16 deg N'
    fo.Longitude_of_launch_location = '59.43 deg W'
    fo.Altitude_of_launch_location_meters = '20.'
    fo.instrument = "radiosonde by Vaisala"
    fo.Number_of_Probe = serial
    fo.Radiosonde_type = sondetype
    fo.resolution = "{:g} sec".format(time_resolution)
    fo.original_file = filelist[ifile]
    fo.git_version = git_module_version

    # add file history, including file creation timestamp, path of python file and of ASCII data file
    fo.history = 'File created ' + time.ctime(time.time()) + ' with ' + os.path.basename(__file__), \
                 ', last program modification on ' + \
        time.ctime(os.path.getmtime(os.path.realpath(__file__)))

    # Define Dimension (record length) from ASCII record counter
    fo.createDimension('levels', num_rows)
    fillval = default_fillvals['f4']

    # Creation of NetCDF Variables, including description and unit
    nc_tindex = fo.createVariable('time', 'f4', ('levels'), fill_value=fillval)
    nc_tindex.long_name = 'time passed since launch'
    nc_tindex.units = 's'
    nc_vvert = fo.createVariable(
        'ascentRate', 'f4', ('levels'), fill_value=fillval)
    nc_vvert.long_name = 'Ascent/descent rate (vertical velocity)'
    nc_vvert.units = 'm/s'
    nc_alti = fo.createVariable(
        'altitude', 'f4', ('levels'), fill_value=fillval)
    nc_alti.long_name = 'altitude'
    nc_alti.units = 'm'
    nc_pres = fo.createVariable(
        'pressure', 'f4', ('levels'), fill_value=fillval)
    nc_pres.long_name = 'pressure'
    nc_pres.units = 'hPa'
    nc_temp = fo.createVariable(
        'temperature', 'f4', ('levels'), fill_value=fillval)
    nc_temp.long_name = 'air temperature'
    nc_temp.units = 'degrees Celsius'
    nc_rh = fo.createVariable('humidity', 'f4', ('levels'), fill_value=fillval)
    nc_rh.long_name = 'relative humidity'
    nc_rh.units = '%'
    nc_dewp = fo.createVariable(
        'dewPoint', 'f4', ('levels'), fill_value=fillval)
    nc_dewp.long_name = 'dew point'
    nc_dewp.units = 'degrees Celsius'
    nc_mix = fo.createVariable(
        'mixingRatio', 'f4', ('levels'), fill_value=fillval)
    nc_mix.long_name = 'Water vapor mixing ratio'
    nc_mix.units = 'g/kg'
    nc_vhori = fo.createVariable(
        'windSpeed', 'f4', ('levels'), fill_value=fillval)
    nc_vhori.long_name = 'wind speed)'
    nc_vhori.units = 'm/s'
    nc_vdir = fo.createVariable(
        'windDirection', 'f4', ('levels'), fill_value=fillval)
    nc_vdir.long_name = 'wind direction'
    nc_vdir.units = 'degrees'
    nc_lat = fo.createVariable(
        'latitude', 'f4', ('levels'), fill_value=fillval)
    nc_lat.long_name = 'latitude'
    nc_lat.units = 'degrees-north'
    nc_long = fo.createVariable(
        'longitude', 'f4', ('levels'), fill_value=fillval)
    nc_long.long_name = 'longitude'
    nc_long.units = 'degrees-east'

    nc_tindex[:] = tindex[:]
    nc_vvert[:] = vvert_m[:]
    nc_alti[:] = alti_m[:]
    nc_pres[:] = pres_m[:]
    nc_temp[:] = temp_m[:]
    nc_rh[:] = rh_m[:]
    nc_dewp[:] = dewp_m[:]
    nc_mix[:] = mix_m[:]
    nc_vhori[:] = vhori_m[:]
    nc_vdir[:] = vdir_m[:]
    nc_lat[:] = lat_m[:]
    nc_long[:] = long_m[:]

    fo.close()
    print('File Done')

time_out = time.time()
print('System time: ', time_out - time_in, ' s')
