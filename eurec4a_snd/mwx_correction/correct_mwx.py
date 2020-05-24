#!/usr/bin/env python
# coding: utf-8

"""
Correct mwx files
- correct station paramters (e.g. barometer offset, station altitutde)
- correct offset of surface meteorology ( can be off by e.g. 30 min)
- correct launch position

ATTENTION: This script does not correct the mwx files. It prepares
them to be corrected. The modified mwx files still need to be
imported in the Vaisala sounding software and recalculated.
"""

import os
import shutil
import tqdm
import glob
import numpy as np
import datetime as dt
from netCDF4 import num2date, date2num
import xarray as xr
import logging
import argparse

from _mwx_helpers import *

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inputfilefmt', required=True,
                    help='Input filename format (mwx files), e.g. METEOR_*.mwx')
parser.add_argument('-o', '--outputfolder', required=False, default='./',
                    help='Output folder of corrected mwx files')
parser.add_argument('-m', '--meteorology', required=True,
                    help='DSHIP data file containing position and meteorology')
args = vars(parser.parse_args())

## Setup
mwx_fn_fmt = args['inputfilefmt']  # '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v3/level0/MET/MWX/METEOR_*.mwx'  # File format that resolves to raw files that need correction
file_dship = args['meteorology']  # 'EUREC4A_Meteor_DSHIP.nc'  # File containing the recorded meteorology
output_dir = args['outputfolder']  # '/scratch/m/m300408/corrected_mwx/'

setup_dict = {'nonDWD1':
		{'StationAltitude':16.9,
		 'GPSAntOffset':2.5,
		 'LaunchSiteOffset':-11.5,
		 'BarometerOffset':-16.9
		},
              'nonDWD2':
	 	{'StationAltitude':16.9,
		 'GPSAntOffset':2.5,
		 'LaunchSiteOffset':-14.2,
		 'BarometerOffset':-16.9
		},
	      'DWD1':
		{'StationAltitude':5.4,
		 'GPSAntOffset':3,
		 'LaunchSiteOffset':0.0,
		 'BarometerOffset':-5.4
		},
	      'DWD2':
		{'StationAltitude':5.4,
		 'GPSAntOffset':3,
		 'LaunchSiteOffset':-2.7,
		 'BarometerOffset':-5.4
		}
	     }

# Start logging
logging.basicConfig(filename="correct_mwx.log", level=logging.INFO)

## Open mwx file
for mwx_file in tqdm.tqdm(sorted(glob.glob(mwx_fn_fmt))):
    logging.info(mwx_file)
    
    # Get operator
    if 'DWD' in mwx_file:
        operator = 'DWD'
    else:
        operator = 'nonDWD'

    # Decompress/open mwx file
    tmpdir, tmpdir_obj = getTmpDir()
    decompressed_files = np.array(decompress(mwx_file, tmpdir+'/'))
    
    ## Find files that need changes
    try:
        snd_mask = [f_snd(file) for file in decompressed_files]
        snd_filename = decompressed_files[snd_mask][0]
        para_mask = [f_para(file) for file in decompressed_files]
        para_filename = decompressed_files[para_mask][0]
        obs_mask = [f_obs(file) for file in decompressed_files]
        obs_filename = decompressed_files[obs_mask][0]
    except:
        logging.error('File not found in {}'.format(decompressed_files))
        continue


    ## Get launch time and look up correct surface values
    # Finding surface values from DSHIP data closest to launch time

    # Read Soundings.xml to get launch time
    itemlist = read_xml(snd_filename)
    for i, item in enumerate(itemlist):
        begin_time = item.attributes['BeginTime'].value
    begin_time_dt = dt.datetime.strptime(begin_time,'%Y-%m-%dT%H:%M:%S.%f')

    # Decide which corrections apply
    if (begin_time_dt < dt.datetime(2020,2,9,18,0)) or (begin_time_dt > dt.datetime(2020,2,20,0,0)):
        mode = '1'
    else:
        mode = '2'

    # Read DSHIP data
    ds_DSHIP = xr.open_dataset(file_dship, decode_times=False)
    ds_DSHIP = ds_DSHIP.sortby('time')

    ds_sel = ds_DSHIP.sel(time=date2num(begin_time_dt, "seconds since 1970-01-01"), method='nearest')
    logging.info('For launch time {}, the closest DSHIP time found is {}'.format(begin_time_dt.strftime('%Y%m%d %H:%M'), num2date(ds_sel.time.values, "seconds since 1970-01-01 00:00:00 UTC").strftime('%Y%m%d %H:%M')))

    ### Surface observations
    itemlist, xml_handle = read_xml(obs_filename, return_handle=True)
    xml_dship_dict = {'Pressure': 'p',
                      'LaunchSitePressure': 'p',
                      'Temperature': 'Tport',
                      'Humidity': 'RHport',
                      'WindDirection': 'DD_true',
                      'WindSpeed': 'FF_true'
                     }

    item = itemlist[0]
    old = item.attributes['Pressure'].value
    new = '{:.2f}'.format(float(ds_sel.p.values))
    logging.info('pressure: {old} --> {new}'.format(old=old,new=new))
    item.attributes['Pressure'].value = '{:.2f}'.format(float(ds_sel.p.values))
    item.attributes['LaunchSitePressure'].value = '{:.2f}'.format(float(ds_sel.p.values))

    old = item.attributes['Temperature'].value
    new = '{:.2f}'.format(np.nanmean([float(ds_sel.Tport.values),float(ds_sel.Tstar.values)])+273.15)
    logging.info('temperature: {old} --> {new}'.format(old=old,new=new))
    item.attributes['Temperature'].value = new

    old = item.attributes['Humidity'].value
    new = '{:.0f}'.format(np.nanmean([float(ds_sel.RHport.values),float(ds_sel.RHstar.values)]))
    logging.info('humidity: {old} --> {new}'.format(old=old,new=new))
    item.attributes['Humidity'].value = new

    old = item.attributes['WindDirection'].value
    new = '{:.0f}'.format(float(ds_sel.DD_true.values))
    logging.info('winddirection: {old} --> {new}'.format(old=old,new=new))
    item.attributes['WindDirection'].value = '{:.0f}'.format(float(ds_sel.DD_true.values))

    old = item.attributes['WindSpeed'].value
    new = '{:.1f}'.format(float(ds_sel.FF_true.values))
    logging.info('windspeed: {old} --> {new}'.format(old=old,new=new))
    item.attributes['WindSpeed'].value = '{:.1f}'.format(float(ds_sel.FF_true.values))

    # Write changed SurfaceObservations.xml back to file
    xml_handle.writexml(open(obs_filename, 'w'))


    # ### Station parameters
    # Sounding.xml
    itemlist, xml_handle = read_xml(snd_filename, return_handle=True)

    item = itemlist[0]
    old = item.attributes['Height'].value
    new = str(setup_dict[operator+mode]['StationAltitude'])
    logging.info('Station altitude: {old} --> {new}'.format(old=old,new=new))
    item.attributes['Height'].value = new
    item.attributes['Altitude'].value = new

    old = item.attributes['AltitudeOffset'].value
    new = str(setup_dict[operator+mode]['LaunchSiteOffset'])
    logging.info('Launch site offset: {old} --> {new}'.format(old=old,new=new))
    item.attributes['AltitudeOffset'].value = new

    old = item.attributes['BarometerOffset'].value
    new = str(setup_dict[operator+mode]['BarometerOffset'])
    logging.info('Barometer offset: {old} --> {new}'.format(old=old,new=new))
    item.attributes['BarometerOffset'].value = new

    old = item.attributes['GpsAntennaOffset'].value
    new = str(setup_dict[operator+mode]['GPSAntOffset'])
    logging.info('GPS Antenna offset: {old} --> {new}'.format(old=old,new=new))
    item.attributes['GpsAntennaOffset'].value = new

#    old = item.attributes['Latitude'].value
#    new = str(float(ds_sel.lat.values))
#    logging.info('Latitude: {old} --> {new}'.format(old=old,new=new))
#    item.attributes['Latitude'].value = new

#    old = item.attributes['Longitude'].value
#    new = str(float(ds_sel.lon.values))
#    logging.info('Longitude: {old} --> {new}'.format(old=old,new=new))
#    item.attributes['Longitude'].value = new

    # Write changed Sounding.xml back to file
    xml_handle.writexml(open(snd_filename, 'w'))

    ## SoundingParameters.xml

    itemlist, xml_handle = read_xml(para_filename, return_handle=True)
    param_new_dict = {
        'Sounding.Station.Type':'2',
        'Sounding.Station.Location.Altitude':str(setup_dict[operator+mode]['StationAltitude']),
        'Sounding.Station.LaunchSiteOffset':str(setup_dict[operator+mode]['LaunchSiteOffset']),
        'Sounding.Station.GpsAntennaOffset':str(setup_dict[operator+mode]['GPSAntOffset']),
        'Sounding.Station.BarometerOffset':str(setup_dict[operator+mode]['BarometerOffset']),
        'Sounding.Station.LaunchSite.Latitude':str(float(ds_sel.lat.values)),
        'Sounding.Station.LaunchSite.Longitude':str(float(ds_sel.lon.values)),
        'Sounding.Station.Location.Latitude':'-32768',
        'Sounding.Station.Location.Longitude':'-32768'
    }

    for i, item in enumerate(itemlist):
        for var in ['Parameter']:
            key = item.attributes['ParameterNamePk'].value
            if key in param_new_dict.keys():
                old = item.attributes['ParameterValue'].value
                new = param_new_dict[key]
                logging.info('{key} {old} --> {new}'.format(key=key, old=old, new=new))
                item.attributes['ParameterValue'].value = new

    # Write changed SurfaceObservations.xml back to file
    xml_handle.writexml(open(para_filename, 'w'))


    ## Compress files again to mwx file
    compress(tmpdir, output_dir+os.path.basename(mwx_file))
    tmpdir_obj.cleanup()

