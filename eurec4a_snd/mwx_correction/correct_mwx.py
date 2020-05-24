#!/usr/bin/env python
# coding: utf-8

"""
Correct mwx files
- correct station paramters (e.g. barometer offset, station altitutde)
- correct offset of surface meteorology ( can be off by e.g. 30 min)
- correct launch position

INPUT
-----
- filename format of mwx files that shall be corrected
- filename of DSHIP data containing meteorological observations at the time of sounding launch
- path of corrected mwx files

EXAMPLE
-------
- python correct_mwx.py 'METEOR_*.mwx' 'EUREC4A_METEOR_DSHIP.nc' '~/corrected_mwx/'

ATTENTION: This script does not correct the mwx files. It prepares
them to be corrected. The modified mwx files still need to be
imported in the Vaisala sounding software and recalculated.
"""

import os, sys
import shutil
import tqdm
import glob
import numpy as np
import datetime as dt
from netCDF4 import num2date, date2num
import xarray as xr
import logging

from _mwx_helpers import *

## Setup
mwx_fn_fmt = sys.argv[0]  # '../EUREC4Asoundings/METEOR_*.mwx'  # File format that resolves to raw files that need correction
file_dship = sys.argv[1]  # 'EUREC4A_Meteor_DSHIP.nc'  # File containing the recorded meteorology
outpath = sys.argv[2]  # Output path for corrected mwx files

# Start logging
logging.basicConfig(filename="correct_mwx.log", level=logging.INFO)

## Open mwx file
for mwx_file in tqdm.tqdm(sorted(glob.glob(mwx_fn_fmt))):
    logging.info(mwx_file)
    # Decompress/open mwx file
    tmpdir, tmpdir_obj = getTmpDir()
    decompressed_files = np.array(decompress(mwx_file, tmpdir+'/'))
    
    ## Find files that need changes
    sync_mask = [f_sync(file) for file in decompressed_files]
    sync_filename = decompressed_files[sync_mask][0]
    snd_mask = [f_snd(file) for file in decompressed_files]
    snd_filename = decompressed_files[snd_mask][0]
    para_mask = [f_para(file) for file in decompressed_files]
    para_filename = decompressed_files[para_mask][0]
    obs_mask = [f_obs(file) for file in decompressed_files]
    obs_filename = decompressed_files[obs_mask][0]


    ## Get launch time and look up correct surface values
    # Finding surface values from DSHIP data closest to launch time

    # Read Soundings.xml to get launch time
    itemlist = read_xml(snd_filename)
    for i, item in enumerate(itemlist):
        begin_time = item.attributes['BeginTime'].value
    begin_time_dt = dt.datetime.strptime(begin_time,'%Y-%m-%dT%H:%M:%S.%f')

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
    new = '{:.2f}'.format(float(ds_sel.Tport.values)+273.15)
    logging.info('temperature: {old} --> {new}'.format(old=old,new=new))
    item.attributes['Temperature'].value = '{:.2f}'.format(float(ds_sel.Tport.values)+273.15)

    old = item.attributes['Humidity'].value
    new = '{:.0f}'.format(float(ds_sel.RHport.values))
    logging.info('humidity: {old} --> {new}'.format(old=old,new=new))
    item.attributes['Humidity'].value = '{:.0f}'.format(float(ds_sel.RHport.values))

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
    new = str(16.9)
    logging.info('Station altitude: {old} --> {new}'.format(old=old,new=new))
    item.attributes['Height'].value = new
    item.attributes['Altitude'].value = new

    old = item.attributes['AltitudeOffset'].value
    new = str(-11.5)
    logging.info('Launch site offset: {old} --> {new}'.format(old=old,new=new))
    item.attributes['AltitudeOffset'].value = new

    old = item.attributes['BarometerOffset'].value
    new = str(-16.9)
    logging.info('Barometer offset: {old} --> {new}'.format(old=old,new=new))
    item.attributes['BarometerOffset'].value = new

    old = item.attributes['GpsAntennaOffset'].value
    new = str(2.5)
    logging.info('GPS Antenna offset: {old} --> {new}'.format(old=old,new=new))
    item.attributes['GpsAntennaOffset'].value = new

    old = item.attributes['Latitude'].value
    new = str(float(ds_sel.lat.values))
    logging.info('Latitude: {old} --> {new}'.format(old=old,new=new))
    item.attributes['Latitude'].value = new

    old = item.attributes['Longitude'].value
    new = str(float(ds_sel.lon.values))
    logging.info('Longitude: {old} --> {new}'.format(old=old,new=new))
    item.attributes['Longitude'].value = new

    # Write changed Sounding.xml back to file
    xml_handle.writexml(open(snd_filename, 'w'))

    ## SoundingParameters.xml

    itemlist, xml_handle = read_xml(para_filename, return_handle=True)
    param_new_dict = {
        'Sounding.Station.Type':'2',
        'Sounding.Station.Location.Altitude':'16.9',
        'Sounding.Station.LaunchSiteOffset':'-11.5',
        'Sounding.Station.GpsAntennaOffset':'2.5',
        'Sounding.Station.BarometerOffset':'-16.9',
        'Sounding.Station.LaunchSite.Latitude':'-32768',
        'Sounding.Station.LaunchSite.Longitude':'-32768',
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
    compress(tmpdir, outpath+os.path.basename(mwx_file))
    tmpdir_obj.cleanup()

