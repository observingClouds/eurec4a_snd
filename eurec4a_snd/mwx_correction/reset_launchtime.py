"""
Script to reset the launch time in a simulated sounding.

Corrections to the station parameters e.g. station altitude,
barometer offset etc. can make it necessary to simulate a
previous sounding as a simple recalculation does not take
changes of the station parameters into account.

Because simulated soundings are treated by the Vaisala MW41
software as new soundings, they loose there original times.

These times are restored to the best possible degree by
executing this script.

USAGE
-----
python reset_launchtime.py simulated*.mwx originial*.mwx output_path
"""

import numpy as np
import os
import glob
import sys
sys.path.append('./')
from _mwx_helpers import *
import datetime as dt
import tqdm

sim_mwx_fmt = sys.argv[0]
orig_mwx_fmt = sys.argv[1]
output_path = sys.argv[2]

# Get time information from original mwx files
sounding_date_serial_dict = {}
for mwx_file in tqdm.tqdm(sorted(glob.glob(orig_mwx_fmt))):
    tmpdir, tmpdir_obj = getTmpDir()
    decompressed_files = np.array(decompress(mwx_file, tmpdir+'/'))
    
    ## Find files that need changes
    snd_mask = [f_snd(file) for file in decompressed_files]
    snd_filename = decompressed_files[snd_mask][0]
    radio_mask = [f_radio(file) for file in decompressed_files]
    radio_filename = decompressed_files[radio_mask][0]

    # Read Soundings.xml to get launch time
    itemlist = read_xml(snd_filename)
    for i, item in enumerate(itemlist):
        begin_time_str = item.attributes['BeginTime'].value
        radio_reset_time_str = item.attributes['RadioResetTime'].value
    itemlist = read_xml(radio_filename)
    for i, item in enumerate(itemlist):
        serial_nb_orig = item.attributes['SerialNbr'].value

    # Input and output files are matched by the sonde serial number
    sounding_date_serial_dict[serial_nb_orig] = {'BeginTime': begin_time_str,
                                                 'RadioResetTime': radio_reset_time_str}
    tmpdir_obj.cleanup()

# Correct time information in simulated sounding and write to outputpath
for mwx_file in tqdm.tqdm(sorted(glob.glob(sim_mwx_fmt))):
    tmpdir, tmpdir_obj = getTmpDir()
    decompressed_files = np.array(decompress(mwx_file, tmpdir+'/'))
    
    ## Find files that need changes
    snd_mask = [f_snd(file) for file in decompressed_files]
    snd_filename = decompressed_files[snd_mask][0]
    radio_mask = [f_radio(file) for file in decompressed_files]
    radio_filename = decompressed_files[radio_mask][0]

    # Read Soundings.xml to get launch time
    itemlist = read_xml(radio_filename)
    for i, item in enumerate(itemlist):
        serial_nb_orig = item.attributes['SerialNbr'].value

    ## Get acctual launch and radio reset time
    begin_time_str = sounding_date_serial_dict[serial_nb_orig]['BeginTime']
    radio_reset_time_str = sounding_date_serial_dict[serial_nb_orig]['RadioResetTime']
    begin_time_dt = dt.datetime.strptime(begin_time_str,'%Y-%m-%dT%H:%M:%S.%f')
    print(begin_time_str, radio_reset_time_str)
    
    itemlist, xml_handle = read_xml(snd_filename, return_handle=True)
    for i, item in enumerate(itemlist):
        item.attributes['BeginTime'].value = begin_time_str
        item.attributes['RadioResetTime'].value = radio_reset_time_str
    xml_handle.writexml(open(snd_filename, 'w'))
    
    ## Update filename with actual launch time
    orig_fn = os.path.basename(mwx_file)
    orig_fn_split = orig_fn.split('_')
    orig_fn_split[1] = begin_time_dt.strftime('%Y%m%d')
    orig_fn_split[2] = begin_time_dt.strftime('%H%M%S')
    updated_fn = '_'.join(orig_fn_split)

    compress(tmpdir, output_path+updated_fn)
    tmpdir_obj.cleanup()
