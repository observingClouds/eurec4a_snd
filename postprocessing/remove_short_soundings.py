"""
Script to remove soundings that contain less than
a specific amount of levels
"""
import os
import glob
import tqdm
import argparse
import numpy as np
import xarray as xr
import logging

logging.basicConfig(filename='remove_short_soundings.log',level=logging.INFO)

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inputfilefmt', required=True,
                    help='Input filename format (level 1 files), e.g. EUREC4A_*soundings_*.nc')
parser.add_argument('-t', '--threshold', required=False, default=30,
                    help='Minimum number of levels that should exist in each sounding file, before it\n'
                    'is marked for removal.'
)
parser.add_argument('-d', '--delete', required=False,
                    help='Removal of sounding files with less than THRESHOLD amount of levels (default: False)',
                    default=False,type=str2bool)
args = vars(parser.parse_args())

level_threshold = int(args['threshold'])
files = sorted(glob.glob(args['inputfilefmt']))

logging.info('Total number of level1 soundings found: {}'.format(len(files)))

short_sounding_idx = []
for f, file in enumerate(tqdm.tqdm(files)):
    ds = xr.open_dataset(file)
    if len(ds.level) < level_threshold:
        short_sounding_idx.append(f)
logging.info('Number of short soundings found: {}'.format(len(short_sounding_idx)))
soundings_short = np.array(files)[short_sounding_idx]
logging.info('The short soundings are: {}'.format(np.array(files)[short_sounding_idx]))
if (args['delete'] is True) and (len(short_sounding_idx)> 0):
    for file in np.array(files)[short_sounding_idx]:
        os.remove(file)
    logging.info('Short soundings have been deleted')
else:
    logging.info('No soundings have been deleted')

