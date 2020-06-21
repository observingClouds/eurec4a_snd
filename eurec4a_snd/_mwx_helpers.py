import os
import tempfile
from xml.dom import minidom
import glob
import numpy as np
import subprocess
import shutil
import pandas as pd
import warnings

class SondeTypeNotImplemented(Exception):
    pass

class VariableNotFoundInSounding(Warning):
    pass

class SondeTypeNotIdentifiable(Warning):
    pass

def custom_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    return str(msg) + '\n'

warnings.formatwarning = custom_formatwarning


def getTmpDir():
    """
    Creates a temporary folder at the systems default location for temporary files.

    Returns:
        Sets Class variables:
        - self.tmpdir_obj: tempfile.TemporaryDirectory
        - self.tmpdir: string containing the path to the folder
    """
    tmpdir_obj = tempfile.TemporaryDirectory()
    tmpdir = tmpdir_obj.name
    return tmpdir, tmpdir_obj


def decompress(file, tmp_folder):
    """
    Decompress file to temporary folder
    """
    decompress_command = ["unzip", file, "-d", tmp_folder]
    _ = subprocess.check_output(decompress_command)

    decompressed_files = sorted(glob.glob(f"{tmp_folder}/*"))
    return decompressed_files


def compress(folder, compressed_file):
    """
    Compress folder to compressed file
    """
    archive = shutil.make_archive(compressed_file, 'zip', folder)
    os.rename(archive, compressed_file)
    return


def open_mwx(mwx_file):
    """
    Open Vaisala MWX41 archive file (.mwx)

    Input
    -----
    mwx_file : str
        Vaisala MW41 archive file

    Returns
    -------
    decompressed_files : list
        List of temporarily decompressed .xml files
        within the archive file
    """
    tmpdir, tmpdir_obj = getTmpDir()
    decompressed_files = np.array(decompress(mwx_file, tmpdir + '/'))
    return decompressed_files

class MWX(object):
    """
        Open Vaisala MWX41 archive file (.mwx)

        Input
        -----
        mwx_file : str
            Vaisala MW41 archive file

        Returns
        -------
        decompressed_files : list
            List of temporarily decompressed .xml files
            within the archive file
        """
    def __init__(self, mwx_file):
        self.tmpdir, self.tmpdir_obj = getTmpDir()
        self.decompressed_files = np.array(decompress(mwx_file, self.tmpdir + '/'))
    def __enter__(self):
        return self
    def __exit__(self, type, value, traceback):
        self.tmpdir_obj.cleanup()
    def get_decompressed_files(self):
        return self.decompressed_files


def check_availability(decomp_files, file, return_name=False):
    """
    Check whether xml file exist in decompressed
    file list

    Returns
    -------
    avail : bool
        Availability of file
    filename : str (optional)
        Full filename of requested file
    """
    basenames = [os.path.basename(decomp_file) for decomp_file in decomp_files]

    # Availability
    availability_mask = np.in1d(basenames, file)

    if np.sum(availability_mask) > 0:
        avail = True
    else:
        avail = False

    if return_name:
        if avail:
            idx = np.where(availability_mask)[0][0]
            fullname = decomp_files[idx]
        else:
            fullname = None
        return avail, fullname
    else:
        return avail


def read_xml(filename, return_handle=False):
    xmldoc = minidom.parse(filename)
    itemlist = xmldoc.getElementsByTagName('Row')
    if return_handle == True:
        return itemlist, xmldoc
    else:
        return itemlist


def get_sounding_profile(file, keys):
    """
    Get sounding profile from provided xml file

    Input
    -----
    file : str
        XML file containing sounding data e.g.
        SynchronizedSoundingData.xml
    keys : list
        list of variables to look for

    Returns
    -------
    pd_snd : pandas.DataFrame
        sounding profile
    """
    itemlist = read_xml(file)
    sounding_dict = {}
    try:
        for i, item in enumerate(itemlist):
            level_dict = {}
            for var in keys:
                level_dict[var] = item.attributes[var].value
            sounding_dict[i] = level_dict
    except:
        warnings.warn('Key {} not found.'.format(var), VariableNotFoundInSounding)
    pd_snd = pd.DataFrame.from_dict(sounding_dict, orient='index', dtype=float)

    # Set missing values to NaN
    pd_snd = pd_snd.replace(-32768, np.nan)
    return pd_snd


def get_sounding_metadata(file, keys):
    itemlist = read_xml(file, keys)
    sounding_meta_dict = {}
    for i, item in enumerate(itemlist):
        assert i == 0, 'further entries were found, meaning soundings meta data could be mixed up'
        for var in keys:
            try:
                sounding_meta_dict[var] = item.attributes[var].value
            except KeyError:
                warnings.warn('Attribute {} could not found and is assumed to be RS41-SGP'.format(var), SondeTypeNotIdentifiable)
                sounding_meta_dict[var] = 'RS41-SGP'
    return sounding_meta_dict

def calc_ascent_rate(sounding):
    """
    Calculate ascent rate
    """
    ascent_rate = np.diff(sounding.Height)/(np.diff(sounding.flight_time.astype(float)/1e9))
    ascent_rate = np.concatenate(([0], ascent_rate))  # 0 at first measurement
    return ascent_rate


def convert_RH_to_dewpoint(T_K, RH):
    """
    Convert T and RH to dewpoint exactly
    following the formula used by the Vaisala
    M41 sounding system
    """
    assert np.any(T_K > 100), ('Temperature seems to be not given in Kelvin')
    K = 15*np.log(100/RH) - 2*(T_K-273.15) + 2711.5
    Tdew = T_K*2*K/(T_K*np.log(100/RH)+2*K)
    
    return Tdew


def calc_vapor_pressure(sounding):
    """
    Calculate water vapor pressure
    """
    if np.any(sounding.Temperature > 100):
        print('Temperature does not seem to be given in Celsius (assume Kelvin and autoconvert to Celsius)')
        t = sounding.Temperature.values -273.15
    else:
        t = sounding.Temperature.values
    vapor_pressure = (sounding.Humidity/100.) * (611.2 * np.exp((17.62 * t)/(243.12 + t)))
    return vapor_pressure


def calc_wv_mixing_ratio(sounding, vapor_pressure):
    """
    Calculate water vapor mixing ratio
    """
    wv_mix_ratio = 1000.*((0.622*vapor_pressure)/(100.*sounding.Pressure - vapor_pressure))
    return wv_mix_ratio


def calc_temporal_resolution(time):
    """
    Calculate temporal resolution of sounding

    Returns the most common temporal resolution
    by calculating the temporal differences
    and returning the most common difference.

    Input
    -----
    time : float
        flight time
        information

    Return
    ------
    temporal_resolution : float
        temporal resolution
    """
    time_differences = np.abs(np.diff(np.ma.compressed(time)))
    time_differences_counts = np.bincount(time_differences.astype(int))
    most_common_diff = time_differences[np.argmax(time_differences_counts)]
    temporal_resolution = most_common_diff
    return temporal_resolution

f_sync = lambda file: 'SynchronizedSoundingData.xml' in file
f_snd = lambda file: 'Soundings.xml' in file
f_para = lambda file: 'SoundingParameters.xml' in file
f_std = lambda file: 'StdPressureLevels.xml' in file
f_obs = lambda file: 'SurfaceObservations.xml' in file
f_radio = lambda file: 'Radiosondes.xml' in file
