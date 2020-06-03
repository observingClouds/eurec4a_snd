import os
import tempfile
from xml.dom import minidom
import glob
import numpy as np
import subprocess
import shutil

class SondeTypeNotImplemented(Exception):
    pass


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


def read_xml(filename, return_handle=False):
    xmldoc = minidom.parse(filename)
    itemlist = xmldoc.getElementsByTagName('Row')
    if return_handle == True:
        return itemlist, xmldoc
    else:
        return itemlist


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
