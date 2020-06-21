"""
Script to interpolate sounding data on common grid
"""
import os
import sys
import glob
import subprocess as sp
import time
import argparse
import logging
import tqdm
import numpy as np
import xarray as xr
import metpy
import metpy.calc as mpcalc
from metpy.units import units
import netCDF4
from netCDF4 import num2date, default_fillvals
import pyproj
import eurec4a_snd

pwd = os.path.dirname(os.path.abspath(__file__))
sys.path.append(pwd)
from postprocessing import *
sys.path.append(pwd+'/../')
from _helpers import get_global_attrs


def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--inputfile", metavar="INPUT_FILE",
                        help="Single converted sonde file (netCDF) or file format\n"
                             "including wildcards",
                        default=None,
                        required=False)

    parser.add_argument("-p", "--inputpath", metavar='/some/example/path',
                        help="Path to the folder containing sonde netCDF files.",
                        default=None,
                        required=False)

    parser.add_argument("-o", "--outputfolder", metavar="/some/example/path/",
                        help="Output folder for interpolated files (netCDF). You can\n"
                             " although not recommended also define an output file\n"
                             "(format).\n"
                             " The following formats can be used:\n"
                             "\t {platform}\t platform name\n"
                             "\t {location}\t platform location\n"
                             "\t {direction}\t sounding direction\n"
                             "\t {date}\t\t date of sounding release\n"
                             "\t\t or self defined date format with\n"
                             "\t\t %%Y%%m%%d %%H%%M and so on\n"
                             "\t\t and others to format the output folder dynamically.",
                        default='./',
                        required=False)

    parser.add_argument('-m', '--method', metavar='METHOD',
                        help="Interpolation method ('bin', 'linear' (default))",
                        default='linear',
                        required=False)
    
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


def setup_logging(verbose):
    assert verbose in ["DEBUG", "INFO", "WARNING", "ERROR"]
    logging.basicConfig(
        level=logging.getLevelName(verbose),
        format="%(levelname)s - %(name)s - %(funcName)s - %(message)s",
        handlers=[
            logging.FileHandler("{}.log".format(__file__)),
            logging.StreamHandler()
        ])


variables_dict = {'launch_time': 'launch_time', 'flight_time': 'flight_time',
                  'pressure':'pressure', 'latitude': 'latitude', 'longitude': 'longitude',
                  'ascentRate':'ascent_rate', 'temperature': 'temperature', 'dewPoint': 'dew_point',
                  'windSpeed': 'wind_speed',
                  'wind_u': 'wind_u', 'wind_v': 'wind_v',
                  'x': 'x', 'y':'y', 'z':'z'
                  }
output_variables = ['altitude', 'temperature', 'pressure',
                    'dew_point', 'wind_u', 'wind_v', 'wind_speed',
                    'longitude', 'latitude', 'mixing_ratio', 'launch_time',
                    'flight_time', 'ascent_rate']
meta_data_dict = {'flight_time': {'long_name': 'time at pressure level',
                                  'standard_name': 'time',
                                  'units': 'seconds since 1970-01-01 00:00:00 UTC',
                                  'coordinates': 'flight_time longitude latitude altitude',
                                  '_FillValue': default_fillvals['f8'],
                                  'ancillary_variables': 'N_ptu',
                                  'cell_method': 'alt: mean (interval=10m comment: N_ptu)'
                                  },
                  'launch_time': {'long_name': 'time at which the sounding started',
                                  'units': 'seconds since 1970-01-01 00:00:00 UTC',
                                  '_FillValue': default_fillvals['f8']
                                  },
                  'ascent_rate': {'long_name': 'ascent/desent rate of measuring device',
                                  'coordinates': 'longitude latitude flight_time sounding_id',
                                  'description': 'calculated from interpolated geopotential height changes',
                                  'units': 'm/s',
                                  '_FillValue': default_fillvals['f4'],
                                  'ancillary_variables': 'N_ptu',
                                  'cell_method': 'alt: mean (interval=10m comment: N_ptu)'
                                  },
                  'altitude': {'long_name': 'geopotential height',
                               'standard_name': 'geopotential_height',
                               'units': 'm',
                               'axis': 'Z',
                               'positive': 'up',
                               # '_FillValue': default_fillvals['f4']
                               },
                  'pressure': {'long_name': 'pressure',
                               'standard_name': 'air_pressure',
                               'units': 'hPa',
                               'coordinates': 'longitude latitude flight_time sounding_id',
                               '_FillValue': default_fillvals['f4'],
                               'ancillary_variables': 'N_ptu',
                               'cell_method': 'alt: mean (interval=10m comment: N_ptu)'
                               },
                  'temperature': {'long_name': 'dry bulb temperature',
                                  'standard_name': 'air_temperature',
                                  'units': 'degrees Celsius',
                                  'coordinates': 'longitude latitude flight_time sounding_id',
                                  '_FillValue': default_fillvals['f4'],
                                  'cell_method': 'alt: point (derived from averaged theta)'
                                  },
                  'theta': {'long_name': 'potential temperature',
                            'standard_name': 'air_potential_temperature',
                            'units': 'K',
                            'coordinates': 'flight_time longitude latitude sounding_id',
                            '_FillValue': default_fillvals['f4'],
                            'ancillary_variables': 'N_ptu',
                            'cell_method': 'alt: mean (interval=10m comment: N_ptu)'
                            }, 
                  'relative_humidity': {'long_name': 'relative_humidity',
                                        'standard_name': 'relative_humidity',
                                        'coordinates': 'flight_time longitude latitude sounding_id',
                                        'units': '%',
                                        '_FillValue': default_fillvals['f4'],
                                        'cell_method': 'alt: point (derived from averaged q, T, p)'
                                        },
                  'specific_humidity': {'long_name': 'specific humidity',
                                        'standard_name': 'specific_humidity',
                                        'units': 'g/kg',
                                        'coordinates': 'flight_time longitude latitude sounding_id',
                                        '_FillValue': default_fillvals['f4'],
                                        'ancillary_variables': 'N_ptu',
                                        'cell_method': 'alt: mean (interval=10m comment: N_ptu)'
                                        },
                  'dew_point': {'long_name': 'dew point temperature',
                                'standard_name': 'dew_point_temperature',
                                'units': 'degrees Celsius',
                                'coordinates': 'flight_time longitude latitude sounding_id',
                                '_FillValue': default_fillvals['f4'],
                                'ancillary_variables': 'N_ptu',
                                'cell_method': 'alt: mean (interval=10m comment: N_ptu)'
                                },
                  'mixing_ratio': {'long_name': 'water vapor mixing ratio',
                                   'coordinates': 'flight_time longitude latitude sounding_id',
                                   'units': 'g/kg',
                                   'standard_name': 'humidity_mixing_ratio',
                                   '_FillValue': default_fillvals['f4'],
                                   'ancillary_variables': 'N_ptu',
                                   'cell_method': 'alt: point (derived from averaged q)'
                                   },
                  'wind_speed': {'long_name': 'wind speed',
                                 'standard_name': 'wind_speed',
                                 'units': 'm/s',
                                 'coordinates': 'flight_time longitude latitude sounding_id',
                                 '_FillValue': default_fillvals['f4'],
                                 'cell_method': 'alt: point (derived from averaged u, v)'
                                 },
                  'wind_direction': {'long_name': 'wind direction',
                                     'coordinates': 'flight_time longitude latitude sounding_id',
                                     'standard_name': 'wind_from_direction',
                                     'units': 'degrees',
                                     '_FillValue': default_fillvals['f4'],
                                     'cell_method': 'alt: point (derived from averaged u, v)'
                                     },
                  'wind_u': {'long_name': 'u-component of the wind',
                             'standard_name': 'eastward_wind',
                             'units': 'm/s',
                             'coordinates': 'flight_time longitude latitude sounding_id',
                             '_FillValue': default_fillvals['f4'],
                             'ancillary_variables': 'N_gps',
                             'cell_method': 'alt: mean (interval=10m comment: N_gps)'
                             },
                  'wind_v': {'long_name': 'v-component of the wind',
                             'standard_name': 'northward_wind',
                             'units': 'm/s',
                             'coordinates': 'flight_time longitude latitude sounding_id',
                             '_FillValue': default_fillvals['f4'],
                             'ancillary_variables': 'N_gps',
                             'cell_method': 'alt: mean (interval=10m comment: N_gps)'
                             },
                  'latitude': {'long_name': 'latitude',
                               'standard_name': 'latitude',
                               'units': 'degrees_north',
                               '_FillValue': default_fillvals['f4'],
                               'ancillary_variables': 'N_gps',
                               'cell_method': 'alt: mean (interval=10m comment: N_gps)'
                               },
                  'longitude': {'long_name': 'longitude',
                                'standard_name': 'longitude',
                                'units': 'degrees_east',
                                '_FillValue': default_fillvals['f4'],
                                'ancillary_variables': 'N_gps',
                                'cell_method': 'alt: mean (interval=10m comment: N_gps)',
                                },
                  'ascent_flag': {'long_name': 'indicator of vertical flight direction',
                                  'flag_values': np.array([1, 0],dtype=np.int8),
                                  'flag_meanings': 'ascending descending',
                                  'valid_range': np.array([0, 1],dtype=np.int8)
                                  },
                  'platform': {'long_name': 'platform identifier',
                               'units': '1',
                               'description': '1: BCO, 2: Meteor, 3: RonBrown, 4: MS-Merian, 5: Atalante'
                               },
                  'N_ptu': {'standard_name': 'number_of_observations',
                         'description': 'number of measurements and method used to derive level 2 PTU data average',
                                  'flag_values':'np.nan, 0 , 1+',
                                  'flag_meanings':'no_data interpolation averaging'
                             },
                  'N_gps': {'standard_name': 'number_of_observations',
                         'description': 'number of measurements and method used to derive level 2 GPS data average',
                                  'flag_values':'np.nan, 0 , 1+',
                                  'flag_meanings':'no_data interpolation averaging'
                             },
                  'x':{'long_name':'WGS84_x'},
                  'y':{'long_name':'WGS84_y'},
                  'z':{'long_name':'WGS84_z'}
                  # 'm_ptu': {'long_name': 'cell_method used for PTU data'
                  #           'standard_name': 'status_flag',
                  #           'flag_values': '0, 1, 2',
                  #           'flag_meanings': 'no_data linear_interpolation averaging'
                  #           },
                  # 'm_gps': {'long_name': 'cell_method used for GPS data'
                  #           'standard_name': 'status_flag',
                  #           'flag_values': '0, 1, 2',
                  #           'flag_meanings': 'no_data linear_interpolation averaging'
                  #           }
                  }
platform_rename_dict = {'Atalante (ATL)': 'Atalante',
                        'Barbados Cloud Observatory (BCO)': 'BCO',
                        'Maria S Merian (MER)': 'MS-Merian',
                        'FS_METEOR (MET)': 'Meteor',
                        'Research Vessel Ronald H. Brown (WTEC) (RHB)': 'RonBrown'}
platform_number_dict = {'Atalante': 5,
                        'BCO': 1,
                        'Merian': 4,
                        'Meteor': 2,
                        'RonBrown': 3}
std_pressures = [1000.00, 925.00, 850.00, 700.00, 500.00, 400.00, 300.00,
                 250.00, 200.00, 150.00, 100.00, 70.00, 50.00]

json_config_fn = pwd+'/../mwx_config.json'


def main(args={}):
    """
    Main function
    """
    try:
        args = get_args()
    except:
        sys.exit()

    setup_logging(args['verbose'])

    logging.debug('Gathering version information')
    try:
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

    logging.debug('Create filelist')
    if args['inputfile'] is None:
        filelist = glob.glob(args['inputpath']+'*.nc')
    else:
        filelist = glob.glob(args['inputfile'])
    filelist = sorted(filelist)

    # Create outputfolder if necessary
    if not os.path.exists(args['outputfolder']):
        os.makedirs(args['outputfolder'])

    for f, file in enumerate(tqdm.tqdm(filelist)):
        logging.info(f'Process file {file}')
        ds = xr.open_dataset(file)
        ds = ds.isel({'sounding': 0})
        ds_input = ds.copy()

        # Check monotonic ascent/descent
        if np.all(np.diff(ds.isel(levels=slice(20,-1)).altitude.values) > 0) or np.all(np.diff(ds.isel(levels=slice(20,-1)).altitude.values) < 0):
            logging.debug('Sounding is monotonic ascending/descending')
        else:
            logging.warning('Sounding is not monotonic ascending/descending. The ascent rate will be artificial')

        # Remove standard pressure levels (extendedVerticalSoundingSignificance)
        # extendedVerticalSoundingSignificance == 65536
        non_std_level_mask = ~np.isin(ds.pressure, std_pressures)
        ds = ds.isel({'levels': non_std_level_mask})
        # extendedVerticalSoundingSignificance == 18432
        non_nan_altitude = ~np.isnan(ds.altitude.values)
        ds = ds.isel({'levels': non_nan_altitude})
        # Geopotential height issue
        # the geopotential height is not a measured coordinate and
        # the same height can occur at different pressure levels
        # here the first occurance is used
        _, uniq_altitude_idx = np.unique(ds.altitude.values, return_index=True)
        ds = ds.isel({'levels': uniq_altitude_idx})

        # Check if platform is known
        # if ds.platform_name not in platform_rename_dict.keys():
        #     logging.error('The platform {} is not known. Please choose one of {}'.format(ds.platform_name, platform_rename_dict.keys()))
        #     sys.exit()

        # Consistent platform test
        if f == 0:
            platform = ds.platform_name
        else:
            assert ds.platform_name == platform, 'The platform seems to change from {} to {}'.format(platform, ds.platform_name)

        # Unique levels test
        if len(ds.altitude) != len(np.unique(ds.altitude)):
            print('Altitude levels are not unique of {}'.format(file))
            break

        # Prepare some data that cannot be linearly interpolated
        wind_dir_rad = np.deg2rad(ds.windDirection.values)
        ds['wind_u'] = xr.DataArray(-1*ds.windSpeed.values * np.sin(wind_dir_rad), dims=['levels'])
        ds['wind_v'] = xr.DataArray(-1*ds.windSpeed.values * np.cos(wind_dir_rad), dims=['levels'])

        if 'altitude_WGS84' in ds.keys():
            # Convert lat, lon, alt to cartesian coordinates
            ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
            lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
            x, y, z = pyproj.transform(lla, ecef,
                                                           ds.longitude.values,
                                                           ds.latitude.values,
                                                           ds.altitude_WGS84.values,
                                                           radians=False)
            for var, val in {'x':x, 'y':y, 'z':z}.items():
                ds[var] = xr.DataArray(val, dims=['levels'])
        else:
            logging.warning('No WGS84 altitude could be found. The averaging of the position might be faulty especially at the 0 meridian and close to the poles')

        dewPoint_K = ds['dewPoint'].values+273.15
        pressure_Pa = ds['pressure'].values*100
        
        theta = calc_theta_from_T(ds['temperature'].values, ds['pressure'].values)
        #dp = mpcalc.dewpoint_from_relative_humidity(
        #ds.temperature.values * units.degC, ds.humidity.values / 100).magnitude
        q = calc_q_from_rh(ds['dewPoint'].values, ds['pressure'].values)
        e_s = calc_saturation_pressure(ds['temperature'].values+273.15)
        w_s = mpcalc.mixing_ratio(e_s*units.Pa, ds['pressure'].values*units.hPa).magnitude
        w = ds['humidity'].values/100.*w_s
        q = w/(1+w)
        #mixing_ratio = np.array(calc_mixing_ratio_hardy(dewPoint_K,
        #                                                pressure_Pa))*1000
        #q = mpcalc.specific_humidity_from_mixing_ratio(mixing_ratio/1000)
        da_w = xr.DataArray([w*1000],
                            dims=['sounding', 'altitude'],
                            coords={'altitude': ds.altitude.values})
        da_theta = xr.DataArray([theta],
                                dims=['sounding', 'altitude'],
                                coords={'altitude': ds.altitude.values})
        da_q = xr.DataArray([q*1000],
                            dims=['sounding', 'altitude'],
                            coords={'altitude': ds.altitude.values})

        ds_new = xr.Dataset()
        for variable_name_in, variable_name_out in variables_dict.items():
            try:
                ds_new[variable_name_out] = xr.DataArray([ds[variable_name_in].values],
                                                         dims=['sounding', 'altitude'],
                                                         coords={'altitude': ds.altitude.values}
                                                         )
            except (KeyError, ValueError):
                ds_new[variable_name_out] = xr.DataArray([ds[variable_name_in].values],
                                                         dims=['sounding']
                                                         )
            ds_new[variable_name_out].attrs = ds[variable_name_in].attrs  # Copy attributes from input
        ds_new['mixing_ratio'] = da_w
        ds_new['theta'] = da_theta
        ds_new['specific_humidity'] = da_q

        flight_time_unix = ds_new.flight_time.astype(float)/1e9
        ds_new['flight_time'].values = flight_time_unix

        # Interpolation
        if args['method'] == 'linear':
            ds_new = ds_new.dropna(dim='altitude',
                                   subset=output_variables,
                                   how='any')
            ds_interp = ds_new.interp(altitude=np.arange(0, 31000, 10))
        elif args['method'] == 'bin':
            ds_interp = ds_new.groupby_bins('altitude',np.arange(-5,31005,10), labels=np.arange(0,31000,10), restore_coord_dims=True).mean()
            ds_interp = ds_interp.rename({'altitude_bins':'altitude'})
            ds_interp['launch_time'] = ds_new['launch_time']

        ## Interpolation NaN
        ds_interp = ds_interp.interpolate_na('altitude', max_gap=50, use_coordinate=True)

        dims_2d = ['sounding', 'altitude']
        coords_1d = {'altitude': ds_interp.altitude.values}

        wind_u = ds_interp.isel({'sounding': 0})['wind_u']
        wind_v = ds_interp.isel({'sounding': 0})['wind_v']
        wind_direction = np.rad2deg(np.arctan2(-1*wind_u, -1*wind_v)) % 360
        wind_speed = np.sqrt(wind_u**2+wind_v**2)
        ds_interp['wind_direction'] = xr.DataArray([np.array(wind_direction.values)],
                                                   dims=dims_2d,
                                                   coords=coords_1d)
        ds_interp['wind_speed'] = xr.DataArray([np.array(wind_speed.values)],
                                                   dims=dims_2d,
                                                   coords=coords_1d)
        if 'altitude_WGS84' in ds.keys():
            lon, lat, alt = pyproj.transform(ecef, lla,
                                             ds_interp['x'].values,
                                             ds_interp['y'].values,
                                             ds_interp['z'].values,
                                             radians=False)

            for var, val in {'latitude':lat, 'longitude':lon, 'altitude_WGS84': alt}.items():
                ds_interp[var] = xr.DataArray(val, dims=dims_2d, coords=coords_1d)

            del ds_interp['x']
            del ds_interp['y']
            del ds_interp['z']
            del ds_interp['altitude_WGS84']

        ds_input = ds_input.sortby('altitude')
        ds_input.altitude.load()
        ds_input.pressure.load()
        interp_pres = pressure_interpolation(ds.pressure.values,
                                             ds.altitude.values,
                                             ds_interp.altitude.values)
        ds_interp['pressure'] = xr.DataArray([np.array(interp_pres)],
                                             dims=dims_2d,
                                             coords=coords_1d)

        # Calculations after interpolation
        #specific_humidity = metpy.calc.specific_humidity_from_mixing_ratio(ds_interp['mixing_ratio']/1000)
        #relative_humidity = metpy.calc.relative_humidity_from_specific_humidity(
        #  specific_humidity, ds_interp.temperature.values * units.degC,
        #  ds_interp.pressure.values * units.hPa)

        #ds_interp['specific_humidity'] = xr.DataArray(np.array(specific_humidity*1000),
        #                                              dims=dims_2d,
        #                                              coords=coords_1d)
        #ds_interp['relative_humidity'] = xr.DataArray(np.array(relative_humidity)*100,
        #                                              dims=dims_2d,
        #                                              coords=coords_1d)
        ds_interp['launch_time'] = xr.DataArray([ds_interp.isel({'sounding': 0}).launch_time.item()/1e9],
                                                dims=['sounding'])
        ds_interp['platform'] = xr.DataArray([platform_number_dict[platform]],
                                             dims=['sounding'])
        #ds_interp['ascent_rate'] = xr.DataArray([calc_ascentrate(ds_interp.isel({'sounding': 0}).altitude.values,
        #                                                         ds_interp.isel({'sounding': 0}).flight_time.values)],
        #                                        dims=dims_2d,
        #                                        coords=coords_1d)

        # Recalculate temperature and relative humidity from theta and q
        temperature = calc_T_from_theta(ds_interp.isel(sounding=0)['theta'].values, ds_interp.isel(sounding=0)['pressure'].values)
        ds_interp['temperature'] = xr.DataArray([np.array(temperature)],
                                                   dims=dims_2d,
                                                   coords=coords_1d)

        w = (ds_interp.isel(sounding=0)['specific_humidity'].values/1000)/(1-ds_interp.isel(sounding=0)['specific_humidity'].values/1000.)
        e_s = calc_saturation_pressure(ds_interp.isel(sounding=0)['temperature'].values+273.15)
        w_s = mpcalc.mixing_ratio(e_s*units.Pa, ds_interp.isel(sounding=0)['pressure'].values*units.hPa).magnitude
        relative_humidity = w/w_s*100
        #relative_humidity = calc_rh_from_q(ds_interp['q'].values, ds_interp['temperature_re'].values, ds_interp['pressure'].values)
        #saturation_p = calc_saturation_pressure(ds_interp['temperature_re'].values[0,:]+273.15)
        #w_s = mpcalc.mixing_ratio(saturation_p * units.Pa, ds_interp['pressure'].values * units.hPa).magnitude
        #relative_humidity = (ds_interp['specific_humidity']/1000)/((1-ds_interp['specific_humidity']/1000)*w_s)*100
 
        ds_interp['relative_humidity'] = xr.DataArray([np.array(relative_humidity)],
                                                   dims=dims_2d,
                                                   coords=coords_1d)

        # Interpolate NaNs
        ## max_gap is the maximum gap of NaNs in meters that will be still interpolated
        #ds_interp = ds_interp.interpolate_na('altitude', max_gap=50, use_coordinate=True)
        ds_interp['N_ptu'] = xr.DataArray(
            ds_new.pressure.groupby_bins('altitude', np.arange(-5, 31005, 10), labels=np.arange(0 ,31000, 10),
                                         restore_coord_dims=True).count().values,
            dims=dims_2d,
            coords=coords_1d)
        data_avail_or_interp = np.where(~np.isnan(ds_interp.pressure), 0, np.nan)
        stacked_data_counts = np.vstack([ds_interp.N_ptu.values[0, :], data_avail_or_interp[0, :]])
        nan_idx_both = np.logical_and(np.isnan(stacked_data_counts[0]), np.isnan(stacked_data_counts[1]))
        data_counts_combined = np.empty(len(stacked_data_counts[0]))
        data_counts_combined.fill(np.nan)
        data_counts_combined[~nan_idx_both] = np.nanmax(stacked_data_counts[:, ~nan_idx_both], axis=0)
        ds_interp['N_ptu'].values[0, :] = data_counts_combined

        ds_interp['N_gps'] = xr.DataArray(
            ds_new.latitude.groupby_bins('altitude', np.arange(-5, 31005, 10), labels=np.arange(0, 31000, 10),
                                         restore_coord_dims=True).count().values,
            dims=dims_2d,
            coords=coords_1d)
        data_avail_or_interp = np.where(~np.isnan(ds_interp.latitude), 0, np.nan)
        stacked_data_counts = np.vstack([ds_interp.N_gps.values[0,:], data_avail_or_interp[0,:]])
        nan_idx_both = np.logical_and(np.isnan(stacked_data_counts[0]), np.isnan(stacked_data_counts[1]))
        data_counts_combined = np.empty(len(stacked_data_counts[0]))
        data_counts_combined.fill(np.nan)
        data_counts_combined[~nan_idx_both] = np.nanmax(stacked_data_counts[:,~nan_idx_both], axis=0)
        ds_interp['N_gps'].values[0,:] = data_counts_combined

        direction = get_direction(ds_interp, ds)
        if direction == 'ascending':
            ds_interp['ascent_flag'] = xr.DataArray([1], dims=['sounding'])
        else:
            ds_interp['ascent_flag'] = xr.DataArray([0], dims=['sounding'])

        # Copy trajectory id from level1 dataset
        ds_interp['sounding_id'] = ds['sounding_id']

        script_basename = os.path.basename(__file__)
        script_modification_time = time.ctime(os.path.getmtime(os.path.realpath(__file__)))
        glob_attrs_dict = {'title': 'EUREC4A interpolated sounding data',
                             'platform_name': platform,
                             'surface_altitude': ds.attrs['surface_altitude'],
                             'instrument': ds.instrument,
                             'doi': 'pending',
                             'created_with': '{file} with its last modifications on {time}'.
                             format(time=script_modification_time,
                                    file=script_basename),
                             'git-version': git_module_version,
                             'python_version': "{} (with numpy:{}, netCDF4:{}, eurec4a_snd:{})".
                             format(sys.version, np.__version__, netCDF4.__version__, __version__),
                             'created_on': str(time.ctime(time.time())),
                             'featureType': 'trajectory',
                             'Conventions': 'CF-1.7'
                             }

        # Overwrite standard attrs with those defined in config file
        # Get global meta data from mwx_config.json
        level='L2'
        glob_attrs_dict2 = get_global_attrs(json_config_fn, f'{ds.campaign_id}_{ds.platform_id}_{ds.instrument_id}_{level}')

        for attrs, value in glob_attrs_dict2.items():
            glob_attrs_dict[attrs] = value

        ds_interp = set_global_attributes(ds_interp, glob_attrs_dict)
        ds_interp = set_additional_var_attributes(ds_interp, meta_data_dict)

        for variable in ['temperature', 'dew_point', 'wind_speed', 'pressure',
                         'wind_u', 'wind_v', 'latitude', 'longitude',
                         'mixing_ratio', 'wind_direction',
                         'specific_humidity', 'relative_humidity',
                         'ascent_rate'
                         ]:
    #             ds_interp[variable].values = np.round(ds_interp[variable].values, 2)
            ds_interp[variable].encoding['dtype'] = 'f4'
        ds_interp['ascent_flag'].encoding['dtype'] = 'int8'
        ds_interp['platform'].encoding['dtype'] = 'uint8'
        import pdb; pdb.set_trace()

        # Transpose dataset if necessary
        for variable in ds_interp.data_vars:
             dims = ds_interp[variable].dims
             if (len(dims) == 2) and (dims[0] != 'sounding'):
                 ds_interp[variable] = ds_interp[variable].T

        time_dt = num2date(ds_interp.isel({'sounding': 0}).launch_time,
                           "seconds since 1970-01-01 00:00:00")
        time_fmt = time_dt.strftime('%Y%m%d%H%M')
        platform_filename = platform  # platform_rename_dict[platform]
        outfile = args['outputfolder']+'{campaign}_{platform}_{instrument_id}_{level}_{date}_{version}.nc'.format(campaign=ds.campaign_id,
                                                                                                                  platform=platform_filename,
                                                                                                                  instrument_id=ds.instrument_id,
                                                                                                                  level=level,
                                                                                                                  date=time_fmt,
                                                                                                                  version=git_module_version
                                                                                                                  )
        logging.info('Write output to {}'.format(outfile))
        write_dataset(ds_interp, outfile)


if __name__ == '__main__':
    main()
