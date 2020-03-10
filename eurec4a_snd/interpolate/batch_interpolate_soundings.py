import glob
import xarray as xr
import numpy as np
import metpy
from metpy import calc
from metpy.units import units
import matplotlib.pyplot as plt
from postprocessing import *
import tqdm
from netCDF4 import num2date, default_fillvals

for subfolder in ['MET','MER','BCO','ATL','RBR']:
    file_fmt = f"/mnt/lustre02/work/mh1119/level2_allsoundings/{subfolder}/*.nc"
    variables_dict = {'temperature': 'temperature', 'dewPoint': 'dew_point',
                      'flight_time': 'flight_time',  'launch_time': 'launch_time',
                      'windSpeed': 'wind_speed', 'pressure': 'pressure',
                      'wind_u': 'wind_u', 'wind_v': 'wind_v',
                      'latitude': 'latitude', 'longitude':'longitude'}
    output_variables = ['altitude','temperature', 'pressure',
                        'dew_point', 'wind_u', 'wind_v', 'wind_speed', 
                        'longitude', 'latitude', 'mixing_ratio', 'launch_time',
                        'flight_time']
    meta_data_dict = {'flight_time': {'long_name': 'time at pressure level',
                                      'units': 'seconds since 1970-01-01 00:00:00 UTC',
                                      'coordinates': 'flight_time longitude latitude altitude',
                                      '_FillValue' : default_fillvals['f8']},
                      'launch_time': {'long_name': 'time at which the sounding started',
                                      'units': 'seconds since 1970-01-01 00:00:00 UTC',
                                      '_FillValue' : default_fillvals['f8']
                                     },
                      'ascent_rate': {'long_name': 'ascent/desent rate of measuring device',
                                      'coordinates': 'flight_time longitude latitude altitude',
                                      'description': 'calculated from interpolated geootential height changes',
                                      'units': 'gpm/s',
                                      '_FillValue' : default_fillvals['f4']},
                      'altitude': {'long_name': 'geopotential height',
                                   'standard_name': 'geopotential_height',
                                   'units': 'gpm',
                                   'axis': 'Y',
                                   '_FillValue' : default_fillvals['f4']
                                  },
                      'pressure': {'long_name': 'pressure',
                                   'coordinates': 'flight_time longitude latitude altitude',
                                   '_FillValue' : default_fillvals['f4']
                                  },
                      'temperature': {'long_name': 'dry bulb temperature',
                                      'coordinates': 'flight_time longitude latitude altitude',
                                      '_FillValue' : default_fillvals['f4']
                                     },
                      'relative_humidity': {'long_name': 'relative_humidity',
                                            'standard_name': 'relative_humidity',
                                            'coordinates': 'flight_time longitude latitude altitude',
                                            'units': '%',
                                            '_FillValue' : default_fillvals['f4']
                                           },
                      'specific_humidity': {'long_name': 'specific humidity',
                                           'standard_name': 'specific_humidity',
                                           'units': 'g/kg',
                                           '_FillValue' : default_fillvals['f4']
                                          },
                      'dew_point': {'long_name': 'dew point temperature',
                                    'coordinates': 'flight_time longitude latitude altitude',
                                    '_FillValue' : default_fillvals['f4']},
                      'mixing_ratio': {'long_name': 'water vapor mixing ratio',
                                       'coordinates': 'flight_time longitude latitude altitude',
                                       'units': 'g/kg',
                                       'standard_name': 'humidity_mixing_ratio',
                                       '_FillValue' : default_fillvals['f4']
                                      },
                      'wind_speed': {'long_name': 'wind speed',
                                     'coordinates': 'flight_time longitude latitude altitude',
                                     '_FillValue' : default_fillvals['f4']
                                    },
                      'wind_direction': {'long_name': 'wind direction',
                                         'coordinates': 'flight_time longitude latitude altitude',
                                         'standard_name': 'wind_from_direction',
                                         'units': 'degrees',
                                         '_FillValue' : default_fillvals['f4']
                                        },
                      'wind_u': {'long_name': 'u-component of the wind',
                                 'standard_name': 'eastward_wind',
                                 'units': 'm/s',
                                 'coordinates':'flight_time longitude latitude altitude',
                                 '_FillValue' : default_fillvals['f4']
                                },
                      'wind_v': {'long_name': 'v-component of the wind',
                                 'standard_name': 'northward_wind',
                                 'units': 'm/s',
                                 'coordinates':'flight_time longitude latitude altitude',
                                 '_FillValue' : default_fillvals['f4']
                                },
                      'latitude': {'long_name': 'latitude',
                                   '_FillValue' : default_fillvals['f4']
                                  },
                      'longitude': {'long_name': 'longitude',
                                    '_FillValue' : default_fillvals['f4']
                                   },
                      'ascent_flag': {'long_name': 'indicator of vertical flight direction',
                                      'flag_values': '1, 0',
                                      'flag_meanings': 'ascending descending',
                                      'valid_range': '0, 1'
                                     },
                      'platform': {'long_name': 'platform identifier',
                                   'units': '-',
                                   'description': '1: BCO, 2: Meteor, 3: RH-Brown, 4: MS-Merian, 5: Atalante'
                                  }
                     }
    platform_rename_dict = {'Atalante (ATL)': 'Atalante',
                            'Barbados Cloud Observatory (BCO)': 'BCO',
                            'Maria S Merian (MER)': 'MS-Merian',
                            'FS_Meteor (M161)': 'Meteor',
                            'Research Vessel Ronald H. Brown (WTEC) (RHB)':'RH-Brown'}
    platform_number_dict = {'Atalante (ATL)': 5,
                            'Barbados Cloud Observatory (BCO)': 1,
                            'Maria S Merian (MER)': 4,
                            'FS_Meteor (M161)': 2,
                            'Research Vessel Ronald H. Brown (WTEC) (RHB)':3}
    std_pressures = [1000.00, 925.00, 850.00, 700.00, 500.00, 400.00, 300.00, 250.00, 200.00, 150.00, 100.00, 70.00, 50.00]

    files = sorted(glob.glob(file_fmt))

    for f, file in tqdm.tqdm(enumerate(files)):
        ds = xr.open_dataset(file)
        ds = ds.isel({'sounding':0})
        ds_input = ds.copy()

        # Remove standard pressure levels (extendedVerticalSoundingSignificance)
        ## extendedVerticalSoundingSignificance == 65536
        non_std_level_mask = ~np.isin(ds.pressure, std_pressures)
        ds = ds.isel({'levels':non_std_level_mask})
        ## extendedVerticalSoundingSignificance == 18432
        non_nan_altitude = ~np.isnan(ds.altitude.values)
        ds = ds.isel({'levels':non_nan_altitude})
        # Geopotential height issue
        ## the geopotential height is not a measured coordinate and
        ## the same height can occur at different pressure levels
        ## here the first occurance is used
        _, uniq_altitude_idx = np.unique(ds.altitude.values, return_index=True)
        ds = ds.isel({'levels':uniq_altitude_idx})

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
        ds['wind_u'] = -1*ds.windSpeed.values * np.sin(wind_dir_rad)
        ds['wind_v'] = -1*ds.windSpeed.values * np.cos(wind_dir_rad)

        mixing_ratio = np.array(calc_mixing_ratio_hardy(ds['dewPoint'].values+273.15,
                                                        ds['pressure'].values*100))*1000
        da_w = xr.DataArray([mixing_ratio],
                            dims = ['sounding', 'altitude'],
                            coords={'altitude':ds.altitude.values})

        ds_new = xr.Dataset()
        for variable_name_in, variable_name_out in variables_dict.items():
            try:
                ds_new[variable_name_out] = xr.DataArray([ds[variable_name_in].values],
                                                         dims = ['sounding', 'altitude'],
                                                         coords={'altitude':ds.altitude.values}
                                                        )
            except (KeyError, ValueError):
                ds_new[variable_name_out] = xr.DataArray([ds[variable_name_in].values],
                                                         dims = ['sounding']
                                                        )
            ds_new[variable_name_out].attrs = ds[variable_name_in].attrs  # Copy attributes from input
        ds_new['mixing_ratio'] = da_w

        ds_new = ds_new.dropna(dim='altitude',subset=output_variables, how='any')

        flight_time_unix = ds_new.flight_time.astype(float)/1e9
        ds_new['flight_time'].values = flight_time_unix

        # Interpolation
        ds_interp = ds_new.interp(altitude=np.arange(0,30000,10))
        wind_direction = np.rad2deg(np.arctan2(-1*ds_interp.isel({'sounding':0})['wind_u'],-1*ds_interp.isel({'sounding':0})['wind_v']))%360
        ds_interp['wind_direction'] = xr.DataArray([np.array(wind_direction.values)],
                                                   dims = ['sounding', 'altitude'],
                                                   coords={'altitude':ds_interp.altitude.values})
        ds_input = ds_input.sortby('altitude')
        ds_input.altitude.load()
        ds_input.pressure.load()
        interp_pres = pressure_interpolation(ds.pressure.values, ds.altitude.values, ds_interp.altitude.values)
        ds_interp['pressure'] = xr.DataArray([np.array(interp_pres)],
                                                   dims = ['sounding', 'altitude'],
                                                   coords={'altitude':ds_interp.altitude.values})

        # Calculations after interpolation
        specific_humidity = metpy.calc.specific_humidity_from_mixing_ratio(ds_interp['mixing_ratio']/1000)
        relative_humidity = metpy.calc.relative_humidity_from_specific_humidity(specific_humidity,
                                ds_interp.temperature.values * units.degC,
                                ds_interp.pressure.values * units['hPa'])

        ds_interp['specific_humidity'] = xr.DataArray(np.array(specific_humidity*1000),
                                                      dims = ['sounding', 'altitude'],
                                                      coords={'altitude':ds_interp.altitude.values})
        ds_interp['relative_humidity'] = xr.DataArray(np.array(relative_humidity)*100,
                                                      dims = ['sounding', 'altitude'],
                                                      coords={'altitude':ds_interp.altitude.values})
        ds_interp['launch_time'] = xr.DataArray([ds_interp.isel({'sounding':0}).launch_time.item()/1e9],
                                                dims=['sounding'])
        ds_interp['platform'] = xr.DataArray([platform_number_dict[platform]],
                                                dims=['sounding'])
        ds_interp['ascent_rate'] = xr.DataArray([calc_ascentrate(ds_interp.isel({'sounding':0}).altitude.values,
                                                                 ds_interp.isel({'sounding':0}).flight_time.values)],
                                                      dims = ['sounding', 'altitude'],
                                                      coords={'altitude':ds_interp.altitude.values})
    #     import pdb; pdb.set_trace()
    #     flight_time_start = num2date(ds_interp.isel({'sounding':0})['flight_time'].values, 'seconds since {}'.format(ds.launch_time.values))
    #     ds_interp['flight_time'] = xr.DataArray(flight_time_start,
    #                                                   dims = ['sounding', 'altitude'],
    #                                                   coords={'altitude':ds_interp.altitude.values})

        direction = get_direction(ds_interp, ds)
        if direction == 'ascending':
            ds_interp['ascent_flag'] = xr.DataArray([1], dims=['sounding'])
        else:
            ds_interp['ascent_flag'] = xr.DataArray([0], dims=['sounding'])

        global_attrs_dict = {'title': 'EUREC4A interpolated sounding data',
                             'platform_name': platform,
                             'surface_altitude': ds.attrs['surface_altitude'],
                             'instrument': ds.instrument,
                             'doi': 'pending',
                             'featureType': 'trajectory',
                             'Conventions': 'CF-1.7'
                            }

        ds_interp = set_global_attributes(ds_interp, global_attrs_dict)
        ds_interp = set_additional_var_attributes(ds_interp, meta_data_dict)
        
        for variable in ['temperature', 'dew_point', 'wind_speed', 'pressure',
                         'wind_u', 'wind_v', 'latitude', 'longitude', 'mixing_ratio', 'altitude', 'wind_direction', 'specific_humidity',
                         'relative_humidity', 'ascent_rate'
                        ]:
#             ds_interp[variable].values = np.round(ds_interp[variable].values, 2)
            ds_interp[variable].encoding = {'dtype': 'f4'}
        ds_interp['ascent_flag'].encoding = {'dtype': 'bool'}
        ds_interp['platform'].encoding = {'dtype':'uint8'}
        

        time_dt = num2date(ds_interp.isel({'sounding':0}).launch_time, "seconds since 1970-01-01 00:00:00")
        time_fmt = time_dt.strftime('%Y%m%d%H%M')
        platform_filename = platform_rename_dict[platform]
        write_dataset(ds_interp, 'EUREC4A_{platform}_soundings_{date}.nc'.format(platform=platform_filename, date=time_fmt))

