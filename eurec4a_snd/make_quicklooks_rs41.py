"""
Script to produce quicklooks of radiosonde type RS41, read from netcdf-files
converted using L1-rs41.py

Original version by: Sabrina Schnitt
"""

import glob
import sys
import getopt
import os
import re
import math
import logging
import netCDF4 as nc
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap


def read_ncfile(ncfile):
    '''
    routine reads variables and attributes of ncfile generated by L1-rs41.py.
    INPUT: path+filename of netcdf-file
    OUTPUT: dictionnary with all variables and attributes.
    '''

    data = {}

    ncf = nc.Dataset(ncfile)

    for k in ncf.variables.keys():
        try:
            data[k] = ncf.variables[k][0, :]
        except ValueError:
            data[k] = ncf.variables[k][:]

    for a in ncf.ncattrs():
        data[a] = ncf.getncattr(a)

    specs = {}

    specs['location'] = ncf.location
    specs['tempres'] = ncf.resolution
    specs['date'] = ncf.date_YYYYMMDD
    specs['time'] = ncf.time_of_launch_HHmmss
    specs['type'] = ncf.instrument
    # Extract platform short name from the global attribute
    specs['platform_short'] = re.search(r"\(([A-Za-z0-9_]+)\)",
                                        ncf.platform_name).group(1)
    # Extract main flight direction from variable ascentRate
    most_common_vertical_movement = np.argmax(
        [np.count_nonzero(data['ascentRate'] > 0),
         np.count_nonzero(data['ascentRate'] < 0)])
    if most_common_vertical_movement == 0:
        # Mostly ascent rates
        specs['direction'] = 'AscentProfile'
    elif most_common_vertical_movement == 1:
        # Mostly descending rates
        specs['direction'] = 'DescentProfile'
    else:
        logging.warning('Main flight direction (ascent/descent) of instrument'
                        ' could not be identified!')
        specs['direction'] = 'Unknow'

    return data, specs


def plot_ptrh(data, specs, outputpath):
    '''
    routine plots vertical profiles of temperature, pressure, rel humidity and
    saves plot as .png
    INPUT:
        - data: dictionnary with data (eg filled by read_ncfile())
        - specs: dictionnary with filename specifications
            (filled by read_ncfile())
        - outputpath: path where png will be stored in.
    OUTPUT: .png file stored in outputpath
    '''
    logging.info('now plotting pressure, temperature, rel humidity sounding.........')

    # define outputname of .png-file:
    variable = 'ptrelh'
    outputname = '{platform}_{instrument}{direction}_{variable}_{date}_{tempres}.png' .format(
        platform=specs['platform_short'],
        instrument=specs['type'].replace(' ', '').replace('_', ''),
        direction=specs['direction'],
        variable=variable,
        date=specs['date']+'_'+specs['time'],
        tempres=specs['tempres'].replace(' ', ''))

    fig, ax = plt.subplots(1, 3, sharey=True, figsize=(8, 6))

    # plot temperature, pressure, humidity in three panels:
    ax[0].plot(data['temperature'], data['altitude'], '.-k', markersize=1)
    ax[1].plot(data['pressure'], data['altitude'], '.-k', markersize=1)
    ax[2].plot(data['humidity'], data['altitude'], '.-k', markersize=1)

    # do some cosmetics regarding the layout, axislabels, etc.:
    for i in range(3):
        # switch off some spines:
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)
        ax[i].spines['left'].set_visible(False)
        ax[i].grid(axis='y', linestyle='-', color='gray')

        # set height axis to start at 0m:
        ax[i].set_ylim(0, ax[i].get_ylim()[-1])
        # major minor ticks:
        ax[i].xaxis.set_minor_locator(AutoMinorLocator())
        ax[i].yaxis.set_minor_locator(AutoMinorLocator())
        ax[i].xaxis.set_major_locator(plt.MaxNLocator(4))
        # switch off major ticks top and right axis:
        ax[i].tick_params(top=False, right=False)
        # switch off minor ticks for top axis:
        ax[i].tick_params(axis='x', which='minor', top=False)
        # make labels larger for all ticks:
        ax[i].tick_params(axis='both', labelsize=14)

    ax[0].spines['left'].set_visible(True)
    ax[0].tick_params(axis='y', right=False, which='minor')
    ax[1].tick_params(left=False)
    ax[1].tick_params(axis='y', right=False, left=False, which='minor')
    ax[2].tick_params(left=False)
    ax[2].spines['right'].set_visible(True)
    ax[2].yaxis.set_ticks_position('right')
    ax[2].yaxis.set_label_position('right')

    # set the relh panel always to values between 0 and 100:
    ax[2].set_xlim(0, 100)
    # and the pressure to max 1100 hPa:
    ax[1].set_xlim(ax[1].get_xlim()[0], 1100)

    # axis labels:
    ax[0].set_ylabel('Altitude [m]', fontsize=14)
    ax[2].set_ylabel('Altitude [m]', fontsize=14)

    ax[0].set_xlabel('Temperature [$^\circ$C]', fontsize=14)
    ax[1].set_xlabel('Pressure [hPa]', fontsize=14)
    ax[2].set_xlabel('Rel Humidity [%]', fontsize=14)

    plt.subplots_adjust(top=0.9, right=0.85, left=0.15)

    fig.suptitle('%s, %s %sUTC' % (specs['location'],
                                   specs['date'],
                                   data['time_of_launch_HHmmss'][:-2]),
                 fontsize=18)

    fig.savefig(outputpath+outputname)

    logging.info('{} profiles saved at {}'.format(variable,
                                                  outputpath+outputname))


def plot_wind(data, specs, outputpath):
    '''
    routine plots vertical profiles of wind speed and direction
    and saves plot as .png
    INPUT:
        - data: dictionnary with data (eg filled by read_ncfile())
        - specs: dictionnary with filename specifications (filled
            by read_ncfile())
        - outputpath: path where png will be stored in.
    OUTPUT: .png file stored in outputpath
    '''
    logging.info('now plotting wind speed and direction sounding.........')
    # define outputname of .png-file:
    variable = 'wind'
    outputname = '{platform}_{instrument}{direction}_{variable}_{date}_{tempres}.png' .format(
        platform=specs['platform_short'],
        instrument=specs['type'].replace(' ', '').replace('_', ''),
        direction=specs['direction'],
        variable=variable,
        date=specs['date']+'_'+specs['time'],
        tempres=specs['tempres'].replace(' ', ''))

    fig, ax = plt.subplots(1, 2, sharey=True, figsize=(8, 6))

    # plot the data into subpanels:
    ax[0].plot(data['windSpeed'], data['altitude'], '.-k', markersize=1)
    ax[1].plot(data['windDirection'], data['altitude'], '.-k', markersize=1)

    # general cosmetics:
    for i in range(2):
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)
        ax[i].spines['left'].set_visible(False)
        ax[i].grid(axis='y', linestyle='-', color='gray')
        ax[i].set_ylim(0, ax[i].get_ylim()[-1])
        ax[i].xaxis.set_minor_locator(AutoMinorLocator())
        ax[i].yaxis.set_minor_locator(AutoMinorLocator())
        ax[i].xaxis.set_major_locator(plt.MaxNLocator(4))
        ax[i].tick_params(top=False, right=False)
        # axis labels:
        ax[i].set_ylabel('Altitude [m]', fontsize=14)
        # switch off minor ticks for top axis:
        ax[i].tick_params(axis='x', which='minor', top=False)
        # make labels larger for all ticks:
        ax[i].tick_params(axis='both', labelsize=14)

    # switch off some ticks and labels and spines manually.
    ax[0].spines['left'].set_visible(True)
    ax[1].tick_params(left=False)
    ax[0].tick_params(right=False, which='minor', axis='y')

    ax[1].spines['right'].set_visible(True)
    ax[1].yaxis.set_ticks_position('right')
    ax[1].yaxis.set_label_position('right')

    # set wind direction axis to valid range:
    ax[1].set_xlim(0, 360)

    ax[0].set_xlabel('Wind Speed [m s$^{-1}$]', fontsize=14)
    ax[1].set_xlabel('Wind Direction [$^\circ$]', fontsize=14)

    plt.subplots_adjust(top=0.9, right=0.85, left=0.15)
    fig.suptitle('%s, %s %sUTC' % (specs['location'],
                                   specs['date'],
                                   data['time_of_launch_HHmmss'][:-2]),
                 fontsize=18)

    fig.savefig(outputpath+outputname)

    logging.info('Wind profile saved at {}'.format(outputpath+outputname))


def plot_map(data, specs, outputpath):
    '''
    routine plots balloon flight on a map.
    INPUT:
        - data: dictionnary with data (eg filled by read_ncfile())
        - specs: dictionnary with filename specifications
            (filled by read_ncfile())
        - outputpath: path where png will be stored in.
    OUTPUT: .png file stored in outputpath.
    REQUIRES: basemap-data-hires package to be installed.
    '''
    logging.info('now plotting map of sounding.........')
    # define outputname of .png-file:
    variable = 'map'
    outputname = '{platform}_{instrument}{direction}_{variable}_{date}_{tempres}.png' .format(
        platform=specs['platform_short'],
        instrument=specs['type'].replace(' ', '').replace('_', ''),
        direction=specs['direction'],
        variable=variable,
        date=specs['date']+'_'+specs['time'],
        tempres=specs['tempres'].replace(' ', ''))

    #fig = plt.figure(figsize=(8, 6))
    fig,ax = plt.subplots(1,figsize=(8,6))
    # determine the boundaries of the map from sounding lon and lat:
    maxlon = math.ceil(np.max(data['longitude'])/0.5)*0.5
    minlon = math.floor(np.min(data['longitude'])/0.5)*0.5

    maxlat = math.ceil(np.max(data['latitude'])/0.5)*0.5
    minlat = math.floor(np.min(data['latitude'])/0.5)*0.5

    # set up basemap projection
    try:
        m = Basemap(projection='cyl', resolution='h', llcrnrlat=minlat,
                    urcrnrlat=maxlat, llcrnrlon=minlon, urcrnrlon=maxlon,
                    area_thresh=1,ax=ax)
    except OSError:
        logging.warning('High resolution map data has not been installed and'
                        ' the low resolution resolution will be used. For the'
                        ' hight resolution install with e.g. conda install -c'
                        ' conda-forge basemap-data-hires')
        m = Basemap(projection='cyl', resolution='l', llcrnrlat=minlat,
                    urcrnrlat=maxlat, llcrnrlon=minlon, urcrnrlon=maxlon,
                    area_thresh=1,ax=ax)
    # plot a topography on top:
    m.etopo(alpha=0.4)

    # coastlines, countries, boundary, background-color, gridlines
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary()
    m.shadedrelief()
    m.fillcontinents(color='#00a500')
    m.drawparallels(np.arange(10, 70, 0.25), labels=[1, 1, 0, 0])
    m.drawmeridians(np.arange(-100, 0, 0.25), labels=[0, 0, 0, 1])

    # plot balloon path:
    x, y = m(data['longitude'], data['latitude'])
    sca = m.scatter(x,y,marker='.',c=data['altitude'],cmap='Reds',vmin=0.,vmax=30000.,zorder=10)
    fig.subplots_adjust(right=0.75,left = 0.1)
    cax = plt.axes([0.85, 0.27, 0.025, 0.45])
    plt.colorbar(sca,cax=cax,label='Altitude [m]')
    #m.plot(x, y, '-k')

    # plot launch position as red square:
    m.plot(float(data['longitude'][0]),
           float(data['latitude'][0]),
           'sb', markersize=5,zorder=15)

    # and the figure title:
    ax.set_title('%s, %s %sUTC' % (specs['location'], specs['date'],
                                data['time_of_launch_HHmmss'][:-2]))

    fig.savefig(outputpath+outputname)

    logging.info('Map saved at {}'.format(outputpath+outputname))


def setup_logging(verbose):
    assert verbose in ["DEBUG", "INFO", "WARNING", "ERROR"]
    logging.basicConfig(
        level=logging.getLevelName(verbose),
        format="%(levelname)s - %(name)s - %(funcName)s - %(message)s",
        handlers=[
            logging.FileHandler("{}.log".format(__file__)),
            logging.StreamHandler()
        ])


def main():
    # options for the skript: date (yymmddhh); inputfilename; outputpath;
    # inputpath
    # (needs to be specified if date is used);
    # defaults: outputpath: './'; inputpath: './'
    # either date+inputpath or inputfilename incl path need to be specified.
    # everything else is optional.

    inputfile = ''
    date = ''
    outputpath = './'  # default setting
    inputpath = './'  # default setting
    ncfile = None

    setup_logging('INFO')

    try:
        opts, args = getopt.getopt(sys.argv[1:],
                                   'd:n:o:i:h',
                                   ['date=', 'inputncfile=', 'outputpath=',
                                    'inputpath=', 'help'])
    except getopt.GetoptError:
        print('usage: python make_quicklooks_rs41.py -p <inputpath> -d <yymmddhh>'
              ' -i <inputncfile> -o <outputpath> ')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', "--help"):  # help option
            print('usage: python make_quicklooks_rs41.py -p <inputpath> -d <yymmddhh> -i <inputncfile> -o <outputpath>')
            print('specify either complete input netcdf filename (-n) or date (-d) of sounding from which file be searched in -p inputpath. ')
            print('default inputpath: current directory. specify with -p option if different. if used together with -d, make sure -p is given first in call.')
            print('default outputpath: current directory. if -o is specified, outputpath is created if not yet existant.')
            sys.exit()

        elif opt in ("-p", "--inputpath"):
            inputpath = arg
        elif opt in ("-d", "--date"):
            try:
                ncfile = glob.glob(inputpath + '*%s*.nc' % arg)[0]
            except IndexError:
                logging.error('couldnt find your specified input: check date'
                              ' or/and inputpath selection.')
                sys.exit()
        elif opt in ("-i", "--inputncfile"):
            ncfile = arg
            if not os.path.isfile(ncfile):
                logging.error('couldnt find your specified inputfile.')
                sys.exit()

        elif opt in ("-o", "--outputpath"):
            outputpath = arg
            # check if there's a backslash after the outputpath-argument:
            if outputpath[-1] != '/':
                outputpath = outputpath+'/'
            if not os.path.isdir(outputpath):
                os.mkdir(outputpath)

        if ncfile is None:
            logging.error('Input file must be defined with either'
                          ' --inputncfile or --inputpath and --date')
            sys.exit()

    logging.info('plotting sounding file %s' % ncfile)

    # read netcdf-variables into dictionnary:
    radiosonde_data, radiosonde_specs = read_ncfile(ncfile)

    # make first quicklook: p, relh, T- profiles.
    #plot_ptrh(radiosonde_data, radiosonde_specs, outputpath)

    # now also plot wind speed and direction:
    #plot_wind(radiosonde_data, radiosonde_specs, outputpath)

    # also plot the sounding onto a map: REQUIRES BASEMAP-DATA-HIRES package
    # to be installed (e.g. through
    # conda install -c conda-forge basemap-data-hires)

    plot_map(radiosonde_data, radiosonde_specs, outputpath)


if __name__ == "__main__":
    main()
