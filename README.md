# eurec4a_snd

| :warning: **This is not yet the final version** Feel free to try and submit issues that you encounter. |
| --- |

This package is ment to help in establishing a common radiosonde output format.
During the EUREC4A/ATOMIC campagne many radiosondes will be launched and an unique fileformat incl. unique variables names will hopefully help there community.

Currently, this package contains a script to convert radiosonde datafile of the `MW41` to netCDF files.

Simply run `sounding_converter -i your_sounding.bfr` after the installation. During the first execution you will be asked to give meta-information which will be included in the output files.

These include e.g. `contact_person` and `platform` name which are probably helpful after the exchange of these files with other scientists.

## Installation

The package can be installed with `conda` (https://www.anaconda.com/distribution/):

The best option is to create a new environment to not run into dependency problems with pre-installed package:

```bash
conda create --name field_campaign
```

Activate the new environment with

```bash
conda activate field_campaign
```
or
```bash
source activate field_campaign
```

The actual package can than be installed with
```bash
conda install -c observingclouds -c conda-forge eurec4a_snd
```

## First execution

During the first execution, the `meta_information.ini` configuration file has to be created. The user is prompted some questions on the command line before the file will be written to the users home directory.

It is also possible to copy the `config/meta_information_template.ini` to one's favorite folder and change the files content accordingly. If the file is not renamed `meta_information.ini` and is in the home directory, the option `-c` has to be used during the call of `sounding_converter` otherwise it will not recognize the configuration.

```python
sounding_converter -i your_sounding.bfr -c /your/path/to/meta_information.ini
```

## Update

Please check for updates at the beginning of the campaign by running
```bash
conda install -c observingclouds -c conda-forge eurec4a_snd
```

## Change meta-information
The meta-information which is used in the netCDF files and in order to create a resonable filename, can always be changed by either

- remove the file `meta_information.ini` in `$HOME` and run `sounding_converter` again, as if it would be your first execution
- edit the file `~/meta_information.ini` directly

## Visualization
### Panoply
The converted dat files are netCDF files which conform to the CF-Conventions as far as possible and make use of the `discrete sampling geometry`. The sounding data can therefore easily be drawn as trajectories without extra efforts. One example is the software [Panoply](https://www.giss.nasa.gov/tools/panoply/)

![Trajectory visualization with panoply](docs/panoply_visualization_traj.png?raw=true "Trajectory visualization with panoply")

### Simple Plotting
The package also includes a few plotting routines that can be called with e.g.

```bash
sounding_visualize -n converted/file/sounding.nc
```

### SkewT Plotting
A skewT diagram can be created with
```bash
sounding_skewT -i converted/file/sounding.nc
```
Further examples on how to create a skewT diagram can be found in `eurec4a_snd/examples/visualizing`. Have a look at [Skew-T examples](eurec4a_snd/examples/visualizing/README.md)

## Example

Examples of the input `.bfr` files and the converted `.nc` file can be found in `examples/data`.

## Trouble shooting

<details>
  <summary>ImportError: DLL load failed: The specific module could not be found</summary>
  <br>
  Windows users might get the above error message when trying to visualize the soundings. The error is caused in pillow. Unfortunately there is not a very good solution yet, but the following might work for you:
  <pre>conda remove --force pillow<br>pip install pillow</pre>
  <br>
</details>
 
<details>
  <summary>Commands `sounding_convert` and `sounding_visualize` cannot be found!</summary>
  <br>
  It seems something went wrong with the installation via `conda`. Although it is recommended to install this package via anaconda because it comes with the benefit that all dependencies should be resolved, you can also download this git repository and run within the `eurec4a_snd` folder:
  <pre>python L1_rs41.py -i your_sounding.bfr</pre>
or for quicklooks
  <pre>python make_quicklooks_rs41.py -n converted/file/sounding.nc</pre>
  <br>
</details>

<details>
  <summary>Slow internet connection: Download failed</summary>
  <br>
  In case of a slow internet connection, the command `conda install some_package` might fail due to connection timeout. In this case it might be a good option to download the failing package manually.

  In this case the `eurec4a_snd` package for OSX is downloaded and installed manually:
  <pre>wget -c https://anaconda.org/observingClouds/eurec4a_snd/v0.0.41/download/osx-64/eurec4a_snd-v0.0.41-py37_0.tar.bz2</pre>
The path needs to be adapted depending on the operating system and the version that should be downloaded. You may actually see the path you need to download in the error message of `conda install eurec4a_snd`.

The installation follows simply with
  <pre>conda install -c observingClouds eurec4a_snd-v0.0.41-py37_0.tar.bz2</pre>
</details>
