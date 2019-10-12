# eurec4a_snd

This package is ment to help in establishing a common radiosonde output format.
During the EUREC4A/ATOMIC campagne many radiosondes will be launched and an unique fileformat incl. unique variables names will hopefully help there community.

Currently, this package contains a script to convert radiosonde datafile of the `MW41` to netCDF files.

Simply run `sounding_converter -i your_sounding.dat` after the installation. During the first execution you will be asked to give meta-information which will be included in the output files.

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

The actual package can than be installed with
```bash
conda install -c observingclouds eurec4a_snd
```

## First execution

During the first execution, the `meta_information.ini` configuration file has to be created. The user is prompted some questions on the command line before the file will be written to the users home directory.

It is also possible to copy the `config/meta_information_template.ini` to ones favorite folder and change the files content accordingly. If the file in not renamed `meta_information.ini` and is in the home directory, the option `-c` has to be used during calling `sounding_converter` otherwise it will not recognize the configuration.

```python
sounding_convert -i your_sounding.dat -c /your/path/to/meta_information.ini
```

## Visualization
### Panoply
The converted dat files are netCDF files which conform to the CF-Conventions as far as possible and make use of the `discrete sampling geometry`. The sounding data can therefore easily be drawn as trajectories without extra efforts.

![Trajectory visualization with panoply](docs/panoply_visualization_traj.png?raw=true "Trajectory visualization with panoply")

### Simple Plotting
The package also includes a few plotting routines that can be called with e.g.

```bash
sounding_visualize -n converted/file/sounding.nc
```

## Example

Examples of the input `.dat` file and the converted `.nc` file can be found in `examples/data`.

## Tipps and tricks

- Commands `sounding_convert` and `sounding_visualize` cannot be found!
  It seems something went wrong with the installation via `conda`. Although it is recommended to install this package via anaconda because it comes with the benefit that all dependencies should be resolved. However, you can also download this git repository and run within the `eurec4a_snd` folder:
  
  ```bash
  python L1_rs41.py -i your_sounding.dat
  ```
  or for the quicklooks
  
  ```bash
  python make_quicklooks_rs41.py -n converted/file/sounding.nc
  ```
