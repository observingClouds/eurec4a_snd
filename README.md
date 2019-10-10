# eurec4a_snd

This package is ment to help in establishing a common radiosonde output format.
During the EUREC4A/ATOMIC campagne many radiosondes will be launched and an unique fileformat incl. unique variables names will hopefully help there community.

Currently, this package contains a script to convert radiosonde datafile of the `MW41` to netCDF files.

Simply run `sounding_converter -i your_sounding.dat` after the installation. During the first execution you will be asked to give meta-information which will be included in the output files.

These include e.g. `contact_person` and `platform` name which are probably helpful after the exchange of these files with other scientists.

## Installation

The package can be installed with `pip`:

```bash
pip install -i https://test.pypi.org/simple/ eurec4a-snd==0.0.5
```

Since this is just a test version, currently the dependencies need to be installed beforehand. The dependencies are stated in `requirements.txt`

## First execution

During the first execution, the `meta_information.ini` configuration file has to be created. The user is prompted some questions on the command line before the file will be written to the users home directory.

It is also possible to copy the `config/meta_information_template.ini` to ones favorite folder and change the files content accordingly. If the file in not renamed `meta_information.ini` and is in the home directory, the option `-c` has to be used during calling `sounding_converter` otherwise it will not recognize the configuration.

```python
sounding_convert -i your_sounding.dat -c /your/path/to/meta_information.ini
```

## Level_1 files
The converted dat files are netCDF files which conform to the CF-Conventions as far as possible and make use of the `discrete sampling geometry`. The sounding data can therefore easily be drawn as trajectories without extra efforts.

![Trajectory visualization with panoply](docs/panoply_visualization_traj.png?raw=true "Trajectory visualization with panoply")
