# Visualization examples

In this folder you can find some experimental additional scripts to visualize the converted sounding data.

## pyMeteo

The following examples needs the `pyMeteo` package installed (`conda install -c cwebster2 pymeteo`) and the `metpy` package for the calculation of `u` and `v`.

The output of `python make_skewT_pymeteo.py converted_sounding.nc` might look like
![Sounding_pymeteo](skewT_pymeteo.png?raw=true "Sounding visualization with pyMeteo")

## metpy

This example needs the `metpy` package installed (`conda install metpy`).

The output of `python make_skewT_metpy.py converted_sounding.nc` might look like
![Sounding_metpy](skewT_metpy.png?raw=true "Sounding visualization with metpy")
