# Post-processing of EUREC4A upper-air soundings

In general three products are available:

| Product | Description                                       |
|---------|---------------------------------------------------|
| Level 0 | Raw sounding files (Vaisala: mwx; MeteoModem: cor)|
| Level 1 | Sounding profiles on geopotential height grid<br/> - separate file for ascent and descent<br/> - netCDF|                                       
| Level 2 | Level 1 profiles gridded on a uniform geopotential<br/>height grid |

The post-processing follows [bash_reprocess.sh](https://github.com/observingClouds/eurec4a_snd/blob/postprocessing/postprocessing/batch_reprocess.sh):

![processing diagram](https://github.com/observingClouds/eurec4a_snd/blob/postprocessing/postprocessing/processing_diagram.png)
