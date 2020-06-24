#!/bin/bash
#SBATCH -J Radiosonde          # Specify job name
#SBATCH -p prepost             # Use partition prepost
#SBATCH --mem-per-cpu=244000     # Specify real memory required per CPU in MegaBytes
#SBATCH -t 12:00:00             # Set a limit on the total run time
#SBATCH -A mh0010              # Charge resources on this project account
#SBATCH -o Radiosonde.o%j          # File name for standard output
#SBATCH -e Radiosonde.e%j          # File name for standard error output

# Script to create Level1 and Level2 data from raw sounding data (level0)

module unload netcdf_c
module load anaconda3
source activate /home/mpim/m300408/conda-envs/new_campaign

cd ~/GITHUB/eurec4a_snd/eurec4a_snd/
git checkout final_changes

rm -r /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/
rm -r /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/

# Correct mwx files from the Meteor
## Get surface values from the DSHIP data
### Download DSHIP data
#mkdir DSHIP
#cd DSHIP
#wget -r --cut-dirs=100 -A dat https://observations.ipsl.fr/aeris/eurec4a-data/SHIPS/RV-METEOR/DSHIP/
#DSHIP2nc -i './observations.ipsl.fr/*.dat' -o EUREC4A_METEOR_DSHIP.nc
#rm ./observations.ipsl.fr/*.dat
#cd ..
### Write corrected surface data (incl. launch position) to mwx files
#mkdir -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level0_tosimulate/
#mkdir -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level0_torecalculate/

#### DWD
#python ../postprocessing/mwx_correction/correct_mwx.py -i '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v3/level0/MET/MWX/DWD/*.mwx' -m './DSHIP/EUREC4A_METEOR_DSHIP.nc' -o '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level0_tosimulate/'
#echo 'Please simulate the files in the Vaisala MW41 software to apply surface value corrections and retrieve the measured 1s resolution'
#echo 'Please put the simulated files into the folder /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level0_simulated/MET/')
##### Reset sounding time back to the original one (simulated soundings get time of simulation)
#python ../postprocessing/mwx_correction/reset_launchtime.py '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level0_simulated/MET/*.mwx' '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v3/level0/MET/MWX/DWD/*.mwx' '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level0_corrected/MET/'
#### non-DWD
#python ../postprocessing/mwx_correction/correct_mwx.py -i '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v3/level0/MET/MWX/*.mwx' -m './DSHIP/EUREC4A_METEOR_DSHIP.nc' -o '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level0_torecalculate/MET/'
#echo 'Please reprocess the files in the Vaisala MW41 software to apply surface value corrections'
#echo 'Please put the reprocessed files into the folder /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level0_corrected/MET/'


# Convert Level0 to Level1
#git checkout postprocessing
## Original METEOR data
#python L1_mwx.py -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v3/level0/MET/MWX/DWD/ -o /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level1_mwx/ --platform_name_long FS_METEOR --platform_name_short MET --platform_location NA
#python L1_mwx.py -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v3/level0/MET/ -o /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level1_mwx/ --platform_name_long FS_METEOR --platform_name_short MET --platform_location NA
## Corrected METEOR data
## Version_5
python L1_mwx.py -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level0_corrected/MET/ -o /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level1_mwx/ --platform_id Meteor --campaign EUREC4A --instrument_id Vaisala-RS
python L1_mwx.py -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v3/level0/MER/ -o /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level1_mwx/ --platform_id MS-Merian --campaign EUREC4A --instrument_id Vaisala-RS
python L1_mwx.py -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v3/level0/ATL/Vaisala/ -o /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level1_mwx/ --platform_id Atalante --campaign EUREC4A --instrument_id Vaisala-RS
python L1_mwx.py -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v3/level0/RHB/ -o /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level1_mwx/ --platform_id RonBrown --campaign EUREC4A --instrument_id Vaisala-RS
python L1_mwx.py -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v3/level0/BCO/ -o /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level1_mwx/ --platform_id BCO --campaign EUREC4A --instrument_id Vaisala-RS
python L1_bufr.py -p "/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v3/level0/ATL/MeteoModem/BFR/" -o /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level1_bufr/ --platform_id Atalante --campaign EUREC4A --instrument_id Meteomodem-RS -s 1 5 6 7 14


# Remove short soundings (less than 30 levels)
python ../postprocessing/remove_short_soundings.py -i '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level1_mwx/EUREC4A*L1*.nc' -t 30 -d True
python ../postprocessing/remove_short_soundings.py -i '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level1_bufr/EUREC4A*L1*.nc' -t 30 -d True

# Convert Level1 to Level2 (interpolate)
mkdir -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level1_mwx/
python ./interpolate/batch_interpolate_soundings.py -m bin -i '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level1_mwx/EUREC4A_Meteor_*Vaisala-RS*.nc' -o '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/'
python ./interpolate/batch_interpolate_soundings.py -m bin -i '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level1_mwx/EUREC4A_RonBrown_*Vaisala-RS*.nc' -o '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/'
python ./interpolate/batch_interpolate_soundings.py -m bin -i '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level1_mwx/EUREC4A_MS-Merian_*Vaisala-RS*.nc' -o '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/'
python ./interpolate/batch_interpolate_soundings.py -m bin -i '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level1_mwx/EUREC4A_BCO_*Vaisala-RS*.nc' -o '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/'
python ./interpolate/batch_interpolate_soundings.py -m bin -i '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level1_mwx/EUREC4A_Atalante_*Vaisala-RS*.nc' -o '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/'
python ./interpolate/batch_interpolate_soundings.py -m bin -i '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level1_bufr/EUREC4A_Atalante_Meteomodem-RS_L1*_20200*.nc' -o '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_bufr/'

# Concatenate Level2 files by platform
module load nco

rm /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/EUREC4A_Meteor_soundings.nc
ncrcat -h /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/*Meteor_*Vaisala-RS*_20*.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/EUREC4A_Meteor_soundings.nc
rm /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/EUREC4A_BCO_soundings.nc
ncrcat -h /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/*BCO_*Vaisala-RS*_20*.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/EUREC4A_BCO_soundings.nc
rm /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/EUREC4A_Atalante_soundings_Vaisala.nc
ncrcat -h /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/*Atalante_*Vaisala-RS*_20*.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/EUREC4A_Atalante_soundings_Vaisala.nc
rm /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/EUREC4A_MS-Merian_soundings.nc
ncrcat -h /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/*Merian_*Vaisala-RS*_20*.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/EUREC4A_MS-Merian_soundings.nc
rm /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/EUREC4A_RonBrown_soundings.nc
ncrcat -h /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/*RonBrown_*Vaisala-RS*_20*.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/EUREC4A_RonBrown_soundings.nc
rm /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_bufr/EUREC4A_Atalante_soundings_Meteomodem.nc
ncrcat -h /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_bufr/EUREC4A_Atalante_Meteomodem-RS*_20*.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_bufr/EUREC4A_Atalante_soundings_Meteomodem.nc


# Checks
#python ./tests/sounding_graphical_comparison.py

# Copy level0 data to export folder
mkdir -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_0/
cp -r /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level0_* /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_0/

# Copy level1 data to export folder
mkdir -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_1/
mkdir -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_1/MET/Vaisala/
cp /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level1_mwx/EUREC4A_Meteor_Vaisala-RS_L1*.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_1/MET/Vaisala/
mkdir -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_1/MER/Vaisala/
cp /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level1_mwx/EUREC4A_MS-Merian_Vaisala-RS_L1*.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_1/MER/Vaisala/
mkdir -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_1/ATL/Vaisala/
cp /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level1_mwx/EUREC4A_Atalante_Vaisala-RS_L1*.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_1/ATL/Vaisala/
mkdir -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_1/RHB/Vaisala/
cp /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level1_mwx/EUREC4A_RonBrown_Vaisala-RS_L1*.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_1/RHB/Vaisala/
mkdir -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_1/BCO/Vaisala/
cp /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level1_mwx/EUREC4A_BCO_Vaisala-RS_L1*.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_1/BCO/Vaisala/
mkdir -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_1/ATL/MeteoModem/
cp /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level1_bufr/EUREC4A_Atalante_Meteomodem-RS_L1*.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_1/ATL/MeteoModem/

# Export level2 data to export folder
mkdir -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/
cp /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/EUREC4A_RonBrown_soundings.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_RonBrown_soundings.nc
cp /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/EUREC4A_Atalante_soundings_Vaisala.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_Atalante_soundings_Vaisala.nc
cp /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/EUREC4A_Meteor_soundings.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_Meteor_soundings.nc
cp /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/EUREC4A_BCO_soundings.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_BCO_soundings.nc
cp /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_mwx/EUREC4A_MS-Merian_soundings.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_MS-Merian_soundings.nc
cp /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v5/level2_bufr/EUREC4A_Atalante_soundings_Meteomodem.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_Atalante_soundings_Meteomodem.nc


python ../postprocessing/remove_lower_levels.py

mv /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_RonBrown_soundings.nc2 /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_RonBrown_soundings.nc
mv /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_Atalante_soundings_Vaisala.nc2 /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_Atalante_soundings_Vaisala.nc
mv /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_Meteor_soundings.nc2 /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_Meteor_soundings.nc
mv /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_BCO_soundings.nc2 /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_BCO_soundings.nc
mv /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_MS-Merian_soundings.nc2 /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_MS-Merian_soundings.nc
mv /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_Atalante_soundings_Meteomodem.nc2 /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/EUREC4A_Atalante_soundings_Meteomodem.nc

python ../postprocessing/rename_variable_names.py
python ../postprocessing/change_units.py -i '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_2/*.nc'
python ../postprocessing/change_units.py -i '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_export2/level_1/*/*/*.nc'

echo "Please include the version number in the filename of the level 2 data. This is not done automatically."

