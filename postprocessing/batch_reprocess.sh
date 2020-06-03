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
source activate new_campaign

cd ~/GITHUB/eurec4a_snd/eurec4a_snd/
git checkout master
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
#python ./mwx_correction/correct_mwx.py -i '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v3/level0/MET/MWX/DWD/*.mwx' -m './DSHIP/EUREC4A_METEOR_DSHIP.nc' -o '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level0_tosimulate/'
#echo 'Please simulate the files in the Vaisala MW41 software to apply surface value corrections and retrieve the measured 1s resolution'
#echo 'Please put the simulated files into the folder /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level0_simulated/MET/')
##### Reset sounding time back to the original one (simulated soundings get time of simulation)
#python ./mwx_correction/reset_launchtime.py '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level0_simulated/MET/*.mwx' '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v3/level0/MET/MWX/DWD/*.mwx' '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level0_corrected/MET/'
#### non-DWD
#python ./mwx_correction/correct_mwx.py -i '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v3/level0/MET/MWX/*.mwx' -m './DSHIP/EUREC4A_METEOR_DSHIP.nc' -o '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level0_torecalculate/MET/'
#echo 'Please reprocess the files in the Vaisala MW41 software to apply surface value corrections'
#echo 'Please put the reprocessed files into the folder /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level0_corrected/MET/'


# Convert Level0 to Level1
## Original METEOR data
#python L1_mwx.py -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v3/level0/MET/MWX/DWD/ -o /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level1_mwx/ --platform_name_long FS_METEOR --platform_name_short MET --platform_location NA
#python L1_mwx.py -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v3/level0/MET/ -o /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level1_mwx/ --platform_name_long FS_METEOR --platform_name_short MET --platform_location NA
## Corrected METEOR data
#python L1_mwx.py -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level0_corrected/MET/ -o /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level1_mwx/ --platform_name_long FS_METEOR --platform_name_short MET --platform_location NA
#python L1_mwx.py -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v3/level0/MER/ -o /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level1_mwx/ --platform_name_long "Maria S Merian" --platform_name_short MER --platform_location NA
#python L1_mwx.py -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v3/level0/ATL/Vaisala/ -o /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level1_mwx/ --platform_name_long "Atalante" --platform_name_short ATL --platform_location NA
#python L1_mwx.py -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v3/level0/RHB/ -o /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level1_mwx/ --platform_name_long "Research Vessel Ronald H. Brown (WTEC)" --platform_name_short RHB --platform_location NA
#python L1_mwx.py -p /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v3/level0/BCO/ -o /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level1_mwx/ --platform_name_long "Barbados Cloud Observatory" --platform_name_short BCO --platform_location "Deebles Point, Barbados, West Indies"

# Convert Level1 to Level2 (interpolate)
git checkout interpolation_method
python ./interpolate/batch_interpolate_soundings.py -m bin -i '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level1_mwx/EUREC4A_MET_sounding*.nc' -o '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/'
python ./interpolate/batch_interpolate_soundings.py -m bin -i '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level1_mwx/EUREC4A_RHB_sounding*.nc' -o '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/'
python ./interpolate/batch_interpolate_soundings.py -m bin -i '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level1_mwx/EUREC4A_MER_sounding*.nc' -o '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/'
python ./interpolate/batch_interpolate_soundings.py -m bin -i '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level1_mwx/EUREC4A_BCO_sounding*.nc' -o '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/'
python ./interpolate/batch_interpolate_soundings.py -m bin -i '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level1_mwx/EUREC4A_ATL_sounding*.nc' -o '/mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/'

# Concatenate Level2 files by platform
module load nco
rm /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/EUREC4A_Meteor_soundings.nc
ncrcat -h /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/*Meteor_soundings_20*.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/EUREC4A_Meteor_soundings.nc
rm /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/EUREC4A_BCO_soundings.nc
ncrcat -h /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/*BCO_soundings_20*.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/EUREC4A_BCO_soundings.nc
rm /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/EUREC4A_Atalante_soundings_Vaisala.nc
ncrcat -h /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/*Atalante_soundings_20*.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/EUREC4A_Atalante_soundings_Vaisala.nc
rm /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/EUREC4A_MS-Merian_soundings.nc
ncrcat -h /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/*Merian_soundings_20*.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/EUREC4A_MS-Merian_soundings.nc
rm /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/EUREC4A_RH-Brown_soundings.nc
ncrcat -h /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/*RH-Brown_soundings_20*.nc /mnt/lustre02/work/mh0010/m300408/EUREC4Asoundings_v4/level2_mwx/EUREC4A_RH-Brown_soundings.nc

# Checks
python ./tests/sounding_graphical_comparison.py

