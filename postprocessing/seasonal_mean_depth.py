##################################
# Take monthly average depth files
# and average seasonally
# ROMS output file names must conform to
# model_name_file_type.Y????M??D??.nc
# Minna Ho, UCLA, March 2018 
##################################
import subprocess

#########################################
# CHANGE THESE INPUTS TO CHANGE YEARS 
# AND MONTHS TO DO CALUCULATION ON
#########################################
# season
season = 'fall'

# freshw or control
run_type = 'control'

# depths
dep_en = '_0_10'

######################
# PATHS
######################
# model name
model_name = 'l2_scb'

# model file types e.g. bgc_flux_avg
#model_types = ['phys_flux','avg','bgc_flux_avg']
model_types = ['avg','bgc_flux_avg']

# path with outputs
#roms_path    = '/data/project5/kesf/ROMS/L2SCB_AP/freshw/monthly/'
roms_path    = '/data/project3/minnaho/freshwater/postprocessing/monthly_avg_depth/depth'+dep_en+'/'+run_type+'/'
# path to save monthly averages
monthly_path = '/data/project3/minnaho/freshwater/postprocessing/seasonal_avg_depth/depth'+dep_en+'/'+run_type+'/'

#########################
# get list of model file names 
# e.g. 'usw42_phys_flux.'
#########################
file_types = []
for i in model_types:
    file_types.append(model_name+'_'+i+'.')

##############################
# find seasonal average
##############################
if season == 'summer':
    print('summer')
    for f in range(len(file_types)):
        subprocess.call('ncea -O '+roms_path+file_types[f]+'Y*M0[6-8]'+dep_en+'.nc '+monthly_path+file_types[f]+'summer'+dep_en+'.nc',shell=True) 
if season == 'winter':
    print('winter')
    for f in range(len(file_types)):
        subprocess.call('ncea -O '+roms_path+file_types[f]+'Y*M12'+dep_en+'.nc '+roms_path+file_types[f]+'Y*M01'+dep_en+'.nc '+roms_path+file_types[f]+'Y*M02'+dep_en+'.nc '+monthly_path+file_types[f]+'winter'+dep_en+'.nc',shell=True) 
if season == 'spring':
    print('spring')
    for f in range(len(file_types)):
        subprocess.call('ncea -O '+roms_path+file_types[f]+'Y*M03'+dep_en+'.nc '+roms_path+file_types[f]+'Y*M04'+dep_en+'.nc '+roms_path+file_types[f]+'Y*M05'+dep_en+'.nc '+monthly_path+file_types[f]+'spring'+dep_en+'.nc',shell=True) 
if season == 'fall':
    print('fall')
    for f in range(len(file_types)):
        subprocess.call('ncea -O '+roms_path+file_types[f]+'Y*M09'+dep_en+'.nc '+roms_path+file_types[f]+'Y*M10'+dep_en+'.nc '+roms_path+file_types[f]+'Y*M11'+dep_en+'.nc '+monthly_path+file_types[f]+'autumn'+dep_en+'.nc',shell=True) 
