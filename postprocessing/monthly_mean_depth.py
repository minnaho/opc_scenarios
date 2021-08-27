##################################
# Take monthly averages
# of average depth files from zslice script 
# ROMS output file names must conform to
# model_name_file_type.Y????M??D??.nc
# Minna Ho, UCLA, March 2018 
##################################
import subprocess

#########################################
# CHANGE THESE INPUTS TO CHANGE YEARS 
# AND MONTHS TO DO CALUCULATION ON
#########################################
# years
st_yr = 2000
en_yr = 2000
# months
st_mon = 4
en_mon = 6

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
model_types = ['avg']

# path with outputs
#roms_path    = '/data/project5/kesf/ROMS/L2SCB_AP/freshw/monthly/'
roms_path    = '/data/project3/minnaho/freshwater/postprocessing/depth_avg/depth'+dep_en+'/'+run_type+'/'
# path to save monthly averages
monthly_path = '/data/project3/minnaho/freshwater/postprocessing/monthly_avg_depth/depth'+dep_en+'/'+run_type+'/'

#########################
# get list of model file names 
# e.g. 'usw42_phys_flux.'
#########################
file_types = []
for i in model_types:
    file_types.append(model_name+'_'+i+'.')

##############################
# find monthly average
##############################
for y in range(st_yr,en_yr+1):
    print('year: '+str(y))
    # if we are on the first year, starts at s_m
    if y == st_yr:
        s_m = st_mon
    else:
        s_m = 1
    # if we are on the last year, end at e_m
    if y == en_yr:
        e_m = en_mon+1
    else:
        e_m = 13
    for m in range(s_m,e_m):
        print('month: '+str(m))
        year_month = 'Y'+str(y)+'M'+'%02d'%m
        # loop through each file type
        subprocess.call('ncea -O '+roms_path+file_types[0]+year_month+'D*.nc '+monthly_path+file_types[0]+year_month+dep_en+'.nc',shell=True) 
