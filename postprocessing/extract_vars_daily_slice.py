##################################
# Find climatology from monthly averages
# of biogeochemical model outputs
# ROMS output file names must conform to
# model_name_file_type.Y????M??.nc
# Minna Ho, UCLA, March 2018 
##################################
import subprocess
from netCDF4 import Dataset

#########################################
# CHANGE THESE INPUTS TO CHANGE YEARS 
# AND MONTHS TO FIND CLIMATOLOGY 
#########################################

start_year = 1999
end_year = 2000

# between 1 and 12
start_month = 7
end_month = 9

# Sed_Flux_POC,Sed_Flux_CaCO3 omitted for zslice because 2D
#bgc_vars = 'N2O_prod'
#avg_vars = 'DIATC,SPC,DIAZC'
#varname = 'biomass'

######################
# PATHS AND FILE NAMES
# change these for different
# model names and model file types
######################
# model name
model_name   = 'l2_scb'

# path with outputs
#daily_path    = '/data/project6/ROMS/L2SCB_AP/fresh/daily/'
#exp = 'fresh'
#daily_path    = '/data/project6/ROMS/L2SCB_AP/nutrients/daily/'
#exp = 'nutri'
#daily_path    = '/data/project6/ROMS/L2SCB_P_1999_2000/daily/'
#exp = 'pipes'
daily_path    = '/data/project3/minnaho/freshwater/postprocessing/fulll_daily/'
exp = 'fulll'
#monthly_path    = '/data/project6/kesf/ROMS/L2SCB_AP/monthly_2012_2017/'

# climatology path - where extractions are saved
#clim_path    = '/data/project5/kesf/ROMS/L2SCB_AP/V3/clim/'
clim_path    = '/data/project3/minnaho/freshwater/postprocessing/surf_sli/'

grid_path = 'roms_grd.nc'
zslice_path = './symlinks/'
#input_path = '/data/project1/minnaho/organization/l2_ap/'
#################################################
# SURFACE seasonal climatology at the surface
# extract variables for each month at the surface
#################################################

'''
# bgc_avg file
# model file types e.g. bgc_flux_avg
model_types = ['bgc_flux_avg']

file_types = []
for i in model_types:
    file_types.append(model_name+'_'+i+'.')

print('calculating climatology')
for m in range(start_month,end_month+1): 
    print('month: '+str(m))
    for f in file_types:
        subprocess.call('ncra -v '+bgc_vars+' -d s_rho,59 '+daily_path+f+'Y????M'+'%02d'%m+'.nc '+clim_path+f+'M'+'%02d'%m+'_'+str(start_year)+'_'+str(end_year)+'_surf.nc',shell=True)
'''

# avg file
# model file types e.g. bgc_flux_avg
model_types = ['avg']

file_types = []
for i in model_types:
    file_types.append(model_name+'_'+i+'.')

months_w_31_days = [1,3,5,7,8,10,12]
leap_years = [1992,1996,2000,2004,2008,2012,2016,2020]

for y in range(start_year,end_year+1):
    print('year: '+str(y))
    # if we are on the first year, starts at s_m
    if y == start_year:
        s_m = start_month
    else:
        s_m = 1
    # if we are on the last year, end at e_m
    if y == end_year:
        e_m = end_month+1
    else:
        e_m = 13
    for m in range(s_m,e_m):
        print('month: '+str(m))
        year_month = 'Y'+str(y)+'M'+'%02d'%m
        # loop through each file type
        for f in file_types:
            if m in months_w_31_days:
                ndays = 31
            if m not in months_w_31_days:
                ndays = 30
                if m == 2 and y in leap_years:
                    ndays = 29
                if m == 2 and y not in leap_years:
                    ndays = 28 
            for d in list(range(1,ndays+1)):
                subprocess.call('ncra -v '+avg_vars+' -d s_rho,59 '+daily_path+f+year_month+'D'+'%02d'%d+'.nc '+clim_path+f+year_month+'D'+'%02d'%d+'_'+exp+'_'+varname+'.nc',shell=True)
                if avg_vars == 'DIATC,SPC,DIAZC':
                    subprocess.call('ncap2 -s "totc=DIATC+SPC+DIAZC" -O '+clim_path+f+year_month+'D'+'%02d'%d+'_'+exp+'_'+varname+'.nc '+clim_path+f+year_month+'D'+'%02d'%d+'_'+exp+'_totbiomass.nc',shell=True)
                
#ncrcat l2_scb_avg.Y*_totbiomass.nc sli_srf_exp_totbiomass.nc

'''

#################################################
# zslice
# extract variables for each month
#################################################

# bgc_avg file
# model file types e.g. bgc_flux_avg
model_types = ['avg','bgc_flux_avg']

file_types = []
for i in model_types:
    file_types.append(model_name+'_'+i+'.')

# make symbolic link for zslice
print('making symbolic link')
for y in range(start_year,end_year+1):
    if y != start_year:
        sm = 1
    else:
        sm = start_month
    if y == end_year:
        em = end_month
    else:
        em = 12
    for m in range(sm,em+1): 
        print('month: '+str(m))
        for f in file_types:
            date_str = 'Y'+str(y)+'M'+'%02d'%m+'.nc'
            subprocess.call('ln -fs '+monthly_path+file_types[0]+date_str+' '+zslice_path+model_name+'_'+model_types[0]+'_'+date_str,shell=True)
            subprocess.call('ln -fs '+monthly_path+file_types[1]+date_str+' '+zslice_path+model_name+'_'+model_types[1]+'_'+date_str,shell=True)

# zslice 0 - 40 m
#bot = -40
#top = 0
#a = list(range(bot,top+1))
#b = ' '.join(str(x) for x in a)
b = '-300'
for y in range(start_year,end_year+1):
    #if y != start_year:
    #    sm = 1
    #else:
    #    sm = start_month
    #if y == end_year:
    #    em = end_month
    #else:
    #    em = 12
    sm = start_month 
    em = end_month 
    for m in range(sm,em+1): 
        print('month: '+str(m))
        for f in file_types:
            date_str = 'Y'+str(y)+'M'+'%02d'%m+'.nc'
            subprocess.call('zslice '+b+' --vars='+avg_vars+' '+grid_path+' '+zslice_path+model_name+'_'+model_types[0]+'_'+date_str,shell=True)
            subprocess.call('zslice '+b+' --vars='+bgc_vars+' '+grid_path+' '+zslice_path+model_name+'_'+model_types[1]+'_'+date_str,shell=True)

ncea z_l2_scb_avg_Y????M06* ../slices/l2_scb_avg.M06_2013_2017_100m.nc
ncea z_l2_scb_avg_Y????M07* ../slices/l2_scb_avg.M07_2013_2017_100m.nc
ncea z_l2_scb_avg_Y????M08* ../slices/l2_scb_avg.M08_2013_2017_100m.nc


ncea z_l2_scb_bgc_flux_avg_Y????M06.nc ../slices/l2_scb_bgc_flux_avg.M06_2013_2017_100m.nc
ncea z_l2_scb_bgc_flux_avg_Y????M07.nc ../slices/l2_scb_bgc_flux_avg.M07_2013_2017_100m.nc
ncea z_l2_scb_bgc_flux_avg_Y????M08.nc ../slices/l2_scb_bgc_flux_avg.M08_2013_2017_100m.nc
'''
