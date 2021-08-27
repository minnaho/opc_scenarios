import glob as glob
import subprocess as subprocess
# use zslice to get values for certain depths
# and then average over to get vertically averaged value

bgc_vars = 'GRAZE_SP,GRAZE_DIAT,GRAZE_TOT,SP_LOSS,DIAT_LOSS,ZOO_LOSS,TOT_PROD,no3_v_sp,nh4_v_sp,no3_v_diat,nh4_v_diat,SP_N_LIM,SP_LIGHT_LIM,DIAT_N_LIM,DIAT_LIGHT_LIM,GRAZE_DIAZ,DIAZ_LIGHT_LIM,NITRIF,Denitrif'

avg_vars = 'temp,salt,SPC,SPCHL,DIATC,DIATCHL,ZOOC,DIAZCHL,O2,DIC,Alk,NO3,PO4,NO2,NH4'

# freshw or control
run_type = 'control'

start_year = 2000
end_year = 2000

start_month = 4
end_month = 6

model_name = 'l2_scb'

model_types = ['avg']
#model_types = ['avg','bgc_flux_avg']

# depths
bot = -10
top = 0

file_types = []
for i in model_types:
    file_types.append(model_name+'_'+i+'.')

# grid
grid_path = '/data/project5/kesf/ROMS/L2SCB_AP/V3/roms_grd.nc'

#input_path = '/data/project3/minnaho/freshwater/postprocessing/roms_files_extract/'+run_type+'/'
input_path = './roms_files_extract/'+run_type+'/'
clim_path_day = '/data/project3/minnaho/freshwater/postprocessing/depth_avg/depth_'+str(top*-1)+'_'+str(bot*-1)+'/'+run_type+'/'
clim_path_mon = '/data/project3/minnaho/freshwater/postprocessing/monthly_avg_depth/depth_'+str(top*-1)+'_'+str(bot*-1)+'/'+run_type+'/'


# time things
months_w_31_days = [1,3,5,7,8,10,12]
leap_years = [1992,1996,2000,2004,2008,2012,2016,2020]

a = list(range(bot,top+1))
b = ' '.join(str(x) for x in a)
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
        for f in range(len(file_types)):
            if model_types[f] == 'avg':
                print('month: '+str(m))
                # days to loop over
                if m in months_w_31_days:
                    ndays = 31
                if m not in months_w_31_days:
                    ndays = 30
                    if m == 2 and y in leap_years:
                        ndays = 29
                    if m == 2 and y not in leap_years:
                        ndays = 28
                for d in list(range(1,ndays+1)):
                    date_str = 'Y'+str(y)+'M'+'%02d'%m+'D'+'%02d'%d
                    subprocess.call('zslice '+b+' '+grid_path+' '+input_path+model_name+'_'+model_types[f]+'_'+date_str+'.nc',shell=True)
                    subprocess.call('ncwa -a depth,'+str(top)+','+str(bot)+' -v '+avg_vars+' '+'z_'+model_name+'_'+model_types[f]+'_'+date_str+'.nc '+clim_path_day+file_types[f]+date_str+'_'+str(top*-1)+'_'+str(bot*-1)+'.nc',shell=True)
                    subprocess.call('rm '+'z_'+model_name+'_'+model_types[0]+'_'+date_str+'.nc',shell=True)         

            if model_types[f] == 'bgc_flux_avg':
                date_str = 'Y'+str(y)+'M'+'%02d'%m
                subprocess.call('zslice '+b+' '+grid_path+' '+input_path+model_name+'_'+model_types[f]+'_'+date_str+'.nc',shell=True)
                subprocess.call('ncwa -a depth,'+str(top)+','+str(bot)+' -v '+bgc_vars+' '+'z_'+model_name+'_'+model_types[f]+'_'+date_str+'.nc '+clim_path_mon+file_types[f]+date_str+'_'+str(top*-1)+'_'+str(bot*-1)+'.nc',shell=True)
                subprocess.call('rm '+'z_'+model_name+'_'+model_types[f]+'_'+date_str+'.nc',shell=True)

