import os
os.environ['QT_QPA_PLATFORM']='offscreen'
import sys
sys.path.append('/data/project3/minnaho/global/')
import l2grid
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pyroms

savepath = '/data/project6/minnaho/opc_scenarios/pdf_npy/'

region_name = 'coast'

dst = 1
den = 100

varnc1 = 'rho'
varstr1 = 'rho'

varnc2 = 'salt'
varstr2 = 'salt'

# choose years
start_year = 2016
end_year = 2016

# choose months between 1 and 12
start_month = 6
end_month = 6

#exp = ['cntrl','l1617','PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
#title_exp = ['CTRL','Loads 16-17','PNDN only','PNDN 50','PNDN 90','FNDN only','FNDN 50','FNDN 90']
#exp = 'cntrl_initap'
#exp = 'loads1617'
#exp = 'PNDN_only'
#exp = 'pndn50'
#exp = 'pndn90'
#exp = 'FNDN_only'
#exp = 'fndn50'
#exp = 'fndn90'
#exp = 'PNDN_only_realistic'
#exp = 'FNDN_only_realistic'
#exp = 'pndn50_realistic'
#exp = 'pndn90_realistic'
#exp = 'fndn50_realistic'
#exp = 'fndn90_realistic'

# scenario path
#roms_path = '/data/project6/ROMS/L2SCB_OPC/'+exp+'/daily/'

# cntrl 1997 1999
#exp = 'cntrl'
#roms_path = '/data/project6/ROMS/L2SCB_1997_2000/daily/'

# cntrl 2012 2017
exp = 'cntrl_2012_2017'
roms_path = '/data/project6/ROMS/L2SCB/daily/'

#exp = 'fulll_2012_2017'
#roms_path = '/data/project6/ROMS/L2SCB_AP/daily/'

filest = 'l2_scb_avg.'

# region masks
region_mask = Dataset('/data/project1/minnaho/make_masks/mask_scb.nc','r')
mask_ssd = np.array(region_mask.variables['mask_ssd'])
mask_nsd = np.array(region_mask.variables['mask_nsd'])
mask_oc = np.array(region_mask.variables['mask_oc'])
mask_sp = np.array(region_mask.variables['mask_sp'])
mask_sm = np.array(region_mask.variables['mask_sm'])
mask_v = np.array(region_mask.variables['mask_v'])
mask_sb = np.array(region_mask.variables['mask_sb'])
# 15 km of coast
mask_cst = np.array(region_mask.variables['mask_coast'])
mask_nc = l2grid.mask_nc # full grid

mask_ssd[mask_ssd==0] = np.nan
mask_nsd[mask_nsd==0] = np.nan
mask_oc[mask_oc==0] = np.nan
mask_sp[mask_sp==0] = np.nan
mask_sm[mask_sm==0] = np.nan
mask_v[mask_v==0] = np.nan
mask_sb[mask_sb==0] = np.nan
mask_nc[mask_nc==0] = np.nan

if region_name == 'ssd':
    mask_mult = mask_ssd
    regtitle = 'South San Diego'
if region_name == 'nsd':
    mask_mult = mask_nsd
    regtitle = 'North San Diego'
if region_name == 'oc':
    mask_mult = mask_oc
    regtitle = 'Orange County'
if region_name == 'sp':
    mask_mult = mask_sp
    regtitle = 'San Pedro'
if region_name == 'sm':
    mask_mult = mask_sm
    regtitle = 'Santa Monica Bay'
if region_name == 'v':
    mask_mult = mask_v
    regtitle = 'Ventura'
if region_name == 'sb':
    mask_mult = mask_sb
    regtitle = 'Santa Barbara'
if region_name == 'grid': # full L2 grid
    mask_mult = mask_nc
    regtitle = 'SCB'
if region_name == 'coast':
    mask_mult = mask_cst
    regtitle = '15 km coast'

nbins = 500

#ext_0_80_O2_Y1998M01_12_l1617.nc
#O2
#amin
#96.9814829270808
#amax
#433.03709566462487
# bin max and min
bmin1 = 1020
bmax1 = 1029

bmin2 = 28
bmax2 = 36

# number of non-nan values 
flt1 = 0
flt2 = 0

# pdf functions
n_p1 = np.ones((nbins+1))*0
n_p2 = np.ones((nbins+1))*0

bin_p1 = np.linspace(bmin1,bmax1,nbins+1)
bin_p2 = np.linspace(bmin2,bmax2,nbins+1)

months_w_31_days = [1,3,5,7,8,10,12]
leap_years = [1992,1996,2000,2004,2008,2012,2016,2020]

for y_i in range(start_year,end_year+1):
    # if we are on the first year, starts at s_m
    if y_i == start_year:
        s_m = start_month
    else:
        s_m = 1
    # if we are on the last year, end at e_m
    if y_i == end_year:
        e_m = end_month+1
    else:
        e_m = 13
    for m_i in range(s_m,e_m):
        # loop through each file type
        if m_i in months_w_31_days:
            ndays = 31
        if m_i not in months_w_31_days:
            ndays = 30
            if m_i == 2 and y_i in leap_years:
                ndays = 29
            if m_i == 2 and y_i not in leap_years:
                ndays = 28
        for d_i in list(range(1,ndays+1)):
            dtstr = 'Y'+str(y_i)+'M'+'%02d'%m_i+'D'+'%02d'%d_i
            print(exp,dtstr)
            filenc = Dataset(roms_path+filest+dtstr+'.nc','r')
            datanc1 = np.squeeze(filenc.variables[varnc1])*mask_mult
            datanc2 = np.squeeze(filenc.variables[varnc2])*mask_mult
            if varstr1 == 'rho':
                datanc1 = datanc1+1027.4 # add rho0

            # initialize file
            varsli1 = np.ones((len(range(dst,den+dst)),mask_nc.shape[0],mask_nc.shape[1]))
            varsli2 = np.ones((len(range(dst,den+dst)),mask_nc.shape[0],mask_nc.shape[1]))

            # get grid file with zeta 
            grdz = pyroms.grid.get_ROMS_grid('L2',zeta=np.squeeze(filenc.variables['zeta']))

            # zslice 1-200m
            for d_i in range(dst,den,dst):
                varsli1[d_i-dst,:,:] = np.squeeze(pyroms.tools.zslice(datanc1,-1*d_i,grdz)[0])
                varsli2[d_i-dst,:,:] = np.squeeze(pyroms.tools.zslice(datanc2,-1*d_i,grdz)[0])

            varsli1[varsli1<=0] = np.nan
            varsli1[varsli1>1E10] = np.nan
            varsli2[varsli2<=0] = np.nan
            varsli2[varsli2>1E10] = np.nan
            
            # get PDF
            # n_d is the count, bin_d is the same for all
            n_d1,bin_d1,patch_d1 = plt.hist(datanc1.flatten(),bins=nbins,range=([bmin1,bmax1]))
            n_d2,bin_d2,patch_d2 = plt.hist(datanc2.flatten(),bins=nbins,range=([bmin2,bmax2]))
            n_d1 = np.append(n_d1,0) # make same shape as bins
            n_d2 = np.append(n_d2,0) # make same shape as bins

            fltd1 = np.where(~np.isnan(datanc1.flatten()))[0].shape[0]
            fltd2 = np.where(~np.isnan(datanc2.flatten()))[0].shape[0]

            np.save(savepath+'n_p_count_'+varstr1+'_'+dtstr+'_'+str(den)+'m_'+exp+'.npy',n_d1)
            np.save(savepath+'flt_nonan_'+varstr1+'_'+dtstr+'_'+str(den)+'m_'+exp+'.npy',fltd1)
            np.save(savepath+'bin_p_'+varstr1+'_'+dtstr+'_'+str(den)+'m_'+exp+'.npy',bin_d1)
            
            np.save(savepath+'n_p_count_'+varstr2+'_'+dtstr+'_'+str(den)+'m_'+exp+'.npy',n_d2)
            np.save(savepath+'flt_nonan_'+varstr2+'_'+dtstr+'_'+str(den)+'m_'+exp+'.npy',fltd2)
            np.save(savepath+'bin_p_'+varstr2+'_'+dtstr+'_'+str(den)+'m_'+exp+'.npy',bin_d2)

            # sum counts and number of values for PDF
            flt1 += np.where(~np.isnan(datanc1.flatten()))[0].shape[0]
            n_p1 += n_d1
            flt2 += np.where(~np.isnan(datanc2.flatten()))[0].shape[0]
            n_p2 += n_d2

#np.save(savepath+'n_p_count_'+varstr1+'_'+str(end_year)+'_'+str(den)+'m_'+exp+'_'+dtstr+'.npy',n_p1)
#np.save(savepath+'flt_nonan_'+varstr1+'_'+str(end_year)+'_'+str(den)+'m_'+exp+'_'+dtstr+'.npy',flt1)
#np.save(savepath+'bin_p_'+varstr1+'_'+str(end_year)+'_'+str(den)+'m_'+exp+'_'+dtstr+'.npy',bin_p1)
#
#np.save(savepath+'n_p_count_'+varstr2+'_'+str(end_year)+'_'+str(den)+'m_'+exp+'_'+dtstr+'.npy',n_p2)
#np.save(savepath+'flt_nonan_'+varstr2+'_'+str(end_year)+'_'+str(den)+'m_'+exp+'_'+dtstr+'.npy',flt2)
#np.save(savepath+'bin_p_'+varstr2+'_'+str(end_year)+'_'+str(den)+'m_'+exp+'_'+dtstr+'.npy',bin_p2)

