import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import l2grid
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

roms_path = '/data/project6/minnaho/opc_scenarios/ext_depth/'
savepath = '/data/project6/minnaho/opc_scenarios/pdf_npy/'

region_name = 'grid'

varnc = 'var'
varstr = 'DIN'

# choose years
start_year = 1998
end_year = 1999

# choose months between 1 and 12
start_month = 1
end_month = 11

#exp = ['cntrl','l1617','PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
#title_exp = ['CTRL','Loads 16-17','PNDN only','PNDN 50','PNDN 90','FNDN only','FNDN 50','FNDN 90']
#exp = ['cntrl']
#exp = ['l1617']
#exp = ['PNDN_only']
#exp = ['pndn50']
#exp = ['pndn90']
#exp = ['FNDN_only']
#exp = ['fndn50']
exp = ['fndn90']

filest = 'ext_0_80_'+varstr+'_'

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
bmin = 0
bmax = 400

# number of non-nan values 
flt = np.ones((len(exp)))*0

# pdf functions
n_p = np.ones((len(exp),nbins+1))*0

bin_p = np.linspace(bmin,bmax,nbins+1)

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
        if y_i == 1999 and m_i == 11:
            ndays = 27 # last day of the simulation
        for d_i in list(range(1,ndays+1)):
            dtstr = 'Y'+str(y_i)+'M'+'%02d'%m_i+'D'+'%02d'%d_i
            print(dtstr)
            for e_i in range(len(exp)):
                print(exp[e_i])
                datanc = np.squeeze(Dataset(roms_path+filest+dtstr+'_'+exp[e_i]+'.nc','r').variables[varnc])*mask_mult
                datanc[datanc<=0] = np.nan
                # n_d is the count, bin_d is the same for all
                n_d,bin_d,patch_d = plt.hist(datanc.flatten(),bins=nbins,range=([bmin,bmax]))
                n_d = np.append(n_d,0) # make same shape as bins
                # sum counts and number of values for PDF
                flt[e_i] += np.where(~np.isnan(datanc.flatten()))[0].shape[0]
                n_p[e_i] += n_d

np.save(savepath+'din_n_p_count_'+str(end_year)+'_'+exp[-1]+'.npy',n_p)
np.save(savepath+'din_flt_nonan_'+str(end_year)+'_'+exp[-1]+'.npy',flt)
np.save(savepath+'din_bin_p_'+str(end_year)+'_'+exp[-1]+'.npy',bin_p)

'''
figw = 12
figh = 8
axisfont = 16
cplt = ['green','blue','orange','orange','orange','gray','gray','gray']
lsty = ['-','-','-','--',':','-','--',':']

fig,ax = plt.subplots(1,1,figsize=[figw,figh])
for e_i in range(len(exp)):
    ax.plot(bins_p,n_p[e_i]/flt[e_i],color=cplt[e_i],linestyle=lsty[e_i],linewidth=1.5,label=title_exp[e_i])

ax.legend(fontsize=axifont)
ax.set_xlabel('O2 mmol/m3',fontsize=axisfont)
ax.set_ylabel('PDF',fontsize=axisfont)
ax.set_title(regtitle+' O2 PDF '+str(start_year)+'/'+'%02d'%start_month+'-'+str(end_year)+'/'+'%02d'%end_month,fontsize=axisfont)

ax.tick_params(axis='both',which='major',labelsize=axisfont)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')

plt.savefig(savepath+savename+'.png',bbox_inches='tight')
'''
