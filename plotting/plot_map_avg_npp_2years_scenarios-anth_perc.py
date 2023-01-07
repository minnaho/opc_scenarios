################################################
# plot maps of opc scenarios - 8 maps total
# 1 month or season
# 7 OPC scenarios
# 1 CTRL scenario
# conda activate cartopy_update 
################################################
import sys
import os
sys.path.append('/data/project3/minnaho/global/')
import l2grid
import numpy as np
from netCDF4 import Dataset,num2date
import glob as glob
import ROMS_depths as depths
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import cmocean
import datetime as datetime
import calendar
import cartopy.crs as ccrs
import cartopy.feature as cpf
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

plt.ion()

timep = 'fullts'
# contour % line
clinemin = -25
clinemax = 25

# avg maps
outpath = '/data/project6/minnaho/opc_scenarios/bgc_flux/'
filename = 'concat_fullts_int_avg_100m_50m_'

savepath = './figs/maps/'

#exp = ['PNDN_only','pndn50','pndn90']
#title_exp = ['PNDN only','PNDN 50','PNDN 90']
#exp = ['l1617','FNDN_only','fndn50','fndn90']
#title_exp = ['Loads 16-17','FNDN only','FNDN 50','FNDN 90']
#exp = ['l1617']
#title_exp = ['Loads 16-17']
#exp = ['cntrl_initap_realistic','fulll_2012_2017','PNDN_only_realistic','pndn50_fixriver','pndn90_fixriver','FNDN_only_realistic','fndn50_fixriver','fndn90_fixriver']
exp = ['cntrl_initap_realistic','fulll_2012_2017','PNDN_only_realistic','pndn50_realistic','pndn90_realistic','FNDN_only_realistic','fndn50_realistic','fndn90_realistic']
title_exp = ['CTRL','ANTH','50% N Red.','50% N Red.\n50% Recy.','50% N Red.\n90% Recy.','85% N Red.','85% N Red.\n50% Recy.','85% N Red.\n90% Recy.']

# roms var
var_nc = 'var_int'
varstr = 'npp'
cblabel = varstr+' % change'

s2d = 86400

# color map
#c_map = cmocean.cm.dense
#c_map = cmocean.cm.thermal
#c_map = cmocean.cm.algae
#c_map = cmocean.cm.deep
#c_map = cmocean.cm.delta
#c_map = 'PRGn'
#c_map = cmocean.cm.balance
#c_map1 = cmocean.cm.algae
c_map = 'PRGn'

# control scenario
cntrl_nc = outpath+filename+'cntrl_initap_realistic'+'_'+varstr+'.nc'
anth_nc = outpath+filename+'fulll_2012_2017'+'_'+varstr+'.nc'
# subtract from cntrl
# convert s^-1 to d^-1
cntrld = np.squeeze(Dataset(cntrl_nc,'r')[var_nc])*s2d
cntrld[cntrld>1E10] = np.nan
cntrld[cntrld<=0] = np.nan

anthd = np.squeeze(Dataset(anth_nc,'r')[var_nc])*s2d
anthd[anthd>1E10] = np.nan
anthd[anthd<=0] = np.nan
# time average
if timep == 'fullts':
    diffmonth = anthd-cntrld
    diffac = np.nanmean(anthd,axis=0)-np.nanmean(cntrld,axis=0)
    cntrlavg = np.nanmean(cntrld,axis=0)
# spring
if timep == 'spring':
    diffmonth = np.concatenate((anthd[5:8],anthd[17:20]))-np.concatenate((cntrld[5:8],cntrld[17:20]))
    diffac = np.nanmean((np.concatenate((anthd[5:8],anthd[17:20]))),axis=0) - np.nanmean((np.concatenate((cntrld[5:8],cntrld[17:20]))),axis=0)
    cntrlavg = np.nanmean((np.concatenate((cntrld[5:8],cntrld[17:20]))),axis=0)
# summer
if timep == 'summer':
    diffmonth = np.concatenate((anthd[8:11],anthd[20:23]))-np.concatenate((cntrld[8:11],cntrld[20:23]))
    diffac = np.nanmean((np.concatenate((anthd[8:11],anthd[20:23]))),axis=0) - np.nanmean((np.concatenate((cntrld[8:11],cntrld[20:23]))),axis=0)
    cntrlavg = np.nanmean((np.concatenate((cntrld[8:11],cntrld[20:23]))),axis=0)
# winter
if timep == 'winter':
    diffmonth = np.concatenate((anthd[2:5],anthd[14:17]))-np.concatenate((cntrld[2:5],cntrld[14:17]))
    diffac = np.nanmean((np.concatenate((anthd[2:5],anthd[14:17]))),axis=0) - np.nanmean((np.concatenate((cntrld[2:5],cntrld[14:17]))),axis=0)
    cntrlavg = np.nanmean((np.concatenate((cntrld[2:5],cntrld[14:17]))),axis=0)


# outputs
ncfiles = []
for e_i in range(len(exp)):
    ncfiles.append(outpath+filename+exp[e_i]+'_'+varstr+'.nc')

# grid path
grid_nc = l2grid.grid_nc
lat_nc = l2grid.lat_nc
lon_nc = l2grid.lon_nc
h_nc = l2grid.h_nc
pm_nc = l2grid.pm_nc
pn_nc = l2grid.pn_nc

# get size of cells
sizex = 1E-3/pm_nc
sizey = 1E-3/pn_nc

# LA/OC region
#region_name = 'laoc'
#lat_min = 33.5
#lat_max = 34.1
#lon_min = -118.9
#lon_max = -117.82

region_name = 'wider'
lat_min = 31.9
lat_max = 34.6
lon_min = -120.7
lon_max = -117

# full grid
#region_name = 'grid'
#lat_min = np.nanmin(lat_nc)
#lat_max = np.nanmax(lat_nc)
#lon_min = np.nanmin(lon_nc)
#lon_max = np.nanmax(lon_nc)

extent = [lon_min,lon_max,lat_min,lat_max]

figw = 18
figh = 8

axfont = 16

# large pipes lat and lon
major_nc = Dataset('/data/project1/minnaho/potw_outfall_data/updated_2013_2017/major_potw_data/major_potw_1971_2017_monthly.nc','r')
lat_potw = np.array(major_nc.variables['latitude'])
lon_potw = np.array(major_nc.variables['longitude'])

coast_10m = cpf.NaturalEarthFeature('physical','coastline','10m')

# max and min of color bar
v_max = 150
v_min = -150


fig,ax = plt.subplots(2,3,figsize=[figw,figh],subplot_kw=dict(projection=ccrs.PlateCarree()))

# only plot scenarios (skip cntrl and anth)
for e_i in range(2,len(exp)):
    print(title_exp[e_i])
    # convert s^-1 to d^-1
    datanc = np.squeeze(Dataset(ncfiles[e_i],'r')[var_nc])*s2d
    datanc[datanc>1E10] = np.nan
    datanc[datanc<=0] = np.nan
    # time average
    # full 2 years
    if timep == 'fullts':
        dataavg = np.nanmean(datanc,axis=0)
    # spring
    if timep == 'spring':
        dataavg = np.nanmean((np.concatenate((datanc[5:8],datanc[17:20]))),axis=0)
    # summer
    if timep == 'summer':
        dataavg = np.nanmean((np.concatenate((datanc[8:11],datanc[20:23]))),axis=0)
    # winter
    if timep == 'winter':
        dataavg = np.nanmean((np.concatenate((datanc[2:5],datanc[14:17]))),axis=0)

    subcntrl = dataavg - cntrlavg
    varplt = ((subcntrl - diffac)/diffac)*100
    #varplt = dataavg
    #p_plot = ax.flat[e_i-2].pcolormesh(lon_nc,lat_nc,varplt,transform=ccrs.PlateCarree(),cmap=c_map,vmin=v_min,vmax=v_max)
    #p_plot = ax.flat[e_i-2].pcolormesh(lon_nc,lat_nc,varplt,transform=ccrs.PlateCarree(),cmap=c_map,norm=mcolors.LogNorm(vmin=v_min,vmax=v_max))

    # interannual variability
    # difference of scenario from cntrl
    expcdiff = datanc-cntrld
    # difference of scenario from anth-cntrl
    expddiff = ((expcdiff - diffmonth)/diffmonth)*100
    

    # calculate area of change > 5 and < -5 mmol m-2 d-1
    # in km2
    npp_pos = np.where(varplt>=clinemax)
    npp_neg = np.where(varplt<=clinemin)
    pos_area = np.nansum(sizex[npp_pos[0],npp_pos[1]]*sizey[npp_pos[0],npp_pos[1]])
    neg_area = np.nansum(sizex[npp_neg[0],npp_neg[1]]*sizey[npp_neg[0],npp_neg[1]])
    print('area >'+str(clinemax)+' % change NPP: '+str(int(pos_area))+' km^2')
    print('area <'+str(clinemin)+' % change NPP: '+str(int(neg_area))+' km^2')
    

    p_plot = ax.flat[e_i-2].pcolormesh(lon_nc,lat_nc,varplt,transform=ccrs.PlateCarree(),cmap=c_map,norm=mcolors.TwoSlopeNorm(vmin=v_min,vcenter=0,vmax=v_max))
    c_plot = ax.flat[e_i-2].contour(lon_nc,lat_nc,varplt,[clinemax],transform=ccrs.PlateCarree(),colors='green')
    c_plot = ax.flat[e_i-2].contour(lon_nc,lat_nc,varplt,[clinemin],transform=ccrs.PlateCarree(),colors='purple')

    
    ax.flat[e_i-2].set_title(title_exp[e_i],fontsize=axfont)

for a_i in range(len(ax.flat)):
    # mark pipe location
    #for l_i in range(len(lon_potw)):
    #    ax.flat[a_i].scatter(lon_potw[l_i],lat_potw[l_i],marker='o',facecolors='none',edgecolors='blue',s=100)
    # other grid stuff
    ax.flat[a_i].add_feature(coast_10m,facecolor='None',edgecolor='k')
    ax.flat[a_i].add_feature(cpf.BORDERS,facecolor='None',edgecolor='k')
    ax.flat[a_i].set_extent(extent)
    gl = ax.flat[a_i].gridlines(draw_labels={'bottom':'x','left':'y'},xlabel_style=dict(size=axfont),ylabel_style=dict(size=axfont))
    if a_i >= 1 and a_i !=3 :
        gl.left_labels=False 
    if a_i >= 0 and a_i < 3 :
        gl.bottom_labels = False

# colorbar
p0 = ax.flat[2].get_position().get_points().flatten()
p1 = ax.flat[5].get_position().get_points().flatten()
cb_ax = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])

cb = fig.colorbar(p_plot,cax=cb_ax,orientation='vertical')
cb.set_label(cblabel,fontsize=axfont)
cb.ax.tick_params(axis='both',which='major',labelsize=axfont)


#savename = 'map_2years_100m_int_'+varstr+'_'+region_name+'_'+exp[-1]+'_abs_fullts.png'
savename = 'map_2years-anth_100m_int_'+varstr+'_'+region_name+'_abs_'+timep+'_perc.png'

#plt.tight_layout()

fig.savefig(savepath+savename,bbox_inches='tight')
print(savename)

