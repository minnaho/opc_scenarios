import seawater as sw
import numpy as np
from netCDF4 import Dataset
import glob as glob
import ROMS_depths as depths
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cmocean as cmocean
import datetime as datetime
import cartopy.crs as ccrs
import cartopy.feature as cpf

plt.ion()

# depth range
dep_en = '_30_45'

# variable
var_nc = 'NO3'

# grid path
grid_path = '/data/project5/kesf/ROMS/L2SCB_AP/V3/roms_grd.nc'
grid_nc = Dataset(grid_path)
lat_nc = np.array(grid_nc.variables['lat_rho'])
lon_nc = np.array(grid_nc.variables['lon_rho'])
h_nc = np.array(grid_nc.variables['h'])

# daily path
fresh_path_day = '/data/project5/kesf/ROMS/L2SCB_AP/freshw/daily/'
contr_path_day = '/data/project5/kesf/ROMS/L2_SCB/DAILY/'

# monthly - cross section
fresh_path_mon_crs = '/data/project3/minnaho/freshwater/postprocessing/monthly_avg_all/freshw/'
contr_path_mon_crs = '/data/project3/minnaho/freshwater/postprocessing/monthly_avg_all/control/'

# seasonal - cross section
fresh_path_sea_crs = '/data/project3/minnaho/freshwater/postprocessing/seasonal_avg_all/freshw/'
contr_path_sea_crs = '/data/project3/minnaho/freshwater/postprocessing/seasonal_avg_all/control/'

# monthly path -depth avg map
fresh_path_mon_map = '/data/project3/minnaho/freshwater/postprocessing/monthly_avg_depth/depth_'+dep_en+'/freshw/'
contr_path_mon_map = '/data/project3/minnaho/freshwater/postprocessing/monthly_avg_depth/depth_'+dep_en+'/control/'

# seasonal path -depth avg map
fresh_path_sea_map = '/data/project3/minnaho/freshwater/postprocessing/seasonal_avg_depth/depth'+dep_en+'/freshw/'
contr_path_sea_map = '/data/project3/minnaho/freshwater/postprocessing/seasonal_avg_depth/depth'+dep_en+'/control/'

# read in data
fresh_map = np.array(Dataset(fresh_path_sea_map+'l2_scb_avg.summer'+dep_en+'.nc','r').variables[var_nc])
contr_map = np.array(Dataset(contr_path_sea_map+'l2_scb_avg.summer'+dep_en+'.nc','r').variables[var_nc])
fresh_map[fresh_map==0] = np.nan
contr_map[contr_map==0] = np.nan

fresh_crs = np.array(Dataset(fresh_path_sea_crs+'l2_scb_avg.summer.nc','r').variables[var_nc])
contr_crs = np.array(Dataset(contr_path_sea_crs+'l2_scb_avg.summer.nc','r').variables[var_nc])


# lat/lon used to make psource file
# choose one lat for cross section
# hyperion
htp_lat = [33.9118,33.9206,33.9017]
htp_lon = [-118.521,-118.529,-118.5267]

# jwpcp
# Y: joint_N N S L: S2 joint_S2
jwp_lat = [33.7008,33.700737,33.697917,33.6892,33.695046] 
jwp_lon = [-118.3381,-118.341962,-118.335836,-118.3167,-118.325734]

# ocsd
ocs_lat = [33.576667,33.575761] # main diffuser is 1 junction is 2
ocs_lon = [-118.01,-118.004022]

# plwtp
plw_lat = [32.665245,32.671671,32.658294]
plw_lon = [-117.323336,-117.325556,-117.324932]

# find i,j
lat_site = 33.9118 # pick one of the above
lon_site = -118.521

min_1D = np.abs( (lat_nc - lat_site)**2 + (lon_nc - lon_site)**2)
y_site, x_site = np.unravel_index(min_1D.argmin(), min_1D.shape)

lon_slice = lon_nc[y_site,:]
h_slice = h_nc[y_site,:]

lon_slice_l = list(lon_slice)*60
lon_reshape = np.array(lon_slice_l).reshape(60,602) # reshape to match z_r


# average seasonal bvfq

# calculate brunt vaisala (buoyancy) frequency
#sw.eos80.pres(depth, lat)
#sw.geostrophic.bfrq
#sw.bfrq(s, t, p, lat)[0]

# plot map and 2 cross sections
lat_min = 33.5
lat_max = 34.1
lon_min = -118.9
lon_max = -117.6

coast_10m = cpf.NaturalEarthFeature('physical','coastline','10m')

axis_font = 16

fig_w = 15
fig_h = 12

# colormap
colormaps = [cmocean.cm.solar]

var_min = 0
var_max = 10

proj_carree = ccrs.PlateCarree()
proj_stereo = ccrs.Stereographic()
proj_ortho = ccrs.Orthographic()
proj_merc = ccrs.Mercator()

proj_ch = proj_carree

extent = [lon_min,lon_max,lat_min,lat_max]

fig,axes = plt.subplots(2,2,figsize=[fig_w,fig_h],subplot_kw=dict(projection=proj_ch))

cm = plt.get_cmap(colormaps[0])
cm.set_bad('gray',alpha=0)


# map
p0 = axes.flat[0].contourf(lon_nc,lat_nc,contr_map,levels=np.linspace(var_min,var_max,100),transform=proj_ch,cmap=colormaps[0],vmin=var_min,vmax=var_max)
p1 = axes.flat[1].contourf(lon_nc,lat_nc,fresh_map,levels=np.linspace(var_min,var_max,100),transform=proj_ch,cmap=colormaps[0],vmin=var_min,vmax=var_max)

axins0 = inset_axes(axes.flat[0],width='60%',height='7%',loc='upper right')
cb = fig.colorbar(p0,cax=axins0,orientation='horizontal',ticks=np.arange(var_min,var_max+1))
cb.set_label('NO3 mmol/m3',fontsize=axis_font)
cb.ax.tick_params(axis='both',which='major',direction='in',labelsize=axis_font)

axes.flat[0].set_title('CTRL',fontsize=axis_font)
axes.flat[1].set_title('FRESH',fontsize=axis_font)
axes.flat[0].set_extent(extent,crs=proj_ch)
axes.flat[1].set_extent(extent,crs=proj_ch)

gl0 = axes.flat[0].gridlines(draw_labels=True,linestyle='--')
gl0.xlabels_top = False
gl0.ylabels_right = False
gl0.xlabel_style = {'size':axis_font}
gl0.ylabel_style = {'size':axis_font}
axes.flat[0].add_feature(coast_10m,facecolor='None',edgecolor='k')

gl1 = axes.flat[1].gridlines(draw_labels=True,linestyle='--')
gl1.xlabels_top = False
gl1.ylabels_right = False
gl1.xlabel_style = {'size':axis_font}
gl1.ylabel_style = {'size':axis_font}
axes.flat[1].add_feature(coast_10m,facecolor='None',edgecolor='k')


