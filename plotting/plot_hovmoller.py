# plot hovmoller
import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from netCDF4 import Dataset,num2date
import cmocean
from scipy.ndimage.filters import gaussian_filter1d
from scipy.ndimage.filters import uniform_filter1d


plt.ion()
#varstr = 'temp'
varstr = 'O2'
#cblabel = 'Omega Aragonite'
cblabel = 'Dissolved Oxygen mg/L'
#cblabel = 'Temperature'

varnc = 'var'

cmapplt = cmocean.cm.dense
#cmapplt = cmocean.cm.thermal
#cmapplt = cmocean.cm.ice

exp = ['PNDN_only','pndn50','pndn90']
title_exp = ['PNDN only','PNDN 50','PNDN 90']


ncpath = '/data/project6/minnaho/opc_scenarios/ext_depth_200/'

savepath = './figs/hovmollers/'

figw = 10
figh = 5
axfont = 16


for n_i in range(len(exp)):
    datanc = Dataset(ncpath+'ext_0_200_'+varstr+'%03d'%n_i+'_interp.nc','r')
    timeunit = datanc.getncattr('start_date')
    timedim = datanc.dimensions['time'].size
    depplt = np.array(datanc.variables['depth'])
    #depplt = range(datanc.dimensions['depth'].size)
    
    dateplt = num2date(np.arange(timedim),timeunit,only_use_cftime_datetimes=False,only_use_python_datetimes=True)
    
    varplt = np.array(datanc.variables[varnc])
    
    
    fig,ax = plt.subplots(1,1,figsize=[figw,figh])
    varplt = varplt.T # transpose data
    
    if 'omega' in varnc:
        v_min = 0.5 # colorbar min
        v_max = 3 # colorbar max
        cbstp = 0.5 # colorbar tick step
        cbfmt = '%.1f' # colorbar number format
        cline = [1] # contour line
    if 'pH' in varnc:
        v_min = 7.6
        v_max = 8.2
        cbstp = 0.1
        cbfmt = '%.1f'
        cline = [7.75] # contour line
    if 'temp' in varnc:
        v_min = 6
        v_max = 20
        cbstp = 2
        cbfmt = '%d'
        cline = [] # contour line
    if 'O2' in varnc:
        v_min = 50*(16/1000)
        v_max = 300*(16/1000)
        cbstp = 0.4
        cbfmt = '%.1f'
        cline = [2] # contour line

    # smooth the values in the vertical to remove jagged lines
    if 'O2' in varnc:
        xsmoothed = gaussian_filter1d(varplt, sigma=3,axis=0)*(16./1000)
    else:
        xsmoothed = gaussian_filter1d(varplt, sigma=3,axis=0)
    #xsmoothed = gaussian_filter1d(varplt, sigma=3,axis=0)*(16/1000)
    #xsmoothed_vert = gaussian_filter1d(varplt, sigma=3,axis=0)

    # do a running mean of 7 days
    #N = 7
    #xsmoothed = uniform_filter1d(xsmoothed_vert, size=N,axis=1)

    
    #p_plt = ax.pcolormesh(dateplt,depplt,varplt,cmap=cmapplt,vmin=v_min,vmax=v_max)
    p_plt = ax.pcolormesh(dateplt,depplt,xsmoothed,cmap=cmapplt,vmin=v_min,vmax=v_max)
    #p_plt = ax.pcolormesh(dateplt,depplt,xsmoothed,cmap=cmapplt)
    #c_plt = ax.contour(dateplt,depplt,varplt,cline,colors='k')
    c_plt = ax.contour(dateplt,depplt,xsmoothed,cline,colors='k')
    
    ax.invert_yaxis()
    ax.set_ylabel('Depth (m)',fontsize=axfont)
    ax.tick_params(axis='both',which='major',labelsize=axfont)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')

    ax.set_xlim([pd.to_datetime('1997-11-01'),pd.to_datetime('1999-11-30')])

    
    tick_spacingy = 50
    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacingy))
    
    p0 = ax.get_position().get_points().flatten()
    cb_ax = fig.add_axes([p0[2]+.015,p0[1],.01,p0[3]-p0[1]])
    cb = fig.colorbar(p_plt,cax=cb_ax,orientation='vertical',format=cbfmt,ticks=np.arange(v_min,v_max+cbstp,cbstp))
    #cb = fig.colorbar(p_plt,cax=cb_ax,orientation='vertical')
    cb.set_label(cblabel,fontsize=axfont)
    cb.ax.tick_params(axis='both',which='major',labelsize=axfont)

    #ax.set_title(mpanames[n_i],fontsize=axfont)
    
    fig.savefig(savepath+'mpa_'+'%03d'%n_i+'_hovmoller_'+varstr+'_smooth.png',bbox_inches='tight')

