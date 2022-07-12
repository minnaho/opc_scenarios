import numpy as np
from netCDF4 import Dataset,num2date
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
import cmocean

#plt.ion()

varstr = 'NH4'
cblabel = 'mmol/m3'
dp0 = 0
dp1 = 320
dpstp = 20
maskn = 8

fname = 'budget_L2_mask'+'%02d'%maskn+'_'+varstr+'_'

exp = ['loads1617','PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
title_exp = ['Loads 16-17','PNDN','PNDN 50','PNDN 90','FNDN','FNDN 50','FNDN 90']
#exp = ['loads1617','PNDN_only','pndn50','pndn90']
#title_exp = ['Loads 16-17','PNDN','PNDN 50','PNDN 90']
#exp = ['loads1617','FNDN_only','fndn50','fndn90']
#title_exp = ['Loads 16-17','FNDN','FNDN 50','FNDN 90']

# experiment to compare to
#compexp = 'loads1617'
compexp = 'cntrl'

ncpath = '/data/project3/minnaho/opc_scenarios/plotting/massbalance_L2/'

savepath = '/data/project3/minnaho/opc_scenarios/plotting/massbalance_L2/figs/'
savename = 'hov_budget_mask'+'%02d'%maskn+'_'+varstr

s2d = 86400

# cntrl
if compexp == 'cntrl':
    # read first file to get size of arrays
    datacomp = Dataset(ncpath+compexp+'/'+fname+str(0)+'_to_'+str(20)+'_'+compexp+'.nc','r')
    if varstr == 'O2': 
        v_min = -25
        v_max = 25
    if varstr == 'NH4': 
        #v_min = -1E-1
        #v_max = 1E-1
        if maskn == 3:
            v_min = -0.5
            v_max = 0.5
        else:
            v_min = -0.1
            v_max = 0.1
    if varstr == 'NO3': 
        #v_min = -1E-1
        #v_max = 1E-1
        v_min = -5
        v_max = 5
    
# l1617
else:
    datacomp = Dataset(ncpath+exp[0]+'/'+fname+str(0)+'_to_'+str(20)+'_'+exp[0]+'.nc','r')
    v_min = -5
    v_max = 5
    
timcomp = np.array(datacomp.variables['time'])
datetimecomp = pd.to_datetime(timcomp-719529,unit='D')
# dimensions are experiment, time and depth (3rd dimension will be variable
# for DIN/TN calculation)
inv_arr = np.ones((len(exp),len(timcomp),len(range(dp0,dp1,dpstp))))

dep_ind = 0
for d_i in range(dp0,dp1,dpstp):
    # exp to compare to
    datacomp = Dataset(ncpath+compexp+'/'+fname+str(d_i)+'_to_'+str(d_i+dpstp)+'_'+compexp+'.nc','r')
    timcomp = np.array(datacomp.variables['time'])
    timcompconv = pd.to_datetime(timcomp-719529,unit='D')
    dttcomp = np.array(datacomp.variables['dt'])
    volcomp = np.array(datacomp.variables['volume'])
    arecomp = np.array(datacomp.variables['area'])
    bgccomp = (np.array(datacomp.variables['bgc'])*s2d)/(dttcomp*arecomp)
    invcomp = np.array(datacomp.variables['invm'])/volcomp
    for e_i in range(len(exp)):
        datanc = Dataset(ncpath+exp[e_i]+'/'+fname+str(d_i)+'_to_'+str(d_i+dpstp)+'_'+exp[e_i]+'.nc','r')
        timnc = np.array(datanc.variables['time'])
        timncconv = pd.to_datetime(timnc-719529,unit='D')
        # find time that matches in exp and compare arrays
        # cntrl
        if compexp == 'cntrl':
            res = [key for key, val in enumerate(timncconv) if val in timcompconv]
        # loads1617 
        else:
            res = [key for key, val in enumerate(timcompconv) if val in timncconv]
        dttnc = np.array(datanc.variables['dt'])
        volnc = np.array(datanc.variables['volume'])
        arenc = np.array(datanc.variables['area'])
        bgcnc = (np.array(datanc.variables['bgc'])*s2d)/(dttnc*arenc)
        invnc = np.array(datanc.variables['invm'])/volnc
        # take difference (of only indices that match time period)
        # cntrl
        if compexp=='cntrl':
            invdif = invnc[res]
        # loads1617
        else:
            invdif = invnc 
        # write to array
        inv_arr[e_i,:,dep_ind] = invdif

    dep_ind += 1

figw = 14
figh = 10
axisfont = 16

width = 0.2

c_map = cmocean.cm.balance

for e_i in range(len(exp)):

    fig,ax = plt.subplots(1,1,figsize=[figw,figh])
    # cntrl
    if compexp=='cntrl':
        p_plt = ax.pcolormesh(timcompconv,list(range(dp0,dp1,dpstp)),np.transpose(inv_arr[e_i]),cmap=c_map,norm=mcolors.DivergingNorm(0),vmin=v_min,vmax=v_max)
        ax.set_title(varstr+' '+title_exp[e_i]+' - CTRL',fontsize=axisfont)
    #loads1617 
    else:
        p_plt = ax.pcolormesh(timncconv,list(range(dp0,dp1,dpstp)),np.transpose(inv_arr[e_i]),cmap=c_map,norm=mcolors.DivergingNorm(0),vmin=v_min,vmax=v_max)
        ax.set_title(varstr+' '+title_exp[e_i]+' - Loads 16-17',fontsize=axisfont)
    ax.invert_yaxis()
    ax.set_ylabel('Depth',fontsize=axisfont)
    ax.tick_params(axis='both',which='major',labelsize=axisfont)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    p0 = ax.get_position().get_points().flatten()
    p1 = ax.get_position().get_points().flatten()
    cb_ax = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])

    cb = fig.colorbar(p_plt,cax=cb_ax,orientation='vertical')
    cb.set_label(cblabel,fontsize=axisfont)
    cb.ax.tick_params(axis='both',which='major',labelsize=axisfont)

    fig.savefig(savepath+savename+'_'+exp[e_i]+'_'+compexp+'.png',bbox_inches='tight')

