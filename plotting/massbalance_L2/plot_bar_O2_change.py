import numpy as np
from netCDF4 import Dataset,num2date
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
import cmocean

#plt.ion()

varstr = 'O2'
cblabel = 'mmol/m3'
dp0 = 0
dp1 = 320
dpstp = 20
maskn = 10

fname = 'budget_L2_mask'+str(maskn)+'_'+varstr+'_'

# experiment to compare to
compexp = 'loads1617'
#compexp = 'cntrl'

if compexp == 'cntrl':
    exp = ['loads1617','PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
    title_exp = ['Loads 16-17','PNDN only','PNDN 50','PNDN 90','FNDN only','FNDN 50','FNDN 90']
else:
    exp = ['PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
    title_exp = ['PNDN only','PNDN 50','PNDN 90','FNDN only','FNDN 50','FNDN 90']


ncpath = '/data/project3/minnaho/opc_scenarios/plotting/massbalance_L2/'

savepath = '/data/project3/minnaho/opc_scenarios/plotting/massbalance_L2/figs/'
savename = 'hov_budget_mask'+str(maskn)+'_'+varstr

s2d = 86400

# cntrl
if compexp == 'cntrl':
    # read first file to get size of arrays
    datacomp = Dataset(ncpath+compexp+'/'+fname+str(0)+'_to_'+str(20)+'_'+compexp+'.nc','r')
    
# l1617
else:
    datacomp = Dataset(ncpath+exp[0]+'/'+fname+str(0)+'_to_'+str(20)+'_'+exp[0]+'.nc','r')
    
timcomp = np.array(datacomp.variables['time'])
datetimecomp = pd.to_datetime(timcomp-719529,unit='D')
# dimensions are time and depth (3rd dimension will be variable
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
            invdif = invnc[res] - invcomp
        # loads1617
        else:
            invdif = invnc - invcomp[res]
        # write to array
        inv_arr[e_i,:,dep_ind] = invdif

    dep_ind += 1

figw = 14
figh = 10
axisfont = 16

width = 0.2

fig,ax = plt.subplots(1,1,figsize=[figw,figh])
x_ind = np.arange(len(exp))

ax.bar(x_ind,np.nansum(np.nansum(inv_arr,axis=1),axis=1)/(inv_arr.shape[1]*inv_arr.shape[2]),width=width)
ax.set_xticks([width,1+width,2+width,3+width,4+width,5+width,6+width])
ax.set_xticklabels(title_exp)

#ax.scatter(x_ind,np.nansum(np.nansum(inv_arr,axis=1),axis=1)/(inv_arr.shape[1]*inv_arr.shape[2]))
#ax.set_xticks(range(len(exp)))
#ax.set_xticklabels(title_exp)


ax.tick_params(axis='both',which='major',labelsize=axisfont)
ax.set_ylabel(varstr+' mmol/m3 change',fontsize=axisfont)
fig.savefig('bar_O2_change'+'_'+compexp+'.png',bbox_inches='tight')
#fig.savefig('scatter_O2_change'+'_'+compexp+'.png',bbox_inches='tight')

