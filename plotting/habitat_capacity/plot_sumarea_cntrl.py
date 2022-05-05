import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset,num2date

plt.ion()

path = '/data/project3/minnaho/opc_scenarios/plotting/habitat_capacity/area/'

#exp = ['l1617','PNDN_only','pndn50','pndn90']
#title_exp = ['Loads 16-17','PNDN only','PNDN 50','PNDN 90']
#exp = ['l1617','FNDN_only','fndn50','fndn90']
#title_exp = ['Loads 16-17','FNDN only','FNDN 50','FNDN 90']
exp = ['l1617','PNDN_only','FNDN_only']
title_exp = ['Loads 16-17','PNDN only','FNDN only']

#exp = ['l1617'] # compare loads 16-17 to cntrl
#title_exp = ['Loads 16-17']

dtstr = 'Y1997M11_Y1999M11'
omth = 1.4

fst = 'area_omega_th_'+str(omth)+'_'+dtstr+'_'

axisfont = 16

figw = 14
figh = 10

fig,ax = plt.subplots(1,1,figsize=[figw,figh])

lsty = ['-','-','--',':']
cpos = 'blue'
cneg = 'red'

for e_i in range(len(exp)):
    datanc = Dataset(path+fst+exp[e_i]+'_cntrl.nc','r')
    timevar = np.array(datanc.variables['time'])
    timeunit = datanc.variables['time'].units
    dt = num2date(timevar,timeunit,only_use_cftime_datetimes=False,only_use_python_datetimes=True)
    posnc = np.array(datanc.variables['posarea'])
    negnc = np.array(datanc.variables['negarea'])
    if exp[e_i] == 'l1617':
        ax.plot(dt,posnc,linewidth=3,color='green',linestyle=lsty[e_i],label=title_exp[e_i])
        ax.plot(dt,negnc,linewidth=3,color='purple',linestyle=lsty[e_i],label=title_exp[e_i])
    else:
        ax.plot(dt,posnc,color=cpos,linestyle=lsty[e_i],label=title_exp[e_i])
        ax.plot(dt,negnc,color=cneg,linestyle=lsty[e_i],label=title_exp[e_i])
    print(title_exp[e_i]+' total habitat expansion '+str(np.nansum(posnc)))
    print(title_exp[e_i]+' total habitat compression '+str(np.nansum(negnc)))
    print(title_exp[e_i]+' 1998 max expansion '+str(np.nanmax(posnc[2:14])))
    print(title_exp[e_i]+' 1998 max compression '+str(np.nanmax(negnc[2:14])))
    print(title_exp[e_i]+' 1999 max expansion '+str(np.nanmax(posnc[14:])))
    print(title_exp[e_i]+' 1999 max compression '+str(np.nanmax(negnc[14:])))

ax.set_ylim([0,60000])
ax.set_ylabel('km$^2$',fontsize=axisfont)
ax.tick_params(axis='both',which='major',labelsize=axisfont)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')

ax.legend(fontsize=axisfont)
fig.savefig('./figs/ts/'+fst+exp[-1]+'_cntrl.png',bbox_inches='tight')
    
    
