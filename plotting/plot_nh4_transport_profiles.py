import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.ion()

figpath = './figs/2years/'

exp = [
       'PNDN_only_realistic',
       'pndn50_fixriver',
       'pndn90_fixriver',
       'FNDN_only_realistic',
       'fndn50_fixriver',
       'fndn90_fixriver']

title_exp = [
             '50% N Red.',
             '50% N Red.\n50% Recy.',
             '50% N Red.\n90% Recy.',
             '85% N Red.',
             '85% N Red.\n50% Recy.',
             '85% N Red.\n90% Recy.']

lstyle = ['-','--',':','-','--',':']
profcols = ['orange','orange','orange','gray','gray','gray']

axfont = 16

fig,ax = plt.subplots(1,2,figsize=[12,8],sharey=True)
for e_i in range(0,3):
    data_eddy = np.array(pd.read_fwf('eddy_'+exp[e_i]+'.txt'))
    data_mean = np.array(pd.read_fwf('mean_'+exp[e_i]+'.txt'))
    data_total = data_eddy+data_mean
    #data_total = data_eddy
    deparr = range(len(data_total))[::-1]
    ax.flat[0].plot(data_total[::2],deparr[::2],color=profcols[e_i],linestyle=lstyle[e_i],linewidth=2,label=title_exp[e_i])

#ax.flat[0].set_xlim([-30,50])
ax.flat[0].set_xlim([-5,60])
ax.flat[0].set_ylim([0,160])
ax.flat[0].invert_yaxis()
#ax.flat[0].legend(loc='lower left',fontsize=axfont-2)
ax.flat[0].legend(loc='center left',fontsize=axfont-2)
ax.flat[0].set_ylabel('Depth (m)',fontsize=axfont)
ax.flat[0].set_xlabel('mmol NH$_4^+$ m$^{-3}$ d$^{-1}$',fontsize=axfont)
ax.flat[0].tick_params(axis='both',which='major',labelsize=axfont)


for e_i in range(3,len(exp)):
    data_eddy = np.array(pd.read_fwf('eddy_'+exp[e_i]+'.txt'))
    data_mean = np.array(pd.read_fwf('mean_'+exp[e_i]+'.txt'))
    data_total = data_eddy+data_mean
    #data_total = data_eddy
    deparr = range(len(data_total))[::-1]
    ax.flat[1].plot(data_total[::2],deparr[::2],color=profcols[e_i],linestyle=lstyle[e_i],linewidth=2,label=title_exp[e_i])

#ax.flat[1].set_xlim([-30,50])
ax.flat[1].set_xlim([-5,60])
ax.flat[1].set_ylim([0,160])
ax.flat[1].invert_yaxis()
#ax.flat[1].legend(loc='lower left',fontsize=axfont-2)
ax.flat[1].legend(loc='center left',fontsize=axfont-2)
ax.flat[1].set_xlabel('mmol NH$_4^+$ m$^{-3}$ d$^{-1}$',fontsize=axfont)
ax.flat[1].tick_params(axis='both',which='major',labelsize=axfont)

fig.savefig(figpath+'nh4_flux_2panel.png',bbox_inches='tight')
