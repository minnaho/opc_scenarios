import numpy as np
import matplotlib.pyplot as plt

mean = np.array([574.9028,508.6095,722.9948,627.3651,773.1115,422.7487])
eddy = np.array([4.2550E3,4.2487e+03,4.1227e+03,3.0296e+03,2.8327e+03,2.5984e+03])

varstr = 'NH4'

title_exp = ['CTRL',
             'ANTH',
             '50% N Red.',
             '50% N Red.\n50% Recy.',
             '50% N Red.\n90% Recy.',
             '85% N Red.',
             '85% N Red.\n50% Recy.',
             '85% N Red.\n90% Recy.']

plt.ion()

axfont = 14

fig,ax = plt.subplots(2,1,figsize=[8,5])
ax.flat[0].bar(range(len(title_exp[2:])),mean)
ax.flat[1].bar(title_exp[2:],eddy)
ax.flat[0].set_title('NH4 mean transport')
ax.flat[1].set_title('NH4 eddy transport')
fig.supylabel(varstr+' mmol m$^{-2}$ d$^{-1}$',fontsize=axfont)
ax.tick_params(axis='both',which='major',labelsize=axfont)
ax.flat[0].set_xticks([])
plt.tight_layout()
fig.savefig('./figs/2years/mean_eddy_bars.png',bbox_inches='tight')

# broken axis 
# https://matplotlib.org/stable/gallery/subplots_axes_and_figures/broken_axis.html

fig,ax = plt.subplots(1,1,figsize=[7,5])
ax.bar(range(len(title_exp[2:])),eddy,label='Eddy')
ax.bar(title_exp[2:],mean,bottom=eddy,label='Mean')
ax.set_title('NH4 transport')
ax.set_ylabel(varstr+' mmol m$^{-2}$ d$^{-1}$',fontsize=axfont)
ax.tick_params(axis='both',which='major',labelsize=axfont)
ax.set_xticks([])
ax.legend(loc='best',fontsize=axfont)
ax.set_yscale('linear')
ax.set_ylim([2000,5000])
#plt.tight_layout()
fig.savefig('./figs/2years/mean_eddy_sum_bars.png',bbox_inches='tight')
