import matplotlib.pyplot as plt
import numpy as np

plt.ion()

#exp = ['l1617','PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
#title_exp = ['Loads 16-17','PNDN only','PNDN 50','PNDN 90','FNDN only','FNDN 50','FNDN 90']
exp = ['PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
title_exp = ['PNDN only','PNDN 50','PNDN 90','FNDN only','FNDN 50','FNDN 90']

figw = 14
figh = 10

axis_font = 16

# compared to control
habcom = [502299,504918,490109,490797,503085,460250,492345]
habexp = [262559,271301,289600,276483,271968,288772,274239]

# compared to loads 16-17
habcom = [97410,104335,88815,95996,52950,92681]
habexp = [92178,91447,100648,93272,104580,76557]


x_ind = np.arange(len(exp))

width = 0.2
delx = 0.01

fig,ax = plt.subplots(1,1,figsize=[figw,figh])

ax.bar(x_ind,habcom,color='red',width=width)
ax.bar(x_ind+width,habexp,color='blue',width=width)
#ax.set_xticks([width+delx,1+width+delx,2+width+delx,3+width+delx,4+width+delx,5+width+delx,6+width+delx])
ax.set_xticks([width+delx,1+width+delx,2+width+delx,3+width+delx,4+width+delx,5+width+delx])
ax.set_xticklabels(title_exp)
ax.tick_params(axis='both',which='major',labelsize=axis_font)
ax.set_ylabel('km$^2$',fontsize=axis_font)
fig.savefig('bar_hab_comp_exp_l1617.png',bbox_inches='tight')


