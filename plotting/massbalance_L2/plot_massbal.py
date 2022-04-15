import numpy as np
from netCDF4 import Dataset,num2date
import matplotlib.pyplot as plt

varstr = 'O2'
dp0 = 0
dp1 = 320
dpstp = 20
maskn = 10

fname = 'budget_L2_mask'+str(maskn)+'_'+varstr+'_'

exp = ['loads1617','PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
title_exp = ['Loads 16-17','PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']

# experiment to compare to
compexp = 'cntrl'

ncpath = '/data/project3/minnaho/opc_scenarios/plotting/massbalance_L2/'

savepath = '/data/project3/minnaho/opc_scenarios/plotting/massbalance_L2/figs'

s2d = 86400

# read first file to get size of arrays
datacomp = Dataset(ncpath+compexp+'/'+fname+str(0)+'_to_'+str(20)+'_'+compexp+'.nc','r')
timcomp = np.array(datacomp.variables['time'])

# dimensions are time and depth (3rd dimension will be variable
# for DIN/TN calculation
compbgc_arr = np.ones((len(timcomp),range(dp0,dp1,dpstp)))

for d_i in range(dp0,dp1,dpstp):
    # exp to compare to
    datacomp = Dataset(ncpath+compexp+'/'+fname+str(d_i)+'_to_'+str(d_i+dpstp)+'_'+compexp+'.nc','r')
    timcomp = np.array(datacomp.variables['time'])
    dttcomp = np.array(datacomp.variables['dt'])
    volcomp = np.array(datacomp.variables['volume'])
    arecomp = np.array(datacomp.variables['area'])
    bgccomp = np.array(datacomp.variables['bgc'])*s2d*dttcomp*arecomp
    invcomp = np.array(datacomp.variables['invm'])*volcomp
    for e_i in range(len(exp)):
        datanc = Dataset(ncpath+exp[e_i]+'/'+fname+str(d_i)+'_to_'+str(d_i+dpstp)+'_'+exp[e_i]+'.nc','r')
        timnc = np.array(datanc.variables['time'])
        dttnc = np.array(datanc.variables['dt'])
        volnc = np.array(datanc.variables['volume'])
        arenc = np.array(datanc.variables['area'])
        bgcnc = np.array(datanc.variables['bgc'])*s2d*dttnc*arenc
        invnc = np.array(datanc.variables['invm'])*volnc


        


