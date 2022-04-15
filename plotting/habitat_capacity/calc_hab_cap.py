import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import l2grid as l2grid
import numpy as np
from netCDF4 import Dataset


# choose years
start_year = 1997
end_year = 1999

# choose months between 1 and 12
start_month = 11
end_month = 11

omth = 1.4

# scenario names
#exp = ['l1617','PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
#title_exp = ['Loads 16-17','PNDN only','PNDN 50','PNDN 90','FNDN only','FNDN 50','FNDN 90']
exp = ['l1617','PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
#title_exp = ['cntrl','PNDN only','PNDN 50','PNDN 90','FNDN only','FNDN 50','FNDN 90']
#exp = ['PNDN_only','pndn50','pndn90']
#title_exp = ['PNDN only','PNDN 50','PNDN 90']
#exp = ['FNDN_only','fndn50','fndn90']
#title_exp = ['FNDN only','FNDN 50','FNDN 90']

ncpath = '/data/project6/friederc/data_products/opc_scenarios/'
fst = 'omega_th_'

# scenario to compare
compstr = 'cntrl'

mask_nc = l2grid.mask_nc
pm_nc = l2grid.pm_nc
pn_nc = l2grid.pn_nc

xisize = 1E-3/pm_nc
etasize = 1E-3/pn_nc

# percentage of habitat capacity change to select as threshold
percnum = 10

nummon = 0
for y_i in range(start_year,end_year+1):
    # if we are on the first year, starts at s_m
    if y_i == start_year:
        s_m = start_month
    else:
        s_m = 1
    # if we are on the last year, end at e_m
    if y_i == end_year:
        e_m = end_month+1
    else:
        e_m = 13
    for m_i in range(s_m,e_m):
        dtstr = 'Y'+str(y_i)+'M'+'%02d'%m_i
        nummon += 1

# area of positive and negative percnum in km^2
#dimensions are experiment and number of time steps
areasumpos = np.ones((len(exp),nummon))*np.nan
areasumneg = np.ones((len(exp),nummon))*np.nan

dtnc = np.ones((nummon))*np.nan

n_i = 0
for y_i in range(start_year,end_year+1):
    # if we are on the first year, starts at s_m
    if y_i == start_year:
        s_m = start_month
    else:
        s_m = 1
    # if we are on the last year, end at e_m
    if y_i == end_year:
        e_m = end_month+1
    else:
        e_m = 13
    for m_i in range(s_m,e_m):
        dtstr = 'Y'+str(y_i)+'M'+'%02d'%m_i
        for e_i in range(len(exp)):
            print(exp[e_i],dtstr)
            # map of +/-10% habitat capacity to save
            #mapcomp = np.ones((mask_nc.shape[0],mask_nc.shape[1]))*np.nan
            # experiment
            datanc = Dataset(ncpath+fst+dtstr+'_'+exp[e_i]+'.nc','r')
            # find threshold index
            tharr = np.array(datanc.variables['threshold'])
            thind = np.where(tharr==omth)[0][0]
            habcap = np.squeeze(datanc.variables['habitat_th'])[thind,:,:]
            # compare
            datacomp = Dataset(ncpath+fst+dtstr+'_'+compstr+'.nc','r')
            habcomp = np.squeeze(datacomp.variables['habitat_th'])[thind,:,:]
            habperc = ((habcap-habcomp)/habcomp)*100
            percp10 = np.where(habperc>=percnum)
            percn10 = np.where(habperc<=-1*percnum)
            #mapcomp[percp10[0],percp10[1]] = 1
            #mapcomp[percn10[0],percn10[1]] = -1
            # get habitat capacity %
            mapcomp = habperc
            areapos = np.nansum((xisize[percp10[0],percp10[1]])*(etasize[percp10[0],percp10[1]]))
            areaneg = np.nansum((xisize[percn10[0],percn10[1]])*(etasize[percn10[0],percn10[1]]))

            # save map
            mapnc = Dataset('./maps/map_'+fst+str(omth)+'_'+dtstr+'_'+exp[e_i]+'_'+compstr+'.nc','w')
            mapnc.createDimension('eta',mask_nc.shape[0])
            mapnc.createDimension('xi',mask_nc.shape[1])
            varhab = mapnc.createVariable('habitat_cap',np.float32,('eta','xi'))
            varhab.description = 'habitat capacity change from '+compstr+' where positive is an increase and negative is decrease'
            varhab[:,:] = mapcomp
            mapnc.close()
        
            # save cumulative area
            areasumpos[e_i,n_i] = areapos
            areasumneg[e_i,n_i] = areaneg

        # increment month
        dtnc[n_i] = np.array(datanc.variables['time'])[0]
        n_i += 1

for e_i in range(len(exp)): 
    sumnc = Dataset('./area/area_'+fst+str(omth)+'_Y'+str(start_year)+'M'+'%02d'%start_month+'_'+'Y'+str(end_year)+'M'+'%02d'%end_month+'_'+exp[e_i]+'_'+compstr+'.nc','w')
    sumnc.createDimension('time',nummon)
    dtvar = sumnc.createVariable('time',np.float32,'time')
    ncpos = sumnc.createVariable('posarea',np.float64,'time')
    ncneg = sumnc.createVariable('negarea',np.float64,'time')
    dtvar[:] = dtnc
    ncpos[:] = areasumpos[e_i]
    ncneg[:] = areasumneg[e_i]
    dtvar.units = 'days since 1997-10-20'
    ncpos.units = 'km^2'
    ncneg.units = 'km^2'
    ncpos.description = 'area of increase of habitat capacity greater than '+str(percnum)+'%'
    ncneg.description = 'area of decrease of habitat capacity greater than -'+str(percnum)+'%'
    sumnc.close() 
            

            
            

