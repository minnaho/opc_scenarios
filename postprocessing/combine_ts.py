import numpy as np
from netCDF4 import Dataset,num2date,date2num

# Y2000M02D02 last time step of biomass pipe int 100 m
varn = 'var'
fout = 'int_100m_pipes_biomass.nc'
shortname = 'int 100m biomass' 

#m01nc = np.array(Dataset('int_100m_pipes_biomass_Y1999M07.nc','r').variables[varn])
#m02nc = np.array(Dataset('int_100m_pipes_biomass_Y1999M08.nc','r').variables[varn])
#m03nc = np.array(Dataset('int_100m_pipes_biomass_Y1999M09.nc','r').variables[varn])
#m04nc = np.array(Dataset('int_100m_pipes_biomass_Y1999M10.nc','r').variables[varn])
#m05nc = np.array(Dataset('int_100m_pipes_biomass_Y1999M11.nc','r').variables[varn])
#m06nc = np.array(Dataset('int_100m_pipes_biomass_Y1999M12.nc','r').variables[varn])
#m07nc = np.array(Dataset('int_100m_pipes_biomass_Y2000M01.nc','r').variables[varn])
m08nc = np.array(Dataset('int_100m_pipes_biomass_Y2000M02.nc','r').variables[varn])
m09nc = np.array(Dataset('int_100m_pipes_biomass_Y2000M03.nc','r').variables[varn])
m10nc = np.array(Dataset('int_100m_pipes_biomass_Y2000M04.nc','r').variables[varn])
m11nc = np.array(Dataset('int_100m_pipes_biomass_Y2000M05.nc','r').variables[varn])
m12nc = np.array(Dataset('int_100m_pipes_biomass_Y2000M06.nc','r').variables[varn])
m13nc = np.array(Dataset('int_100m_pipes_biomass_Y2000M07.nc','r').variables[varn])
m14nc = np.array(Dataset('int_100m_pipes_biomass_Y2000M08.nc','r').variables[varn])
m15nc = np.array(Dataset('int_100m_pipes_biomass_Y2000M09.nc','r').variables[varn])
#m16nc = np.array(Dataset('int_100m_pipes_biomass_Y2000M10.nc','r').variables[varn])

allnc = np.array(Dataset('int_100m_pipes_biomass_part.nc','r').variables['var'])

# ends at 2/1/2000, take back to 1/31/2000
ncslice = allnc[:-1,:,:]

fullnc0 = np.append(np.append(ncslice,m08nc,axis=0),m09nc,axis=0)
fullnc1 = np.append(np.append(fullnc0,m10nc,axis=0),m11nc,axis=0)
fullnc2 = np.append(np.append(fullnc1,m12nc,axis=0),m13nc,axis=0)
fullnc3 = np.append(np.append(fullnc2,m14nc,axis=0),m15nc,axis=0)
#fullnc4 = np.append(fullnc2,m16nc,axis=0)
#fullnc4 = np.append((np.append(np.append(fullnc3,m08nc,axis=0),m09nc,axis=0)),m10nc,axis=0)

#fullnc0 = np.append((np.append((np.append(np.append(m01nc,m02nc,axis=0),m03nc,axis=0)),m04nc,axis=0)),m05nc,axis=0)
#fullnc1 = np.append((np.append((np.append(np.append(fullnc0,m06nc,axis=0),m07nc,axis=0)),m08nc,axis=0)),m09nc,axis=0)
#fullnc2 = np.append((np.append((np.append(np.append(fullnc1,m10nc,axis=0),m11nc,axis=0)),m12nc,axis=0)),m13nc,axis=0)
#fullnc3 = np.append((np.append(np.append(fullnc2,m14nc,axis=0),m15nc,axis=0)),m16nc,axis=0)

outnc = Dataset(fout,'w')

timedim = outnc.createDimension('ocean_time',None)
latitude = outnc.createDimension('latitude',fullnc3.shape[1])
longitude = outnc.createDimension('longitude',fullnc3.shape[2])

var = outnc.createVariable(varn,'float32',('ocean_time','latitude','longitude'))
var.units = 'mmol/m2'
var.shortname = shortname
var.longname = 'days since 1999-07-01'

var[:,:,:] = fullnc3

outnc.close()
