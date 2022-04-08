%if ~exist(fout)
% create the ncfile
ncid = netcdf.create(fout,'CLOBBER');
% extract the dimensions
  dimtime = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED')) ;
% create the variables

%% time vector
D    =  netcdf.defVar(ncid,'time', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','days');
netcdf.putAtt(ncid,D,'long_name','date at the calculation of the flux (days from 0000-00-00)');

D    =  netcdf.defVar(ncid,'dt', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','seconds');
netcdf.putAtt(ncid,D,'long_name','time step');

%% volume
D    =  netcdf.defVar(ncid,'volume', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','m3');
netcdf.putAtt(ncid,D,'long_name','volume of the region');

D    =  netcdf.defVar(ncid,'area', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','m2');
netcdf.putAtt(ncid,D,'long_name','area of the region');

%% total horizontal fluxes
D    =  netcdf.defVar(ncid,'horizontal_flux', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol');
netcdf.putAtt(ncid,D,'long_name','total horizontal flux');

D    =  netcdf.defVar(ncid,'zonaln', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol');
netcdf.putAtt(ncid,D,'long_name','E-W zonal flux');

D    =  netcdf.defVar(ncid,'zonals', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol');
netcdf.putAtt(ncid,D,'long_name','E-W zonal flux');

D    =  netcdf.defVar(ncid,'zonal', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol');
netcdf.putAtt(ncid,D,'long_name','E-W zonal flux');

D    =  netcdf.defVar(ncid,'meridional', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol');
netcdf.putAtt(ncid,D,'long_name','S-N meriodional flux');

D    =  netcdf.defVar(ncid,'east', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol');
netcdf.putAtt(ncid,D,'long_name','E-flux');

D    =  netcdf.defVar(ncid,'west', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol');
netcdf.putAtt(ncid,D,'long_name','W-flux');

D    =  netcdf.defVar(ncid,'north', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol');
netcdf.putAtt(ncid,D,'long_name','N-flux');

D    =  netcdf.defVar(ncid,'south', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol');
netcdf.putAtt(ncid,D,'long_name','S-flux');

%% vertical flux
D    =  netcdf.defVar(ncid,'vertical_flux', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol');
netcdf.putAtt(ncid,D,'long_name','total vertical flux');

D    =  netcdf.defVar(ncid,'adv_vert_flux', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol');
netcdf.putAtt(ncid,D,'long_name','total advective vertical flux');

D    =  netcdf.defVar(ncid,'diffusive_flux', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol');
netcdf.putAtt(ncid,D,'long_name','total diffusive vertical flux');


%% dNdT
D    =  netcdf.defVar(ncid,'dNdT', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol');
netcdf.putAtt(ncid,D,'long_name','dN/dT between two time step (dt, usually = 1 month)');

D    =  netcdf.defVar(ncid,'inv1', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol');
netcdf.putAtt(ncid,D,'long_name','Inventory at 1');

D    =  netcdf.defVar(ncid,'inv2', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol');
netcdf.putAtt(ncid,D,'long_name','Inventory at 2');

D    =  netcdf.defVar(ncid,'invm', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol');
netcdf.putAtt(ncid,D,'long_name','mean concentration');

%% bgc
D    =  netcdf.defVar(ncid,'nitrification', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol');
netcdf.putAtt(ncid,D,'long_name','total nitrification between two time steps (dt, usually = 1 month)');

D    =  netcdf.defVar(ncid,'denitrification', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol');
netcdf.putAtt(ncid,D,'long_name','total denitrification between two time steps (dt, usually = 1 month)');

D    =  netcdf.defVar(ncid,'uptake_diat', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol');
netcdf.putAtt(ncid,D,'long_name','total diatom uptake between two time steps (dt, usually = 1 month)');

D    =  netcdf.defVar(ncid,'uptake_sp', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol');
netcdf.putAtt(ncid,D,'long_name','total small phytoplankton uptake between two time steps (dt, usually = 1 month)');

D    =  netcdf.defVar(ncid,'bgc', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol');
netcdf.putAtt(ncid,D,'long_name','total biogeochemical flux between two time steps (dt, usually = 1 month)');

D    =  netcdf.defVar(ncid,'meanNup', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol/m3');
netcdf.putAtt(ncid,D,'long_name','mean NO3 above 50m');

D    =  netcdf.defVar(ncid,'meanNdn', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmol/m3');
netcdf.putAtt(ncid,D,'long_name','mean NO3 below 50m');

D    =  netcdf.defVar(ncid,'NpPp', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmolN');
netcdf.putAtt(ncid,D,'long_name','NO3 prime Phyto prime');

D    =  netcdf.defVar(ncid,'PpZp', 'float', [dimtime]);
netcdf.putAtt(ncid,D,'units','mmolN');
netcdf.putAtt(ncid,D,'long_name','Phyto prime Zoo prime');

%% close the file
netcdf.endDef(ncid);
netcdf.close(ncid)

%end

