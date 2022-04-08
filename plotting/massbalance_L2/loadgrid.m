function [pm pn lon_rho lat_rho lon_psi lat_psi f mask_rho h angle NY NX NZ] = loadgrid(Simu); 
if Simu==1
grd = '/data/project3/kesf/ROMS/USSW1/grid/roms_grd.nc';
elseif Simu==0
grd = '/data/project4/kesf/ROMS/USW4/grid/roms_grd.nc';
elseif Simu==2
grd = '/data/project5/kesf/ROMS/L2_SCB/organization/L2_SCB/grid_3/roms_grd.nc';
elseif Simu==3
grd = '/data/project4/kesf/ROMS/L3_LAOC/grid/roms_grd.nc';
elseif Simu==4
grd = '/data/project4/kesf/ROMS/L4_OC/grid/roms_grd.nc';
elseif Simu==5
grd = '/data/project4/kesf/ROMS/L5_OC/grid/roms_grd.nc';
elseif Simu==11
grd= '/data/project6/kesf/ROMS/USNW1/organization/Other_files/usw1_grd.nc' ;
elseif Simu==7
%grd= '/data/project3/kesf/ROMS/L2_SFM/grid/sf03_grd.nc';
grd= '/data/project3/kesf/ROMS/L2_SFM/grid/roms_grd.nc';
%grd= '/data/project3/kesf/ROMS/L2_SFM/grid/dfm_grd.nc';
end

    in_dir = './';
    idx = 0;

pm  = ncread(grd, 'pm')';
pn  = ncread(grd, 'pn')';
lon_rho = ncread(grd, 'lon_rho')';
lat_rho = ncread(grd, 'lat_rho')';
lon_psi = ncread(grd, 'lon_psi')';
lat_psi = ncread(grd, 'lat_psi')';
f = ncread(grd,'f')';
angle = ncread(grd,'angle')';
[NY,NX] = size(lon_rho);
NZ=60 ;

mask_rho= ncread(grd, 'mask_rho')';
h   = ncread(grd, 'h')';
