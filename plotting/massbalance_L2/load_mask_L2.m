addpath(genpath('/data/project3/kesf/tools_matlab/matlab_paths/'))
%% run masks
load mask_gridL2.mat

mask0 = double(mask_gridL2.mask0) ;
mask1 = double(mask_gridL2.mask1) ;
mask2 = double(mask_gridL2.mask2) ;
mask3 = double(mask_gridL2.mask3) ;
mask4 = double(mask_gridL2.mask4) ;
mask5 = double(mask_gridL2.mask5) ;
mask6 = double(mask_gridL2.mask6) ;
mask7 = double(mask_gridL2.mask7) ;
mask8 = double(mask_gridL2.mask8) ;
mask9 = double(mask_gridL2.mask9) ;


disp('read masks..>> done')

return

%% LOAD GRID
Simu =2 ; % 1 for L1 , 0 for L0, 2 for L2-SCB
[pm pn lon_rho lat_rho lon_psi lon_psi f mask_rho h angle NY NX NZ] = loadgrid(Simu);
Simu =1 ; % 1 for L1 , 0 for L0, 2 for L2-SCB
[pm pn lon_rho1 lat_rho1 lon_psi lon_psi f mask_rho h angle NY NX NZ] = loadgrid(Simu);
mask0_L1 = griddata(lon_rho,lat_rho,mask0,lon_rho1,lat_rho1);
mask1_L1 = griddata(lon_rho,lat_rho,mask1,lon_rho1,lat_rho1);
mask2_L1 = griddata(lon_rho,lat_rho,mask2,lon_rho1,lat_rho1);
mask3_L1 = griddata(lon_rho,lat_rho,mask3,lon_rho1,lat_rho1);
mask4_L1 = griddata(lon_rho,lat_rho,mask4,lon_rho1,lat_rho1);
mask5_L1 = griddata(lon_rho,lat_rho,mask5,lon_rho1,lat_rho1);
mask6_L1 = griddata(lon_rho,lat_rho,mask6,lon_rho1,lat_rho1);
mask7_L1 = griddata(lon_rho,lat_rho,mask7,lon_rho1,lat_rho1);
mask8_L1 = griddata(lon_rho,lat_rho,mask8,lon_rho1,lat_rho1);
mask9_L1 = griddata(lon_rho,lat_rho,mask9,lon_rho1,lat_rho1);


