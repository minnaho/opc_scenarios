% extract integrated variables from model

addpath(genpath('/data/project3/kesf/tools_matlab/matlab_paths/'))
%call_ncviewcolors

depth = 100; % depth to integrate from surface down to
%var_name = 'NO3'
var_name = 'biomass'
% freshw
%rep = '/data/project6/ROMS/L2SCB_AP/fresh/daily/';
%str = 'fresh'
% nutrient
%rep = '/data/project6/ROMS/L2SCB_AP/nutrients/daily/';
%str = 'nutri'
% full
rep = '/data/project3/minnaho/freshwater/postprocessing/fulll_daily/';
str = 'fulll'
% control
%rep = '/data/project3/minnaho/freshwater/postprocessing/cntrl_daily/';
%str = 'cntrl'
% pipes only
%rep = '/data/project6/ROMS/L2SCB_P_1999_2000/daily/';
%str = 'pipes'

%% LOAD GRID
Simu =2 ; % 1 for L1 , 0 for L0, 2 for L2-SCB
[pm pn lon_rho lat_rho lon_psi lon_psi f mask_rho h angle NY NX NZ] = loadgrid(Simu);
    theta_s = 6.0;
    theta_b = 3.0;
    hc = 250;
    sc_type = 'new2012'; % for my zlevs4!!

fout = ['int_',int2str(depth),'m_',str,'_',var_name,'.nc'];
ncvar = 'var';
shortname = ['integrated ',var_name] ;
longname = ['integrated ',var_name] ;
unit = 'mmol/m2';


        create_netcdf3D_L2(fout,ncvar,shortname,longname,unit)

repavg = dir([rep,'l2_scb_avg.Y1999*.nc']) ;

cpt = 1;
for fr = 1:length(repavg)
 disp(repavg(fr,1).name)
 file = [rep,'/',repavg(fr,1).name] 

zeta  = ncread(file, 'zeta')' ;
[z_w,Cw1] = zlevs4(h, zeta, theta_s, theta_b, hc, NZ, 'w',sc_type);
dz = diff(z_w);
        zbot = flipdim(cumsum(flipdim(dz,1)),1);
        ztop = [zbot(2:end,:,:);zeros(1,NY,NX)];
        z_rho = (zbot+ztop)./2 ;
   dataout  = ncread(file, 'DIATC') + ncread(file, 'SPC') + ncread(file, 'DIAZC');
   %dataout  = ncread(file, var_name);
   dataout = permute(dataout, [3 2 1]);

for i=1:depth
        VAR1(:,:,i)  = vinterp ( dataout, -(abs(z_rho)) ,  -abs(i) ) ;
end
   VAR = squeeze(nansum(VAR1(:,:,1:depth),3));

ncwrite(fout, 'var', VAR' , [1 1 cpt]);

   cpt = cpt+1
%  end % kk

end % fr


