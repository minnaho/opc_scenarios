addpath(genpath('/data/project3/kesf/tools_matlab/matlab_paths/'))
%call_ncviewcolors
% carbonate system variables - 
% pH,omega from CO2SYS and omega from Juranek

depth = 50; % depth to slice
var_name1 = 'omega_juranek'
var_name2 = 'omega_co2sys'
var_name3 = 'pH'
%var_name = 'biomass'
% freshw
%rep = '/data/project6/ROMS/L2SCB_AP/fresh/daily/';
%str = 'fresh'
% nutrient
%rep = '/data/project5/kesf/ROMS/L2SCB_AP/nutrients/daily/';
% full
%rep = '/data/project6/kesf/ROMS/L2SCB_AP/daily/';
% control
%rep = '/data/project5/kesf/ROMS/L2_SCB/DAILY/';
%str = 'cntrl'
% pipes only
rep = '/data/project6/ROMS/L2SCB_P_1999_2000/daily/';
str = 'pipes'

%dt = 'Y2000M10'


%% LOAD GRID
Simu =2 ; % 1 for L1 , 0 for L0, 2 for L2-SCB
[pm pn lon_rho lat_rho lon_psi lon_psi f mask_rho h angle NY NX NZ] = loadgrid(Simu);
    theta_s = 6.0;
    theta_b = 3.0;
    hc = 250;
    sc_type = 'new2012'; % for my zlevs4!!

fout1 = ['sli_50m_',str,'_',var_name1,'.nc'];
fout2 = ['sli_50m_',str,'_',var_name2,'.nc'];
fout3 = ['sli_50m_',str,'_',var_name3,'.nc'];
%fout1 = ['sli_50m_',str,'_',var_name1,'_',dt,'.nc'];
%fout2 = ['sli_50m_',str,'_',var_name2,'_',dt,'.nc'];
%fout3 = ['sli_50m_',str,'_',var_name3,'_',dt,'.nc'];
shortname1 = ['slice 50m ',var_name1] ;
longname1  = ['slice 50m ',var_name1] ;

shortname2 = ['slice 50m ',var_name2] ;
longname2  = ['slice 50m ',var_name2] ;

shortname3 = ['slice 50m ',var_name3] ;
longname3  = ['slice 50m ',var_name3] ;
unit = 'mmol/m2';


        create_netcdf3D_L2(fout1,var_name1,shortname1,longname1,unit)
        create_netcdf3D_L2(fout2,var_name2,shortname2,longname2,unit)
        create_netcdf3D_L2(fout3,var_name3,shortname3,longname3,unit)

repavg = dir([rep,'l2_scb_avg.Y*.nc']) ;
%repavg = dir([rep,'l2_scb_avg.',dt,'D*.nc']) ;


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
   %dataout  = ncread(file, 'DIATC') + ncread(file, 'SPC') + ncread(file, 'DIAZC');

   %dataout  = ncread(file, var_name);
   %dataout = permute(dataout, [3 2 1]);
   extract_pH_omega

   %VAR = squeeze(VAR1(:,:));

ncwrite(fout1, var_name1, om_juranek' , [1 1 cpt]);
ncwrite(fout2, var_name2, om' , [1 1 cpt]);
ncwrite(fout3, var_name3, pH' , [1 1 cpt]);

   cpt = cpt+1
%  end % kk

end % fr


