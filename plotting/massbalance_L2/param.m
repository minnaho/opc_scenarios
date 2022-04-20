%%%% NO MODIFICATION
tname = 'ocean_time'; % for original files (on sigma levels)
displaynumbers =0 ; % display diagnostics
k=1; % index of the mask
Simu=2; % 1 or 4 or 2
region= 1; % dont change
    theta_s = 6.0;
    theta_b = 3.0;
    hc = 250;
    sc_type = 'new2012'; % for my zlevs4!!
bio=1 ; % 1 to activate biogeochemical rates calculation
vname = 'O2';  % name of the variable, O2, Alk, DIC , DON, ORG for sum of particulate matter
%vname = 'C' ;
%vname = 'NO3' ;
%vname = 'NH4' ;



%%%%%%%%%%%%%% USER CAN modify HERE ONLY
rep_out = './';
%mkdir (rep_out)
addpath(rep_out)

msk=10  % YOUR MASK , goes from 1 to 10: "10 is full bight"
%% READ FILES
%% MODEL OUTPUTS DIRECTOR
if scenario==1
repstr = 'loads1617';
rep_mod  = ['/data/project6/ROMS/L2SCB_OPC/',repstr,'/monthly/'] ;
%repstr = 'cntrl';
%rep_mod  = ['/data/project3/minnaho/postprocessing/cntrl_monthly'] ;
rep_his    = rep_mod ;
rep_out = [rep_out,'/',repstr,'/'];
mkdir(rep_out)
elseif scenario==2
repstr = 'l1617';
rep_mod  = ['/data/project6/ROMS/L2SCB_OPC/',repstr,'/monthly/'] ;
rep_his  = rep_mod ;
rep_out = [rep_out,'/',repstr,'/'];
mkdir(rep_out)
end
        repavg = dir([rep_mod,'/l2_scb_avg.Y*.nc']) ;
        repp = dir([rep_mod,'/l2_scb_phys_flux.*.nc']) ;
        repb = dir([rep_mod,'/l2_scb_bgc_flux_avg.*.nc']) ;
        rstrep = dir([rep_mod,'/l2_scb_his.*.nc']) ;

                final=length(repavg) ;

