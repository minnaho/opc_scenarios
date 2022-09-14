% runpfd_chll.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run PDF (chll) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: PDF to compliment Minna's analysis. Plots dentisty, salinity, and chll 
% 
% POSEIDON
% Last modified: 9/6/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add paths to things
addpath(genpath('/data/project3/kesf/tools_matlab/matlab_paths/'))
addpath(genpath('/data/project1/minnaho/make_masks/'))
figuresdir = '/data/project5/paigehoel/matlabscripts/'; %Locations to save figures
tabledir = '/data/project5/paigehoel/matlabscripts/'; %Locations to save figures
tabledir_s = '/data/project5/paigehoel/matlabscripts/tables/'; %Locations to save figures

%% Model Locations
fresh = '/data/project6/ROMS/L2SCB_AP/fresh/daily/';
con = '/data/project5/kesf/ROMS/L2_SCB/daily/';
nutri = '/data/project6/ROMS/L2SCB_AP/nutrients/daily/' ;
outfall= '/data/project6/ROMS/L2SCB_P_1999_2000/daily/';
full = '/data/project6/kesf/ROMS/L2SCB_AP/daily/';
nitrate = '/data/project6/ROMS/L2SCB_AP/nitrate/daily/';

theta_s = 6; % streching parameters (surface)
theta_b = 3; % streching parameters (bottom)
hc = 250;
NZ = 60;
sc_type ='new2012';

%% LOAD GRID
grd = '/data/project5/kesf/ROMS/L2_SCB/organization/L2_SCB/grid_3/roms_grd.nc';
h = ncread(grd,'h')';
lon = ncread(grd,'lon_rho')';
lat = ncread(grd,'lat_rho')';
pm  = ncread(grd, 'pm')';
pn  = ncread(grd, 'pn')';
[NY, NX] = size(h) ;

%% Mask file
maskfile = '/data/project1/minnaho/make_masks/mask_scb.nc';
coast_area = ncread(maskfile, 'mask_coast'); % 15 km
coast_area = permute(coast_area, [2 1]);
bight_area = ncread(maskfile, 'mask_bight');
bight_area = permute(bight_area, [2 1]);
mask = coast_area;
mask_name = 'Coast';  
area = coast_area;

% numbers correspond to days within the daily file
fall  = 92:183; % Oct 1 to Dec 31
winter = 184:276; % Jan 1 to Mar 31
spring = 277:367; % Apr 1 to Jun 30
summer = 368:457; % Jul 1 to Sept 30
year = 92:457; % Oct 1 to Sept 30
season_dates = {fall winter spring summer year};
season_title = {'Fall','Winter','Spring', 'Summer','year'};

S = 5;
seas = season_dates{S}; 
season = char(season_title(S));

%% Nutrient model
run = 'Nutrient';
repd = nutri;
for n = seas 
repavg = dir([repd,'l2_scb_avg*.nc']) ;
avg = [repd,'/',repavg(n,1).name] ;
z = ncread(avg,'zeta')';

% convert z to actual height
[z_w,Cw] = zlevs4(h, z, theta_s, theta_b, hc, NZ, 'w',sc_type);
dz = diff(z_w);
     zbot = flipdim(cumsum(flipdim(dz,1)),1);
     ztop = [zbot(2:end,:,:);zeros(1,NY,NX)];
     z_rho = (zbot+ztop)./2 ;

% read all of the variables

var_Chll = ncread(avg,'SPCHL')+ncread(avg,'DIATCHL')+ncread(avg,'DIAZCHL');% mg Chl-a m-3
 var_Chll = permute(var_Chll, [3 2 1]);
VARChl = squeeze(var_Chll(60,:,:)); % extract surface value

file = VARChl;
l = histogram(file);
l.BinEdges = [0:.1:50]; % chll
counts_Chl(:,n) = l.Values;

end

Nutri = {counts_Chl};
disp('nutri done')


%% Outfall model
run = 'Outfall';
repd = outfall;
for n = seas 
repavg = dir([repd,'l2_scb_avg*.nc']) ;
avg = [repd,'/',repavg(n,1).name] ;
z = ncread(avg,'zeta')';

% convert z to actual height
[z_w,Cw] = zlevs4(h, z, theta_s, theta_b, hc, NZ, 'w',sc_type);
dz = diff(z_w);
     zbot = flipdim(cumsum(flipdim(dz,1)),1);
     ztop = [zbot(2:end,:,:);zeros(1,NY,NX)];
     z_rho = (zbot+ztop)./2 ;

% read all of the variables

var_Chll = ncread(avg,'SPCHL')+ncread(avg,'DIATCHL')+ncread(avg,'DIAZCHL');
 var_Chll = permute(var_Chll, [3 2 1]);
VARChl = squeeze(var_Chll(60,:,:));
a = max(VARChl,[],'all')

file = VARChl;
l = histogram(file);
l.BinEdges = [0:.1:50]; % chll
counts_Chl(:,n) = l.Values;

end
% take 60
Outfall = {counts_Chl};
disp('outfall done')

%% Fresh model
run = 'Fresh';
repd = fresh;
for n = seas 
repavg = dir([repd,'l2_scb_avg*.nc']) ;
avg = [repd,'/',repavg(n,1).name] ;
z = ncread(avg,'zeta')';

% convert z to actual height
[z_w,Cw] = zlevs4(h, z, theta_s, theta_b, hc, NZ, 'w',sc_type);
dz = diff(z_w);
     zbot = flipdim(cumsum(flipdim(dz,1)),1);
     ztop = [zbot(2:end,:,:);zeros(1,NY,NX)];
     z_rho = (zbot+ztop)./2 ;

% read all of the variables
var_Chll = ncread(avg,'SPCHL')+ncread(avg,'DIATCHL')+ncread(avg,'DIAZCHL');
 var_Chll = permute(var_Chll, [3 2 1]);
VARChl = squeeze(var_Chll(60,:,:));

file = VARChl;
l = histogram(file);
l.BinEdges = [0:.1:50]; % chll
counts_Chl(:,n) = l.Values;
end

Fresh = {counts_Chl};
disp('fresh done')

%% Nitrate model
run = 'Nitrate';
repd = nitrate;
for n = seas 
repavg = dir([repd,'l2_scb_avg*.nc']) ;
avg = [repd,'/',repavg(n,1).name] ;
z = ncread(avg,'zeta')';

% convert z to actual height
[z_w,Cw] = zlevs4(h, z, theta_s, theta_b, hc, NZ, 'w',sc_type);
dz = diff(z_w);
     zbot = flipdim(cumsum(flipdim(dz,1)),1);
     ztop = [zbot(2:end,:,:);zeros(1,NY,NX)];
     z_rho = (zbot+ztop)./2 ;

% read all of the variables
var_Chll = ncread(avg,'SPCHL')+ncread(avg,'DIATCHL')+ncread(avg,'DIAZCHL');
 var_Chll = permute(var_Chll, [3 2 1]);
VARChl = squeeze(var_Chll(60,:,:));

file = VARChl;
l = histogram(file);
l.BinEdges = [0:.1:50]; % chll
counts_Chl(:,n) = l.Values;
end

Nitrate = {counts_Chl};
disp('nitrate done')

%% Full model
run = 'Full';
repd = full;
for n = seas 
repavg = dir([repd,'l2_scb_avg*.nc']) ;
avg = [repd,'/',repavg(912+n,1).name] ;
z = ncread(avg,'zeta')';

% convert z to actual height
[z_w,Cw] = zlevs4(h, z, theta_s, theta_b, hc, NZ, 'w',sc_type);
dz = diff(z_w);
     zbot = flipdim(cumsum(flipdim(dz,1)),1);
     ztop = [zbot(2:end,:,:);zeros(1,NY,NX)];
     z_rho = (zbot+ztop)./2 ;

% read all of the variables

var_Chll = ncread(avg,'SPCHL')+ncread(avg,'DIATCHL')+ncread(avg,'DIAZCHL');
 var_Chll = permute(var_Chll, [3 2 1]);
VARChl = squeeze(var_Chll(60,:,:));
a = max(VARChl,[],'all')


file = VARChl;
l = histogram(file);
l.BinEdges = [0:.1:50]; % chll
counts_Chl(:,n) = l.Values;

end

Full = {counts_Chl};
disp('full done')

%% Control model
run = 'Control';
repd = con;
for n = seas 
repavg = dir([repd,'l2_scb_avg*.nc']) ;
avg = [repd,'/',repavg((n+912),1).name] ;
z = ncread(avg,'zeta')';

% convert z to actual height
[z_w,Cw] = zlevs4(h, z, theta_s, theta_b, hc, NZ, 'w',sc_type);
dz = diff(z_w);
     zbot = flipdim(cumsum(flipdim(dz,1)),1);
     ztop = [zbot(2:end,:,:);zeros(1,NY,NX)];
     z_rho = (zbot+ztop)./2 ;

% read all of the variables
var_Chll = ncread(avg,'SPCHL')+ncread(avg,'DIATCHL')+ncread(avg,'DIAZCHL');
 var_Chll = permute(var_Chll, [3 2 1]);
VARChl = squeeze(var_Chll(60,:,:));

file = VARChl;
l = histogram(file);
l.BinEdges = [0:.1:50]; % chll
counts_Chl(:,n) = l.Values;
end
Control = {counts_Chl};
disp('control done')

names = {'Chl'};

for t = 1;
variable = char(names(t));
% Save value table
T{t} = [Control(:,t) Fresh(:,t) Nutri(:,t) Nitrate(:,t) Outfall(:,t) Full(:,t)];
T_new = array2table(T{t},'VariableNames', {'Control','Fresh','Nutrient','Nitrate','Outfall','Full'});
output_name = [tabledir, variable, season]; %change with season loop
writetable(T_new, output_name);
clear T T_new d D_new
end %% table save

disp('now do an actual pdf and make some figures!')

%save('pdf_minna_surf_chll.mat', 'Control','Fresh','Nutri','Nitrate','Outfall','Full')
%
%keyboard
%
%var = 1;
%
%file = Nutri(var);
%mu = mean(file);
%sigma = std(file);
%pdf_Nutri = pdf('Normal',file,mu,sigma);
%
%file = Nitrate(var);
%mu = mean(file);
%sigma = std(file);
%pdf_Nitrate = pdf('Normal',file,mu,sigma);
%
%file = Outfall(var);
%mu = mean(file);
%sigma = std(file);
%pdf_Outfall = pdf('Normal',file,mu,sigma);
%
%file = Control(var);
%mu = mean(file);
%sigma = std(file);
%pdf_Control = pdf('Normal',file,mu,sigma);
%
%file = Fresh(var);
%mu = mean(file);
%sigma = std(file);
%pdf_Fresh = pdf('Normal',file,mu,sigma);
%
%figure
%
%return
