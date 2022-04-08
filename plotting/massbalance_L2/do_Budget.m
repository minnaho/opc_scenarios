%% =================================================================== %%
%% this program calculae a total budget of a biogeochemical tracer     %%
%% in a certain area of the coastal ocean			       %%
%% Faycal Kessouri - SCCWRP/UCLA (faycalk@sccwrp.org)	 - 	       %%
%% version 1 on 06/2017						       %%
%% modified on 05/2018						       %%
%% adapted for oxygen and ammonium on 04/2022			       %%
%% =================================================================== %%
addpath(genpath('/data/project3/kesf/tools_matlab/matlab_paths/'))

nbs=2 ; % NUMBER OF SCENARIOS
for scenario=1:nbs
param
%%%%%%%%%%%%%
%%% GRID FILES
%% load the grid
[pm pn lon lat lon_psi lon_psi f mask_rho h angle NY NX NZ] = loadgrid(Simu);
%%%%%%%%%%%%%
%% DO THE MASK

list1 = 0:20:300
%list1 = 0 ;
for dd = 1:length(list1)
depthmin = list1(dd) ; % shallower limit
depthmax = depthmin+20 ; % deeper limit

clear mask_box3dr mask_box mask_boxr indz indx indy
%domask
load_mask_L2

if msk==1
mask_box=mask1;
mask_boxr = mask1 ;
elseif msk==2
mask_box=mask2;
mask_boxr = mask2 ;
elseif msk==3
mask_box=mask3;
mask_boxr = mask3 ;
elseif msk==4
mask_box=mask4;
mask_boxr = mask4 ;
elseif msk==5
mask_box=mask5;
mask_boxr = mask5 ;
elseif msk==6
mask_box=mask6;
mask_boxr = mask6 ;
elseif msk==7
mask_box=mask7;
mask_boxr = mask7 ;
elseif msk==8
mask_box=mask8;
mask_boxr = mask8 ;
elseif msk==9
mask_box=mask9;
mask_boxr = mask9 ;
elseif msk==10
mask_box=(mask1.*0)+1; mask_box(mask_rho==0)=NaN;
mask_box(1:20,:) = NaN ;
mask_box(1392:1412,:) = NaN ;
mask_box(:,1:20) = NaN ;
mask_box(:,582:602) = NaN ;
mask_boxr = mask_box ;
end


mask_box3dr = repmat(mask_boxr,1,1,60);
mask_box3dr = permute(mask_box3dr, [3 1 2]);

fout =  [rep_out,'budget_L2_mask',num2str(msk,'%.2d'),'_',vname,'_',num2str(depthmin),'_to_',num2str(depthmax),'_',repstr,'.nc'];

%% CREATE OUTPUT FILES
if (strcmp(vname,'O2')==1)
filesOrg_O2
elseif (strcmp(vname,'NO3')==1)
filesOrg_nitrate
elseif (strcmp(vname,'NH4')==1)
filesOrg_ammonium
end




%filesOrg_nit_2D
%%%%%%%%%%%%%%
%% estimate distance for horizontal eddy calculation
distForVelocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% STARTING POINT OF THE LOOP %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cpt = 1;
disp('loop about to start')

for x=2:final-1
		tic
	readBudgetFile   %% HERE CELLS ARE CHOSEN
	dNdt
	horizontalFlux
	verticalFlux

if bio==1
if (strcmp(vname,'O2')==1)
cycle_oxygen
elseif (strcmp(vname,'NO3')==1)
cycle_nitrate
elseif (strcmp(vname,'NH4')==1)
cycle_ammonium
end
end % bio

		cpt = cpt+1 ;
		toc
	disp([num2str(x),'/',num2str(final)])
end % x
depthmax
end % dd

end % SCENARIO
