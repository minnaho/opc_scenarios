

    VAR = ncread(bfile,'NITRIF');
    VAR   = permute(VAR, [3 2 1]);
%% remove below depth limit
VAR(zbot2b>depthmax | zbot2b<=depthmin)=NaN ;
clear VAR1
	for i=1:length(indz) ; VAR1(i) = VAR(indz(i),indy(i),indx(i)) ; i; end
	NITRIF = squeeze(nansum(nansum(nansum(VAR1 .* d2b .* area3d2)))) .*dt ; % mmol/m3/s * m2 * m * s = mmol per time step (dt, usually = 1 month)

    VAR = ncread(bfile,'nh4_v_sp');
    VAR   = permute(VAR, [3 2 1]);
%% remove below depth limit
VAR(zbot2b>depthmax | zbot2b<=depthmin)=NaN ;
clear VAR1
        for i=1:length(indz) ; VAR1(i) = VAR(indz(i),indy(i),indx(i)) ; i; end
        nh4_v_sp = squeeze(nansum(nansum(nansum(VAR1 .* d2b .* area3d2)))) .*dt ;

    VAR = ncread(bfile,'nh4_v_diat');
    VAR   = permute(VAR, [3 2 1]);
%% remove below depth limit
VAR(zbot2b>depthmax | zbot2b<=depthmin)=NaN ;
clear VAR1
        for i=1:length(indz) ; VAR1(i) = VAR(indz(i),indy(i),indx(i)) ; i; end
        nh4_v_diat = squeeze(nansum(nansum(nansum(VAR1 .* d2b .* area3d2)))) .*dt ;

QCaCO3 = 0.4 ;
parm_labile_ratio = 0.7 ;
Q = 0.137 ;

    VAR = ncread(bfile,'ZOO_LOSS');
    VAR   = permute(VAR, [3 2 1]);
%% remove below depth limit
VAR(zbot2b>depthmax | zbot2b<=depthmin)=NaN ;
clear VAR1
        for i=1:length(indz) ; VAR1(i) = VAR(indz(i),indy(i),indx(i)) ; i; end
        ZOO_LOSS = squeeze(nansum(nansum(nansum(VAR1 .* d2b .* area3d2)))) .*dt ;

    VAR = ncread(bfile,'SP_LOSS');
    VAR   = permute(VAR, [3 2 1]);
%% remove below depth limit
VAR(zbot2b>depthmax | zbot2b<=depthmin)=NaN ;
clear VAR1
        for i=1:length(indz) ; VAR1(i) = VAR(indz(i),indy(i),indx(i)) ; i; end
        SP_LOSS = squeeze(nansum(nansum(nansum(VAR1 .* d2b .* area3d2)))) .*dt ;

    VAR = ncread(bfile,'GRAZE_SP');
    VAR   = permute(VAR, [3 2 1]);
%% remove below depth limit
VAR(zbot2b>depthmax | zbot2b<=depthmin)=NaN ;
clear VAR1
        for i=1:length(indz) ; VAR1(i) = VAR(indz(i),indy(i),indx(i)) ; i; end
        GRAZE_SP = squeeze(nansum(nansum(nansum(VAR1 .* d2b .* area3d2)))) .*dt ;

    VAR = ncread(bfile,'DIAT_LOSS');
    VAR   = permute(VAR, [3 2 1]);
%% remove below depth limit
VAR(zbot2b>depthmax | zbot2b<=depthmin)=NaN ;
clear VAR1
        for i=1:length(indz) ; VAR1(i) = VAR(indz(i),indy(i),indx(i)) ; i; end
        DIAT_LOSS = squeeze(nansum(nansum(nansum(VAR1 .* d2b .* area3d2)))) .*dt ;

    VAR = ncread(bfile,'GRAZE_DIAT');
    VAR   = permute(VAR, [3 2 1]);
%% remove below depth limit
VAR(zbot2b>depthmax | zbot2b<=depthmin)=NaN ;
clear VAR1
        for i=1:length(indz) ; VAR1(i) = VAR(indz(i),indy(i),indx(i)) ; i; end
        GRAZE_DIAT = squeeze(nansum(nansum(nansum(VAR1 .* d2b .* area3d2)))) .*dt ;

    VAR = ncread(bfile,'DIAZ_LOSS_AGG');
    VAR   = permute(VAR, [3 2 1]);
%% remove below depth limit
VAR(zbot2b>depthmax | zbot2b<=depthmin)=NaN ;
clear VAR1
        for i=1:length(indz) ; VAR1(i) = VAR(indz(i),indy(i),indx(i)) ; i; end
        DIAZ_LOSS = squeeze(nansum(nansum(nansum(VAR1 .* d2b .* area3d2)))) .*dt ;

    VAR = ncread(bfile,'GRAZE_DIAZ');
    VAR   = permute(VAR, [3 2 1]);
%% remove below depth limit
VAR(zbot2b>depthmax | zbot2b<=depthmin)=NaN ;
clear VAR1
        for i=1:length(indz) ; VAR1(i) = VAR(indz(i),indy(i),indx(i)) ; i; end
        GRAZE_DIAZ = squeeze(nansum(nansum(nansum(VAR1 .* d2b .* area3d2)))) .*dt ;

    VAR = ncread(bfile,'POC_REMIN');
    VAR   = permute(VAR, [3 2 1]);
%% remove below depth limit
VAR(zbot2b>depthmax | zbot2b<=depthmin)=NaN ;
clear VAR1
        for i=1:length(indz) ; VAR1(i) = VAR(indz(i),indy(i),indx(i)) ; i; end
        POC_REMIN = squeeze(nansum(nansum(nansum(VAR1 .* d2b .* area3d2)))) .*dt ;

    VAR = ncread(bfile,'DON_REMIN');
    VAR   = permute(VAR, [3 2 1]);
%% remove below depth limit
VAR(zbot2b>depthmax | zbot2b<=depthmin)=NaN ;
clear VAR1
        for i=1:length(indz) ; VAR1(i) = VAR(indz(i),indy(i),indx(i)) ; i; end
        DON_REMIN = squeeze(nansum(nansum(nansum(VAR1 .* d2b .* area3d2)))) .*dt ;


    VAR = ncread(bfile,'Sed_Flux_POC')';
	remin_sed_poc = nan(NZ,NY,NX);
	remin_sed_poc(1,:,:) = VAR;
%% remove below depth limit
remin_sed_poc(zbot2b>depthmax | zbot2b<=depthmin)=NaN ;
clear VAR1
        for i=1:length(indz) ; VAR1(i) = remin_sed_poc(indz(i),indy(i),indx(i)) ; i; end
        REMIN_SED_POC = squeeze(nansum(nansum(nansum(VAR1 .* d2b .* area3d2)))) .*dt ;

sp_loss_poc = QCaCO3 * SP_LOSS ;
epsC      = 1.00e-8 ; % small C concentration (mmol C/m^3)
epsTinv   = 3.17e-8 ; % small inverse time scale (1/year) (1/sec)
f_zoo_detr = (0.1333 * (GRAZE_DIAT + epsC * epsTinv) + ...
         0.0333 * (GRAZE_SP + epsC * epsTinv) + ...
         0.0 * (GRAZE_DIAZ + epsC * epsTinv)) / ...
         (GRAZE_DIAT + GRAZE_SP + GRAZE_DIAZ + 3.0 * epsC * epsTinv) ;

zoo_loss_dic = (1-parm_labile_ratio).*(1-f_zoo_detr).*ZOO_LOSS ;
sp_loss_dic = parm_labile_ratio .* (SP_LOSS - sp_loss_poc) ;
graze_sp_dic = 0.36 .* GRAZE_SP ;
diat_loss_dic = parm_labile_ratio .* 0.95 .* DIAT_LOSS ;
graze_diat_dic = 0.31 .* GRAZE_DIAT ;

diaz_loss_dic = parm_labile_ratio * DIAZ_LOSS ;
graze_diaz_dic = 0.55 * GRAZE_DIAZ ;
parm_kappa_nitrif   = 0.06 ;
%ammox = parm_kappa_nitrif * NH4 ;

bgc_flux = -(nh4_v_diat+nh4_v_sp) + ((zoo_loss_dic + sp_loss_dic + graze_sp_dic + diat_loss_dic + graze_diat_dic + POC_REMIN + diaz_loss_dic + graze_diaz_dic).*Q) + DON_REMIN - NITRIF + (Q*REMIN_SED_POC) ;

air_sea_flux = 0 ; % ./nansum(area(mask_box==k)) ; % always applies, could be all zeros

%bgc_flux = -(hor_flux+vert_flux) + DnDt_int_box ;
amm_prod = bgc_flux - (NITRIF+nh4_v_diat+nh4_v_sp);

%% SAVE OUTPUTS
%        ncwrite(fout_air_sea_flux, 'var',air_sea_flux , [cpt]);
        ncwrite(fout, 'nitrification',NITRIF , [cpt]);
        ncwrite(fout, 'amm_prod',amm_prod , [cpt]);
        ncwrite(fout, 'uptake_diat',nh4_v_diat , [cpt]);
        ncwrite(fout, 'uptake_sp',nh4_v_sp , [cpt]);
        ncwrite(fout, 'bgc', bgc_flux , [cpt]);

return


