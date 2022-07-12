

    VAR = ncread(bfile,'NITRIF');
    VAR   = permute(VAR, [3 2 1]);
clear VAR1
	for i=1:length(indz) ; VAR1(i) = VAR(indz(i),indy(i),indx(i)) ; i; end
	NITRIF = squeeze(nansum(nansum(nansum(VAR1 .* d2b .* area3d2)))) .*dt ; % mmol/m3/s * m2 * m * s = mmol per time step (dt, usually = 1 month)

    VAR = ncread(bfile,'Denitrif');
    VAR   = permute(VAR, [3 2 1]);
clear VAR1
        for i=1:length(indz) ; VAR1(i) = VAR(indz(i),indy(i),indx(i)) ; i; end
        Denitrif = squeeze(nansum(nansum(nansum(VAR1 .* d2b .* area3d2)))) .*dt ; % mmol/m3/s * m2 * m * s = mmol per time step (dt, usually = 1 month)

    VAR = ncread(bfile,'no3_v_sp');
    VAR   = permute(VAR, [3 2 1]);
clear VAR1
        for i=1:length(indz) ; VAR1(i) = VAR(indz(i),indy(i),indx(i)) ; i; end
        no3_v_sp = squeeze(nansum(nansum(nansum(VAR1 .* d2b .* area3d2)))) .*dt ;

    VAR = ncread(bfile,'no3_v_diat');
    VAR   = permute(VAR, [3 2 1]);
clear VAR1
        for i=1:length(indz) ; VAR1(i) = VAR(indz(i),indy(i),indx(i)) ; i; end
        no3_v_diat = squeeze(nansum(nansum(nansum(VAR1 .* d2b .* area3d2)))) .*dt ;

d2bot = squeeze(dz2b(1,:,:)) ;

    if (min(indz(:)) == 1) % box includes surface layer
        % need to consider air-sea O2 flux
        VAR = ncread(bfile,'Sed_denitr')';
%	d2bot = squeeze(dz2(1,:,:)) ;
%	air_sea_flux = FG_O2(mask_box==1) .* dt.*area(mask_box==1)./d2surf(mask_box==1) ;
	Sed_denitr =nansum( VAR(mask_box==k) .* dt.*area(mask_box==k).* d2bot(mask_box==k) ); 
	else 
	Sed_denitr =0 ;
    end


air_sea_flux = 0 ; % ./nansum(area(mask_box==k)) ; % always applies, could be all zeros

bgc_flux = NITRIF - Denitrif - no3_v_diat - no3_v_sp - Sed_denitr ; % ./nansum(area(mask_box==k)) ;

%% SAVE OUTPUTS
%        ncwrite(fout_air_sea_flux, 'var',air_sea_flux , [cpt]);
        ncwrite(fout, 'nitrification',NITRIF , [cpt]);
        ncwrite(fout, 'denitrification',Sed_denitr+Denitrif , [cpt]);
        ncwrite(fout, 'uptake_diat',no3_v_diat , [cpt]);
        ncwrite(fout, 'uptake_sp',no3_v_sp , [cpt]);
        ncwrite(fout, 'bgc', bgc_flux , [cpt]);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EDDY BGC %%
%% %% %% %% %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Simu==1
sbfile = [rep_clim,'ussw1_bgc_flux_avg.M',num2str(month, '%2.2d'),'_1997_2007.nc'];
elseif Simu==4
sbfile = [rep_clim,'usw42_bgc_flux_avg.M',num2str(month, '%2.2d'),'_1997_2007.nc'];
end
    VAR = ncread(sbfile,'NITRIF');
    VAR   = permute(VAR, [3 2 1]);
clear VAR1
        for i=1:length(indz) ; VAR1(i) = VAR(indz(i),indy(i),indx(i)) ; i; end
        SNITRIF = squeeze(nansum(nansum(nansum(VAR1 .* d2 .* area3d2)))) .*dt ;

    VAR = ncread(sbfile,'Denitrif');
    VAR   = permute(VAR, [3 2 1]);
clear VAR1
        for i=1:length(indz) ; VAR1(i) = VAR(indz(i),indy(i),indx(i)) ; i; end
        SDenitrif = squeeze(nansum(nansum(nansum(VAR1 .* d2 .* area3d2)))) .*dt ;

    VAR = ncread(sbfile,'no3_v_sp');
    VAR   = permute(VAR, [3 2 1]);
clear VAR1
        for i=1:length(indz) ; VAR1(i) = VAR(indz(i),indy(i),indx(i)) ; i; end
        Sno3_v_sp = squeeze(nansum(nansum(nansum(VAR1 .* d2 .* area3d2)))) .*dt ;

    VAR = ncread(sbfile,'no3_v_diat');
    VAR   = permute(VAR, [3 2 1]);
clear VAR1
        for i=1:length(indz) ; VAR1(i) = VAR(indz(i),indy(i),indx(i)) ; i; end
        Sno3_v_diat = squeeze(nansum(nansum(nansum(VAR1 .* d2 .* area3d2)))) .*dt ;

    if (min(indz(:)) == 1) % box includes surface layer
        % need to consider air-sea O2 flux
        VAR = ncread(sbfile,'Sed_denitr')';
        d2bot = squeeze(dz2(1,:,:)) ;
%       air_sea_flux = FG_O2(mask_box==1) .* dt.*area(mask_box==1)./d2surf(mask_box==1) ;
        SSed_denitr =nansum( VAR(mask_box==k) .* dt.*area(mask_box==k).* d2bot(mask_box==k) );
        else
        SSed_denitr =0 ;
    end

Sbgc_flux = SNITRIF - SDenitrif - Sno3_v_diat - Sno3_v_sp - SSed_denitr ; % ./nansum(area(mask_box==k)) ;

ebgc_flux= bgc_flux-Sbgc_flux ;
eNITRIF = NITRIF-SNITRIF ;
eDenitrif = Denitrif-SDenitrif ;
eno3_v_diat = no3_v_diat-Sno3_v_diat ;
eno3_v_sp = no3_v_sp-Sno3_v_sp ;
eSed_denitr = Sed_denitr-SSed_denitr ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% end EDDY BGC PART %% %% %% %% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% SAVE OUTPUTS
%        ncwrite(fout_air_sea_flux, 'var',air_sea_flux , [cpt]);
        ncwrite(fout, 'nitrification',NITRIF , [cpt]);
        ncwrite(fout, 'denitrification',Sed_denitr+Denitrif , [cpt]);
        ncwrite(fout, 'uptake_diat',no3_v_diat , [cpt]);
        ncwrite(fout, 'uptake_sp',no3_v_sp , [cpt]);
        ncwrite(fout, 'bgc', bgc_flux , [cpt]);


