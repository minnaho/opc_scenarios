
    VAR = ncread(bfile,'J_O2');
    VAR   = permute(VAR, [3 2 1]);
%% remove below depth limit
VAR(zbot2b>depthmax | zbot2b<=depthmin)=NaN ;
clear VAR1
        for i=1:length(indz) ; VAR1(i) = VAR(indz(i),indy(i),indx(i)) ; i; end
        OXYGB = squeeze(nansum(nansum(nansum(VAR1 .* d2b .* area3d2)))) .*dt ; % mmol/m3/s * m2 * m * s = mmol per time step (dt, usually = 1 month)

    if (max(indz(:)) == NZ) % box includes surface layer
d2surf = squeeze(dz2b(60,:,:)) ;

        % need to consider air-sea O2 flux
        VAR = ncread(bfile,'FG_O2')';
%	air_sea_flux = FG_O2(mask_box==1) .* dt.*area(mask_box==1).*d2surf(mask_box==1) ;
	air_sea_flux =nansum( VAR(mask_box==k) .* dt.*area(mask_box==k)); 
        else air_sea_flux =0 ;
    end
bgc_flux = OXYGB +air_sea_flux; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ncwrite(fout, 'J_O2',OXYGB , [cpt]);
        ncwrite(fout, 'air_sea_flux',air_sea_flux , [cpt]);
        ncwrite(fout, 'bgc', bgc_flux , [cpt]);


