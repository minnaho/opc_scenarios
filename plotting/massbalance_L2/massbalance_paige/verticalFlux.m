    % hf: vertical fluxes are in mmol m-2 s-1

    name_va = sprintf('VertAdvFlux_%s',vname) ; 
    va = ncread(pfile, name_va) .* dt; %mmol/m2/dt
    name_vd = sprintf('VertDiffFlux_%s',vname) ; 
    vd = ncread(pfile, name_vd) .* dt; %mmol/m2/dt

        va = permute(va, [3 2 1]);
        vd = permute(vd, [3 2 1]);

%% remove below depth limit
%va(zbot2b>depthmax | zbot2b<=depthmin)=NaN ;
%vd(zbot2b>depthmax | zbot2b<=depthmin)=NaN ;

clear vacbot vdcbot vactop vdctop
for i=1:length(indz) ; vacbot(i) = va(indz(i),indy(i),indx(i)) ; i; end
for i=1:length(indz) ; vdcbot(i) = vd(indz(i),indy(i),indx(i)) ; i; end
for i=1:length(indz) ; vactop(i) = va(indz(i)+1,indy(i),indx(i)) ; i; end
for i=1:length(indz) ; vdctop(i) = vd(indz(i)+1,indy(i),indx(i)) ; i; end

vacb = nansum(vacbot(:).*area3d2(:)) ;
vact = nansum(vactop(:).*area3d2(:)) ;
vdcb = nansum(vdcbot(:).*area3d2(:)) ;
vdct = nansum(vdctop(:).*area3d2(:)) ;
vac = nansum([-vact+vacb]) ;
vdc = nansum([-vdct+vdcb]) ;

vac = nansum( [vacbot.*area3d2 -vactop.*area3d2] ) ;
vdc = nansum( [vdcbot.*area3d2 -vdctop.*area3d2] ) ;
    vert_flux = (vac + vdc) ; % ./nansum(area(mask_box==k))  ; % mmol/dt
    phys_flux = hor_flux + vert_flux ; 

printvolume = nansum(d2 .* area3d2);
printarea = nansum(area(mask_box==k).* d2bot(mask_box==k)) ;
%% SAVE OUTPUTS
% physics
%        ncwrite(fout, 'vert_flux_eddy', weddy , [cpt]);
        ncwrite(fout, 'adv_vert_flux', vac , [cpt]);
        ncwrite(fout, 'diffusive_flux', vdc , [cpt]);
        ncwrite(fout, 'vertical_flux', vert_flux , [cpt]);

