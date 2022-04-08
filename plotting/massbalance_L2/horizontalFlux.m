    % hf: horizontal fluxes are in mmol/s
    name_hx = sprintf('HorXAdvFlux_%s',vname);
    hx=ncread(pfile, name_hx) .* dt; % mmol/dt
    name_hy = sprintf('HorYAdvFlux_%s',vname);
    hy=ncread(pfile, name_hy) .* dt; % mmol/dt
        hx = permute(hx, [3 2 1]);
        hy = permute(hy, [3 2 1]);

%% hf
clear Thx_w Thx_e Thy_s Thy_n
for i=1:length(indz) ; Thx_w(i) = hx(indz(i),indy(i),indx(i)) ; i; end
for i=1:length(indz) ; Thx_e(i) = hx(indz(i),indy(i),indx(i)+1) ; i; end
for i=1:length(indz) ; Thy_s(i) = hy(indz(i),indy(i),indx(i)) ; i; end
for i=1:length(indz) ; Thy_n(i) = hy(indz(i),indy(i)+1,indx(i)) ; i; end

hx_w = nansum(Thx_w(:)) ;
hx_e = nansum(Thx_e(:)) ;
hy_s = nansum(Thy_s(:)) ;
hy_n = nansum(Thy_n(:)) ;

zonal = nansum([-hx_e+hx_w]) ;
merid = nansum([-hy_n+hy_s]) ;

    hor_flux = zonal+merid ;
 

%along_cross

%% SAVE OUTPUTS
% physics
        ncwrite(fout, 'horizontal_flux', hor_flux , [cpt]);

        ncwrite(fout, 'zonal', zonal , [cpt]);
        ncwrite(fout, 'meridional', merid , [cpt]);

        ncwrite(fout, 'east', hx_e , [cpt]);
        ncwrite(fout, 'west', hx_w , [cpt]);
        ncwrite(fout, 'north', hy_n , [cpt]);
        ncwrite(fout, 'south', hy_s , [cpt]);


