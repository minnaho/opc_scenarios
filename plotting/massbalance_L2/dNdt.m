
clear c1 c2 cm
for i=1:length(indz) ; c1(i) = conc1(indz(i),indy(i),indx(i)) ; i; end
for i=1:length(indz) ; c2(i) = conc2(indz(i),indy(i),indx(i)) ; i; end
for i=1:length(indz) ; cm(i) = concm(indz(i),indy(i),indx(i)) ; i; end

m1 = squeeze(nansum(nansum(nansum(c1 .* d1 .* area3d2)))) ; % full column, result is in mmol/m2
m2 = squeeze(nansum(nansum(nansum(c2 .* d2 .* area3d2)))) ;
mm = squeeze(nansum(nansum(nansum(cm .* d2b .* area3d2)))) ;

DnDt_int_box = (m2 - m1) ; %./nansum(area(mask_box==k)) ; % mmol/dt
printvolume = nansum(d2 .* area3d2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE OUTPUTS

        ncwrite(fout, 'dNdT', DnDt_int_box, [cpt]);
        ncwrite(fout, 'inv1',m1 , [cpt]);
        ncwrite(fout, 'inv2',m2 , [cpt]);
        ncwrite(fout, 'invm',mm , [cpt]);

