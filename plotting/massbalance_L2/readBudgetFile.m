
avgfile1 =  [rep_mod,'/',repavg(x,1).name] ;
pfile  =  [rep_mod,'/',repp(x,1).name] ;
bfile  =  [rep_mod,'/',repb(x,1).name] ;

clear RSTY2 RSTM2 RSTY1 RSTM1

testname = repavg(x,1).name ;

if Simu==2
RSTY2 = testname(13:16) ;
RSTM2 = testname(18:19) ;
else
RSTY2 = testname(12:15) ;
RSTM2 = testname(17:18) ;
end

if str2num(RSTM2)==1
RSTY1 = str2num(RSTY2)-1;
RSTM1 = 12;
else
RSTY1 = str2num(RSTY2);
RSTM1 = str2num(RSTM2)-1;
end

rstrep1 = dir([rep_his,'/l2_scb_his.Y',num2str(RSTY1),'M',num2str(RSTM1, '%2.2d'),'*.nc']) ;
rstrep2 = dir([rep_his,'/l2_scb_his.Y',num2str(RSTY2),'M',num2str(RSTM2, '%2.2d'),'*.nc']) ;

hfile1 =  [rep_his,'/',rstrep1(1,1).name] ;
hfile2 =  [rep_his,'/',rstrep2(1,1).name] ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfile1 = hfile1 ;
cfile2 = hfile2 ;

%end
%%
          disp(['Reading ', cfile2])

%    theta_s = 6.0;
%    theta_b = 3.0;
%    hc = 250;
%    sc_type = 'new2012'; % for my zlevs4!!
%mask_rho= ncread(grd, 'mask_rho')';
%h   = ncread(grd, 'h')';

% hf: these should be in the grid file, but may not always be
if (exist('mask_u','var') ~= 1)
    mask_u = rho2u(mask_rho);
end
if (exist('mask_v','var') ~= 1)
    mask_v = rho2v(mask_rho);
end

%conc1 = ncread(cfile1,vname) ;
%conc2 = ncread(cfile2,vname) ;

if (strcmp(vname,'PON')==1)
conc1a = ncread(cfile1,'ZOOC').*c2n ;
conc1b = ncread(cfile1,'DIATC').*c2n ;
conc1c = ncread(cfile1,'DIAZC').*c2n ;
conc1d = ncread(cfile1,'SPC').*c2n ;
conc1e = ncread(cfile1,'DOC').*c2n ;
conc1 =  conc1a+conc1b+conc1c+conc1d+conc1e  ;
conc2a = ncread(cfile2,'ZOOC').*c2n ;
conc2b = ncread(cfile2,'DIATC').*c2n ;
conc2c = ncread(cfile2,'DIAZC').*c2n ;
conc2d = ncread(cfile2,'SPC').*c2n ;
conc2e = ncread(cfile2,'DOC').*c2n ;
conc2 =  conc2a+conc2b+conc2c+conc2d+conc2e  ;
concma = ncread(avgfile1,'ZOOC').*c2n ;
concmb = ncread(avgfile1,'DIATC').*c2n ;
concmc = ncread(avgfile1,'DIAZC').*c2n ;
concmd = ncread(avgfile1,'SPC').*c2n ;
concme = ncread(avgfile1,'DOC').*c2n ;
concm =  concma+concmb+concmc+concmd+concme  ;
else
conc1 = ncread(cfile1,vname) ;
conc2 = ncread(cfile2,vname) ;
concm = ncread(avgfile1,vname) ;
end

[ac bc cc dc] = size(conc1) ;
if dc>1
conc1 = squeeze(conc1(:,:,:,2)) ;
end
[ac bc cc dc] = size(conc2) ;
if dc>1
conc2 = squeeze(conc2(:,:,:,2)) ;
end

%conc1 = squeeze(conc1(:,:,:,2)) ;
%conc2 = squeeze(conc2(:,:,:,2)) ;

conc1 = permute(conc1, [3 2 1]);
conc2 = permute(conc2, [3 2 1]);
concm = permute(concm, [3 2 1]);

zeta1 = ncread(cfile1,'zeta') ;
zeta2 = ncread(cfile2,'zeta') ;
zetab2 = ncread(avgfile1,'zeta') ;
%zeta1 = squeeze(zeta1(:,:,2))' ;
%zeta2 = squeeze(zeta2(:,:,2))' ;

[zetaa zetab zetac] = size(zeta1) ;
if zetac>1
zeta1 = squeeze(zeta1(:,:,2))' ;
else
zeta1=zeta1';
end
[zetaa zetab zetac] = size(zeta2) ;
if zetac>1
zeta2 = squeeze(zeta2(:,:,2))' ;
else
zeta2=zeta2';
end

[z_w1,Cw1] = zlevs4(h, zeta1, theta_s, theta_b, hc, NZ, 'w',sc_type);
[z_w2,Cw2] = zlevs4(h, zeta2, theta_s, theta_b, hc, NZ, 'w',sc_type);
[z_w2b,Cw2b] = zlevs4(h, zetab2', theta_s, theta_b, hc, NZ, 'w',sc_type);
dz1 = diff(z_w1);
dz2 = diff(z_w2);
dz2b = diff(z_w2b);

        zbot1 = flipdim(cumsum(flipdim(dz1,1)),1);
        ztop1 = [zbot1(2:end,:,:);zeros(1,NY,NX)];
        zbot2 = flipdim(cumsum(flipdim(dz2,1)),1);
        ztop2 = [zbot2(2:end,:,:);zeros(1,NY,NX)];
        zbot2b = flipdim(cumsum(flipdim(dz2b,1)),1);
        ztop2b = [zbot2b(2:end,:,:);zeros(1,NY,NX)];
	z2b = (ztop2b + zbot2b)./2 ;

%mask_box3drup = mask_box3dr ;
%mask_box3drdn = mask_box3dr ;
%mask_box3dr(zbot2b>=depthmax)=NaN ;
mask_box3dr(zbot2b>depthmax | ztop2b<=depthmin)=NaN ;
[indz indy indx] = ind2sub(size(mask_box3dr),find(mask_box3dr == k));

%mask_box3drdn(zbot2b<50 | zbot2b>60)=NaN ;
%mask_box3drup(zbot2b<40 | zbot2b>50)=NaN ;
%[indz60 indy60 indx60] = ind2sub(size(mask_box3drdn),find(mask_box3drdn == k));
%[indz40 indy40 indx40] = ind2sub(size(mask_box3drup),find(mask_box3drup == k));

%[indx,indy] = FindCloestPoint_ROMS( lon, lat, -118.62, 33.8, mask )
%indz = 60 ;

%for i=1:length(indz) ; list(i) = zbot2b(indz(i),indy(i),indx(i)) ; end
%for i=1:length(indz40) ; list40(i) = zbot2b(indz40(i),indy40(i),indx40(i)) ; end
%for i=1:length(indz60) ; list60(i) = zbot2b(indz60(i),indy60(i),indx60(i)) ; end

[indy2 indx2] = find(mask_box==1) ;
listy2 = unique(indy2) ;
listx2 = unique(indx2) ;



%[indzr indyr indxr] = ind2sub(size(mask_box3d_red),find(mask_box3d_red == k));

%list = find(indz1==60 & indy1==min(indy1)) ; 
%
%indz = indz1(list(end)) % max(indz1) max(indz1)] ;
%indx = indx1(list(end))
%indx = max(indx1)-5 : max(indx1)-3 ;
%indx = max(indx1)-3 : max(indx1)-1 ;
%indy = indy1(list(end))% min(indy1) min(indy1)];

%c1 = conc1(indz,indy,indx) ;
%c2 = conc2(indz,indy,indx) ;
% missing values are 1e+33 - replace with nan so that nansum etc. can be used later
% WARN make sure that's the case for the run under investigation!

%% TIME
t1 = ncread(cfile1,tname) ;
t2 = ncread(cfile2,tname) ;
sizet1=size(t1);
sizet2=size(t2);
if sizet1(1)>1 ; t1 = t1(2) ;end ;
if sizet2(1)>1 ; t2 = t2(2) ;end ;
%t1 = t1(2) ;
%t2 = t2(2) ;
dt = t2 -t1 ; % in seconds
month=RSTM2;
        ncwrite(fout, 'time',(t1/86400+datenum(1994,1,1)) , [cpt]);
        ncwrite(fout, 'dt',dt , [cpt]);
%d1 = dz1(indz,indy,indx) ;
%d2 = dz2(indz,indy,indx) ;

%% DEPTH
clear d1 d2 d2b
for i=1:length(indz) ; d1(i) = dz1(indz(i),indy(i),indx(i)) ; i; end
for i=1:length(indz) ; d2(i) = dz2(indz(i),indy(i),indx(i)) ; i; end
for i=1:length(indz) ; d2b(i) = dz2b(indz(i),indy(i),indx(i)) ; i; end
d2bot = squeeze(dz2b(1,:,:)) ;

%% VOLUME
%area = 1./(pm(indy,indx).*pn(indy,indx)); % m2
area = 1./(pm.*pn); % m2
area3d = repmat(area,1,1,NZ);
area3d = permute(area3d, [3 1 2]);
area3dv = repmat(area,1,1,NZ+1);
area3dv = permute(area3dv, [3 1 2]);
clear area3d2
for i=1:length(indz) ; area3d2(i) = area3d(indz(i),indy(i),indx(i)) ; i; end
vol2 =  nansum(d2b .* area3d2) ;
printvolume = nansum(d2b .* area3d2);
printarea = nansum(area(mask_box==k).* d2bot(mask_box==k)) ;

%% ZONAL AREA
pm3d = repmat(1./pm,1,1,NZ); % m
pm3d = permute(pm3d, [3 1 2]);
clear pm3d2
for i=1:length(indz) ; pm3d2(i) = pm3d(indz(i),indy(i),indx(i)) ; i; end
surf_x = pm3d2.*d2;
surfx = nansum(surf_x);
surf_dx = dz2.*pm3d ;

%% MERIDIONAL AREA
pn3d = repmat(1./pn,1,1,NZ); % m
pn3d = permute(pn3d, [3 1 2]);
clear pn3d2
for i=1:length(indz) ; pn3d2(i) = pn3d(indz(i),indy(i),indx(i)) ; i; end
surf_y = pn3d2.*d2b;
surfy = nansum(surf_y);
surf_dy = dz2b.*pn3d ;

%% SAVE DIMENSIONS
        ncwrite(fout, 'volume',printvolume , [cpt]);
        ncwrite(fout, 'area',printarea , [cpt]);


