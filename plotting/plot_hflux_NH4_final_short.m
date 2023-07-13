addpath(genpath('/data/project3/kesf/tools_matlab/matlab_paths/'))
call_ncviewcolors

[lon,lat,NX,NY,h,area,mask_rho]=loadgrid(2) ;

dir_gr = '/data/project3/minnaho/opc_scenarios/plotting/figs/2years/';

% PNDN only realistic
load('/data/project3/kesf/sharing/budget_ww/CoupeFlux_Shelf200m_NH4_PNDN_only_realistic.mat')

dist = CoupeFlux.dist ; 
distc = CoupeFlux.distc ;
area = CoupeFlux.Area ;
depth = CoupeFlux.dpth ;
onln = 86400.*CoupeFlux.FluxDIN_onln./area ;
nmean = 86400.*CoupeFlux.FluxDIN_mean./area ;
eddy = onln-nmean;
date= CoupeFlux.date ;

clear mean_int eddy_int tot_int
[a b c] = size(nmean);
mean_int=zeros(a,b,999) ;
eddy_int=zeros(a,b,999) ;
tot_int=zeros(a,b,999) ;
for i=1:b
for t=1:a
aa1=squeeze(nmean(t,i,:)) ;
aa2=squeeze(eddy(t,i,:)) ;
aa3=squeeze(onln(t,i,:)) ;
bb=squeeze(depth(t,i,:)) ;

if isnan(bb(1))
mean_int(t,i,:) = [-1000:1:-2].*0 + NaN;
eddy_int(t,i,:) = [-1000:1:-2].*0 + NaN;
tot_int(t,i,:) = [-1000:1:-2].*0 + NaN;
else
mean_int(t,i,:) = interp1(bb,aa1,[-1000:1:-2]);
eddy_int(t,i,:) = interp1(bb,aa2,[-1000:1:-2]);
tot_int(t,i,:) = interp1(bb,aa3,[-1000:1:-2]);
end

end
end

mean_int(mean_int==0)=NaN;
eddy_int(eddy_int==0)=NaN;
tot_int(tot_int==0)=NaN;

tmean1 = squeeze(nanmean(mean_int,1));
tzmean1 = squeeze(nanmean(mean_int,[1 3]));
tpmean1 = squeeze(nanmean(mean_int,[1 2]));
zpmean1 = squeeze(nanmean(mean_int,[2 3]));

teddy1 = squeeze(nanmean(eddy_int,1));
tzeddy1 = squeeze(nanmean(eddy_int,[1 3]));
tpeddy1 = squeeze(nanmean(eddy_int,[1 2]));
zpeddy1 = squeeze(nanmean(eddy_int,[2 3]));

ttot1 = squeeze(nanmean(tot_int,1));
tztot1 = squeeze(nanmean(tot_int,[1 3]));
tptot1 = squeeze(nanmean(tot_int,[1 2]));
zptot1 = squeeze(nanmean(tot_int,[2 3]));

mmonth = str2num(datestr(date,'mm'));

list_win = find(mmonth>=9 | mmonth<=3);
list_sum = find(mmonth>=4 & mmonth<=8);
tptotwin1 = squeeze(nanmean(tot_int(list_win,:,:),[1 2]));
tptotsum1 = squeeze(nanmean(tot_int(list_sum,:,:),[1 2]));
tpmeanwin1 = squeeze(nanmean(mean_int(list_win,:,:),[1 2]));
tpmeansum1 = squeeze(nanmean(mean_int(list_sum,:,:),[1 2]));
tpeddywin1 = squeeze(nanmean(eddy_int(list_win,:,:),[1 2]));
tpeddysum1 = squeeze(nanmean(eddy_int(list_sum,:,:),[1 2]));

% PNDN 50 fixriver 
load('/data/project3/kesf/sharing/budget_ww/CoupeFlux_Shelf200m_NH4_pndn50_fixriver.mat')

dist = CoupeFlux.dist ; 
distc = CoupeFlux.distc ;
area = CoupeFlux.Area ;
depth = CoupeFlux.dpth ;
onln = 86400.*CoupeFlux.FluxDIN_onln./area ;
nmean = 86400.*CoupeFlux.FluxDIN_mean./area ;
eddy = onln-nmean;
date= CoupeFlux.date ;

clear mean_int eddy_int tot_int
[a b c] = size(nmean);
mean_int=zeros(a,b,999) ;
eddy_int=zeros(a,b,999) ;
tot_int=zeros(a,b,999) ;
for i=1:b
for t=1:a
aa1=squeeze(nmean(t,i,:)) ;
aa2=squeeze(eddy(t,i,:)) ;
aa3=squeeze(onln(t,i,:)) ;
bb=squeeze(depth(t,i,:)) ;

if isnan(bb(1))
mean_int(t,i,:) = [-1000:1:-2].*0 + NaN;
eddy_int(t,i,:) = [-1000:1:-2].*0 + NaN;
tot_int(t,i,:) = [-1000:1:-2].*0 + NaN;
else
mean_int(t,i,:) = interp1(bb,aa1,[-1000:1:-2]);
eddy_int(t,i,:) = interp1(bb,aa2,[-1000:1:-2]);
tot_int(t,i,:) = interp1(bb,aa3,[-1000:1:-2]);
end

end
end

mean_int(mean_int==0)=NaN;
eddy_int(eddy_int==0)=NaN;
tot_int(tot_int==0)=NaN;

tmean2 = squeeze(nanmean(mean_int,1));
tzmean2 = squeeze(nanmean(mean_int,[1 3]));
tpmean2 = squeeze(nanmean(mean_int,[1 2]));
zpmean2 = squeeze(nanmean(mean_int,[2 3]));

teddy2 = squeeze(nanmean(eddy_int,1));
tzeddy2 = squeeze(nanmean(eddy_int,[1 3]));
tpeddy2 = squeeze(nanmean(eddy_int,[1 2]));
zpeddy2 = squeeze(nanmean(eddy_int,[2 3]));

ttot2 = squeeze(nanmean(tot_int,1));
tztot2 = squeeze(nanmean(tot_int,[1 3]));
tptot2 = squeeze(nanmean(tot_int,[1 2]));
zptot2 = squeeze(nanmean(tot_int,[2 3]));

mmonth = str2num(datestr(date,'mm'));

list_win = find(mmonth>=9 | mmonth<=3);
list_sum = find(mmonth>=4 & mmonth<=8);
tptotwin2 = squeeze(nanmean(tot_int(list_win,:,:),[1 2]));
tptotsum2 = squeeze(nanmean(tot_int(list_sum,:,:),[1 2]));
tpmeanwin2 = squeeze(nanmean(mean_int(list_win,:,:),[1 2]));
tpmeansum2 = squeeze(nanmean(mean_int(list_sum,:,:),[1 2]));
tpeddywin2 = squeeze(nanmean(eddy_int(list_win,:,:),[1 2]));
tpeddysum2 = squeeze(nanmean(eddy_int(list_sum,:,:),[1 2]));

% PNDN 90 fixriver 
load('/data/project3/kesf/sharing/budget_ww/CoupeFlux_Shelf200m_NH4_pndn90_fixriver.mat')

dist = CoupeFlux.dist ; 
distc = CoupeFlux.distc ;
area = CoupeFlux.Area ;
depth = CoupeFlux.dpth ;
onln = 86400.*CoupeFlux.FluxDIN_onln./area ;
nmean = 86400.*CoupeFlux.FluxDIN_mean./area ;
eddy = onln-nmean;
date= CoupeFlux.date ;

clear mean_int eddy_int tot_int
[a b c] = size(nmean);
mean_int=zeros(a,b,999) ;
eddy_int=zeros(a,b,999) ;
tot_int=zeros(a,b,999) ;
for i=1:b
for t=1:a
aa1=squeeze(nmean(t,i,:)) ;
aa2=squeeze(eddy(t,i,:)) ;
aa3=squeeze(onln(t,i,:)) ;
bb=squeeze(depth(t,i,:)) ;

if isnan(bb(1))
mean_int(t,i,:) = [-1000:1:-2].*0 + NaN;
eddy_int(t,i,:) = [-1000:1:-2].*0 + NaN;
tot_int(t,i,:) = [-1000:1:-2].*0 + NaN;
else
mean_int(t,i,:) = interp1(bb,aa1,[-1000:1:-2]);
eddy_int(t,i,:) = interp1(bb,aa2,[-1000:1:-2]);
tot_int(t,i,:) = interp1(bb,aa3,[-1000:1:-2]);
end

end
end

mean_int(mean_int==0)=NaN;
eddy_int(eddy_int==0)=NaN;
tot_int(tot_int==0)=NaN;

tmean3 = squeeze(nanmean(mean_int,1));
tzmean3 = squeeze(nanmean(mean_int,[1 3]));
tpmean3 = squeeze(nanmean(mean_int,[1 2]));
zpmean3 = squeeze(nanmean(mean_int,[2 3]));

teddy3 = squeeze(nanmean(eddy_int,1));
tzeddy3 = squeeze(nanmean(eddy_int,[1 3]));
tpeddy3 = squeeze(nanmean(eddy_int,[1 2]));
zpeddy3 = squeeze(nanmean(eddy_int,[2 3]));

ttot3 = squeeze(nanmean(tot_int,1));
tztot3 = squeeze(nanmean(tot_int,[1 3]));
tptot3 = squeeze(nanmean(tot_int,[1 2]));
zptot3 = squeeze(nanmean(tot_int,[2 3]));

mmonth = str2num(datestr(date,'mm'));

list_win = find(mmonth>=9 | mmonth<=3);
list_sum = find(mmonth>=4 & mmonth<=8);
tptotwin3 = squeeze(nanmean(tot_int(list_win,:,:),[1 2]));
tptotsum3 = squeeze(nanmean(tot_int(list_sum,:,:),[1 2]));
tpmeanwin3 = squeeze(nanmean(mean_int(list_win,:,:),[1 2]));
tpmeansum3 = squeeze(nanmean(mean_int(list_sum,:,:),[1 2]));
tpeddywin3 = squeeze(nanmean(eddy_int(list_win,:,:),[1 2]));
tpeddysum3 = squeeze(nanmean(eddy_int(list_sum,:,:),[1 2]));

% FNDN only realistic
load('/data/project3/kesf/sharing/budget_ww/CoupeFlux_Shelf200m_NH4_FNDN_only_realistic.mat')

dist = CoupeFlux.dist ;
distc = CoupeFlux.distc ;
area = CoupeFlux.Area ;
depth = CoupeFlux.dpth ;
onln = 86400.*CoupeFlux.FluxDIN_onln./area ;
nmean = 86400.*CoupeFlux.FluxDIN_mean./area ;
eddy = onln-nmean;
date= CoupeFlux.date ;

clear mean_int eddy_int tot_int
[a b c] = size(nmean);
mean_int=zeros(a,b,999) ;
eddy_int=zeros(a,b,999) ;
tot_int=zeros(a,b,999) ;
for i=1:b
for t=1:a
aa1=squeeze(nmean(t,i,:)) ;
aa2=squeeze(eddy(t,i,:)) ;
aa3=squeeze(onln(t,i,:)) ;
bb=squeeze(depth(t,i,:)) ;

if isnan(bb(1))
mean_int(t,i,:) = [-1000:1:-2].*0 + NaN;
eddy_int(t,i,:) = [-1000:1:-2].*0 + NaN;
tot_int(t,i,:) = [-1000:1:-2].*0 + NaN;
else
mean_int(t,i,:) = interp1(bb,aa1,[-1000:1:-2]);
eddy_int(t,i,:) = interp1(bb,aa2,[-1000:1:-2]);
tot_int(t,i,:) = interp1(bb,aa3,[-1000:1:-2]);
end

end
end

mean_int(mean_int==0)=NaN;
eddy_int(eddy_int==0)=NaN;
tot_int(tot_int==0)=NaN;

tmean4 = squeeze(nanmean(mean_int,1));
tzmean4 = squeeze(nanmean(mean_int,[1 3]));
tpmean4 = squeeze(nanmean(mean_int,[1 2]));
zpmean4 = squeeze(nanmean(mean_int,[2 3]));

teddy4 = squeeze(nanmean(eddy_int,1));
tzeddy4 = squeeze(nanmean(eddy_int,[1 3]));
tpeddy4 = squeeze(nanmean(eddy_int,[1 2]));
zpeddy4 = squeeze(nanmean(eddy_int,[2 3]));

ttot4 = squeeze(nanmean(tot_int,1));
tztot4 = squeeze(nanmean(tot_int,[1 3]));
tptot4 = squeeze(nanmean(tot_int,[1 2]));
zptot4 = squeeze(nanmean(tot_int,[2 3]));

mmonth = str2num(datestr(date,'mm'));

list_win = find(mmonth>=9 | mmonth<=3);
list_sum = find(mmonth>=4 & mmonth<=8);
tptotwin4 = squeeze(nanmean(tot_int(list_win,:,:),[1 2]));
tptotsum4 = squeeze(nanmean(tot_int(list_sum,:,:),[1 2]));
tpmeanwin4 = squeeze(nanmean(mean_int(list_win,:,:),[1 2]));
tpmeansum4 = squeeze(nanmean(mean_int(list_sum,:,:),[1 2]));
tpeddywin4 = squeeze(nanmean(eddy_int(list_win,:,:),[1 2]));
tpeddysum4 = squeeze(nanmean(eddy_int(list_sum,:,:),[1 2]));

% FNDN 50 fixriver
load('/data/project3/kesf/sharing/budget_ww/CoupeFlux_Shelf200m_NH4_fndn50_fixriver.mat')

dist = CoupeFlux.dist ;
distc = CoupeFlux.distc ;
area = CoupeFlux.Area ;
depth = CoupeFlux.dpth ;
onln = 86400.*CoupeFlux.FluxDIN_onln./area ;
nmean = 86400.*CoupeFlux.FluxDIN_mean./area ;
eddy = onln-nmean;
date= CoupeFlux.date ;

clear mean_int eddy_int tot_int
[a b c] = size(nmean);
mean_int=zeros(a,b,999) ;
eddy_int=zeros(a,b,999) ;
tot_int=zeros(a,b,999) ;
for i=1:b
for t=1:a
aa1=squeeze(nmean(t,i,:)) ;
aa2=squeeze(eddy(t,i,:)) ;
aa3=squeeze(onln(t,i,:)) ;
bb=squeeze(depth(t,i,:)) ;

if isnan(bb(1))
mean_int(t,i,:) = [-1000:1:-2].*0 + NaN;
eddy_int(t,i,:) = [-1000:1:-2].*0 + NaN;
tot_int(t,i,:) = [-1000:1:-2].*0 + NaN;
else
mean_int(t,i,:) = interp1(bb,aa1,[-1000:1:-2]);
eddy_int(t,i,:) = interp1(bb,aa2,[-1000:1:-2]);
tot_int(t,i,:) = interp1(bb,aa3,[-1000:1:-2]);
end

end
end

mean_int(mean_int==0)=NaN;
eddy_int(eddy_int==0)=NaN;
tot_int(tot_int==0)=NaN;

tmean5 = squeeze(nanmean(mean_int,1));
tzmean5 = squeeze(nanmean(mean_int,[1 3]));
tpmean5 = squeeze(nanmean(mean_int,[1 2]));
zpmean5 = squeeze(nanmean(mean_int,[2 3]));

teddy5 = squeeze(nanmean(eddy_int,1));
tzeddy5 = squeeze(nanmean(eddy_int,[1 3]));
tpeddy5 = squeeze(nanmean(eddy_int,[1 2]));
zpeddy5 = squeeze(nanmean(eddy_int,[2 3]));

ttot5 = squeeze(nanmean(tot_int,1));
tztot5 = squeeze(nanmean(tot_int,[1 3]));
tptot5 = squeeze(nanmean(tot_int,[1 2]));
zptot5 = squeeze(nanmean(tot_int,[2 3]));

mmonth = str2num(datestr(date,'mm'));

list_win = find(mmonth>=9 | mmonth<=3);
list_sum = find(mmonth>=4 & mmonth<=8);
tptotwin5 = squeeze(nanmean(tot_int(list_win,:,:),[1 2]));
tptotsum5 = squeeze(nanmean(tot_int(list_sum,:,:),[1 2]));
tpmeanwin5 = squeeze(nanmean(mean_int(list_win,:,:),[1 2]));
tpmeansum5 = squeeze(nanmean(mean_int(list_sum,:,:),[1 2]));
tpeddywin5 = squeeze(nanmean(eddy_int(list_win,:,:),[1 2]));
tpeddysum5 = squeeze(nanmean(eddy_int(list_sum,:,:),[1 2]));

% FNDN 90 fixriver
load('/data/project3/kesf/sharing/budget_ww/CoupeFlux_Shelf200m_NH4_fndn90_fixriver.mat')

dist = CoupeFlux.dist ;
distc = CoupeFlux.distc ;
area = CoupeFlux.Area ;
depth = CoupeFlux.dpth ;
onln = 86400.*CoupeFlux.FluxDIN_onln./area ;
nmean = 86400.*CoupeFlux.FluxDIN_mean./area ;
eddy = onln-nmean;
date= CoupeFlux.date ;

clear mean_int eddy_int tot_int
[a b c] = size(nmean);
mean_int=zeros(a,b,999) ;
eddy_int=zeros(a,b,999) ;
tot_int=zeros(a,b,999) ;
for i=1:b
for t=1:a
aa1=squeeze(nmean(t,i,:)) ;
aa2=squeeze(eddy(t,i,:)) ;
aa3=squeeze(onln(t,i,:)) ;
bb=squeeze(depth(t,i,:)) ;

if isnan(bb(1))
mean_int(t,i,:) = [-1000:1:-2].*0 + NaN;
eddy_int(t,i,:) = [-1000:1:-2].*0 + NaN;
tot_int(t,i,:) = [-1000:1:-2].*0 + NaN;
else
mean_int(t,i,:) = interp1(bb,aa1,[-1000:1:-2]);
eddy_int(t,i,:) = interp1(bb,aa2,[-1000:1:-2]);
tot_int(t,i,:) = interp1(bb,aa3,[-1000:1:-2]);
end

end
end

mean_int(mean_int==0)=NaN;
eddy_int(eddy_int==0)=NaN;
tot_int(tot_int==0)=NaN;

tmean6 = squeeze(nanmean(mean_int,1));
tzmean6 = squeeze(nanmean(mean_int,[1 3]));
tpmean6 = squeeze(nanmean(mean_int,[1 2]));
zpmean6 = squeeze(nanmean(mean_int,[2 3]));

teddy6 = squeeze(nanmean(eddy_int,1));
tzeddy6 = squeeze(nanmean(eddy_int,[1 3]));
tpeddy6 = squeeze(nanmean(eddy_int,[1 2]));
zpeddy6 = squeeze(nanmean(eddy_int,[2 3]));

ttot6 = squeeze(nanmean(tot_int,1));
tztot6 = squeeze(nanmean(tot_int,[1 3]));
tptot6 = squeeze(nanmean(tot_int,[1 2]));
zptot6 = squeeze(nanmean(tot_int,[2 3]));

mmonth = str2num(datestr(date,'mm'));

list_win = find(mmonth>=9 | mmonth<=3);
list_sum = find(mmonth>=4 & mmonth<=8);
tptotwin6 = squeeze(nanmean(tot_int(list_win,:,:),[1 2]));
tptotsum6 = squeeze(nanmean(tot_int(list_sum,:,:),[1 2]));
tpmeanwin6 = squeeze(nanmean(mean_int(list_win,:,:),[1 2]));
tpmeansum6 = squeeze(nanmean(mean_int(list_sum,:,:),[1 2]));
tpeddywin6 = squeeze(nanmean(eddy_int(list_win,:,:),[1 2]));
tpeddysum6 = squeeze(nanmean(eddy_int(list_sum,:,:),[1 2]));

%% PLOT

col1 = [0.4 0.8 0.4];
col2 = [0.2 0.4 0.8];
col3 = [1.0 0.2 0.4];
col4 = [0.6 0.8 1.0];
col5 = [0.2 0.2 0.2];
col6 = [0.9 0.7 0.2];

for i=1:999
tmean1sm(:,i) = smoothdata(tmean1(:,i),'gaussian',20) ;
tmean2sm(:,i) = smoothdata(tmean2(:,i),'gaussian',20) ;
ttot1sm(:,i) = smoothdata(ttot1(:,i),'gaussian',20) ;
ttot2sm(:,i) = smoothdata(ttot2(:,i),'gaussian',20) ;
teddy1sm(:,i) = smoothdata(teddy1(:,i),'gaussian',20) ;
teddy2sm(:,i) = smoothdata(teddy2(:,i),'gaussian',20) ;
end
[asm bsm]=size(tmean1sm);
y = repmat([-1000:1:-2] , asm,1) ;
x = repmat(distc,bsm,1)' ;

%% FIGURES
link_cpt='/data/project3/kesf/tools_matlab/matlab_paths/cpt_all/gist/';
cptc = 'earth' ;

%fig = figure('position',[50 50 1100 300])
fig = figure();
%ax1=subplot (1,5,1) ;
hold on
plot(tpeddy1,[-1000:1:-2] ,'color',col3 ,'lineStyle','-','linewidth',2)
plot(tpeddy2,[-1000:1:-2] ,'color',col3 ,'lineStyle','--','linewidth',2)
plot(tpeddy3,[-1000:1:-2] ,'color',col3 ,'lineStyle',':','linewidth',2)

plot(tpeddy4,[-1000:1:-2] ,'color',col4 ,'lineStyle','-','linewidth',2)
plot(tpeddy5,[-1000:1:-2] ,'color',col4 ,'lineStyle','--','linewidth',2)
plot(tpeddy6,[-1000:1:-2] ,'color',col4 ,'lineStyle',':','linewidth',2)

%plot( tpmean1-tpmean2 , [-1000:1:-2] ,'color',col4 ,'lineStyle','-','linewidth',2)
ylabel('z (m)')
xlabel('mmol m^{-2} d^{-1}')
title('A) NH_4^+ decomp.','fontsize',14,'fontweight','normal')
set(gca,'fontsize',14)
xlim([-30 50])
ylim([-160 0])
box on
%legend('\Delta eddy','\Delta mean','location','southeast','fontname','courier','fontsize',9)
legend('50% N Red. eddy','50% N Recy. 50% Recy. eddy','50% N Recy. 90% Recy. eddy','85% N Red. eddy','85% N Red. 50% Recy. eddy','85% N Red. 90% Recy. eddy','location','best','fontsize',14)

fig = figure();
%ax1=subplot (1,5,1) ;
hold on
plot(tpmean1,[-1000:1:-2] ,'color',col3 ,'lineStyle','-','linewidth',2)
plot(tpmean2,[-1000:1:-2] ,'color',col3 ,'lineStyle','--','linewidth',2)
plot(tpmean3,[-1000:1:-2] ,'color',col3 ,'lineStyle',':','linewidth',2)

plot(tpmean4,[-1000:1:-2] ,'color',col4 ,'lineStyle','-','linewidth',2)
plot(tpmean5,[-1000:1:-2] ,'color',col4 ,'lineStyle','--','linewidth',2)
plot(tpmean6,[-1000:1:-2] ,'color',col4 ,'lineStyle',':','linewidth',2)

%plot( tpmean1-tpmean2 , [-1000:1:-2] ,'color',col4 ,'lineStyle','-','linewidth',2)
ylabel('z (m)')
xlabel('mmol m^{-2} d^{-1}')
title('A) NH_4^+ decomp.','fontsize',14,'fontweight','normal')
set(gca,'fontsize',14)
xlim([-20 40])
ylim([-160 0])
box on
%legend('\Delta eddy','\Delta mean','location','southeast','fontname','courier','fontsize',9)
legend('50% N Red. mean','50% N Recy. 50% Recy. mean','50% N Recy. 90% Recy. mean','85% N Red. mean','85% N Red. 50% Recy. mean','85% N Red. 90% Recy. mean','location','best','fontsize',14)

tpmean1_sum = sum(tpmean1(end-160:end));
tpmean2_sum = sum(tpmean2(end-160:end));
tpmean3_sum = sum(tpmean3(end-160:end));
tpmean4_sum = sum(tpmean4(end-160:end));
tpmean5_sum = sum(tpmean5(end-160:end));
tpmean6_sum = sum(tpmean6(end-160:end));

tpeddy1_sum = sum(tpeddy1(end-160:end));
tpeddy2_sum = sum(tpeddy2(end-160:end));
tpeddy3_sum = sum(tpeddy3(end-160:end));
tpeddy4_sum = sum(tpeddy4(end-160:end));
tpeddy5_sum = sum(tpeddy5(end-160:end));
tpeddy6_sum = sum(tpeddy6(end-160:end));

tptotal1 = tpmean1_sum+tpeddy1_sum
tptotal2 = tpmean2_sum+tpeddy2_sum
tptotal3 = tpmean3_sum+tpeddy3_sum
tptotal4 = tpmean4_sum+tpeddy4_sum
tptotal5 = tpmean5_sum+tpeddy5_sum
tptotal6 = tpmean6_sum+tpeddy6_sum

tpmean1_sum
tpmean2_sum
tpmean3_sum
tpmean4_sum
tpmean5_sum
tpmean6_sum
tpeddy1_sum
tpeddy2_sum
tpeddy3_sum
tpeddy4_sum
tpeddy5_sum
tpeddy6_sum

writematrix(tpeddy1,'eddy_PNDN_only_realistic.txt')
writematrix(tpeddy2,'eddy_pndn50_fixriver.txt')
writematrix(tpeddy3,'eddy_pndn90_fixriver.txt')
writematrix(tpeddy4,'eddy_FNDN_only_realistic.txt')
writematrix(tpeddy5,'eddy_fndn50_fixriver.txt')
writematrix(tpeddy6,'eddy_fndn90_fixriver.txt')

writematrix(tpmean1,'mean_PNDN_only_realistic.txt')
writematrix(tpmean2,'mean_pndn50_fixriver.txt')
writematrix(tpmean3,'mean_pndn90_fixriver.txt')
writematrix(tpmean4,'mean_FNDN_only_realistic.txt')
writematrix(tpmean5,'mean_fndn50_fixriver.txt')
writematrix(tpmean6,'mean_fndn90_fixriver.txt')

legend('50% N Red. eddy','50% N Recy. 50% Recy. eddy','50% N Recy. 90% Recy. eddy','85% N Red. eddy','85% N Red. 50% Recy. eddy','85% N Red. 90% Recy. eddy','location','best','fontsize',14)

tpmean1_05 = prctile(tpmean1,5)
tpmean2_05 = prctile(tpmean2,5)
tpmean3_05 = prctile(tpmean3,5)
tpmean4_05 = prctile(tpmean4,5)
tpmean5_05 = prctile(tpmean5,5)
tpmean6_05 = prctile(tpmean6,5)

tpmean1_95 = prctile(tpmean1,95)
tpmean2_95 = prctile(tpmean2,95)
tpmean3_95 = prctile(tpmean3,95)
tpmean4_95 = prctile(tpmean4,95)
tpmean5_95 = prctile(tpmean5,95)
tpmean6_95 = prctile(tpmean6,95)

tpeddy1_05 = prctile(tpeddy1,5)
tpeddy2_05 = prctile(tpeddy2,5)
tpeddy3_05 = prctile(tpeddy3,5)
tpeddy4_05 = prctile(tpeddy4,5)
tpeddy5_05 = prctile(tpeddy5,5)
tpeddy6_05 = prctile(tpeddy6,5)

tpeddy1_95 = prctile(tpeddy1,95)
tpeddy2_95 = prctile(tpeddy2,95)
tpeddy3_95 = prctile(tpeddy3,95)
tpeddy4_95 = prctile(tpeddy4,95)
tpeddy5_95 = prctile(tpeddy5,95)
tpeddy6_95 = prctile(tpeddy6,95)

return

ax3=subplot (1,5,2) ;
hold on
plot( tptot1-tptot2 , [-1000:1:-2] ,'color',col5 ,'lineStyle','-','linewidth',2)
ylabel('z (m)')
xlabel('mmol m^{-2} d^{-1}')
title('B) total','fontname','courier','fontsize',9,'fontweight','normal')
set(gca,'fontname','courier','fontsize',9)
xlim([-30 170])
ylim([-250 0])
box on
legend('\Delta total','location','southeast','fontname','courier','fontsize',9)

%% CLIMATOLOGY OF THE EDDY CHANGE ASSESSMENT
ax4=subplot (1,5,[3 4 5]) ;
yyaxis left
yticks([])
yyaxis right
hold on
zpeddy1cl = climato_ts(zpeddy1,date)  ;
zpeddy2cl = climato_ts(zpeddy2,date)  ;
zpmean1cl = climato_ts(zpmean1,date)  ;
zpmean2cl = climato_ts(zpmean2,date)  ;
br = bar([(zpeddy1cl-zpeddy2cl)' (zpmean1cl-zpmean2cl)'],'stacked')
br(1).FaceColor=col3;
br(2).FaceColor=col4;
xlim([0 13])
ylim([-10 65])
set(gca,'xtick',[1:12])
ylabel('mmol m^{-2} d^{-1}')
set(gca,'xtick',[1:12],'xticklabel',{'jan' 'feb' 'mar' 'apr' 'may' 'jun' 'jul' 'aug' 'sep' 'oct' 'nov' 'dec'})
xtickangle(45)
title('C) Seasonality','fontname','courier','fontsize',9,'fontweight','normal')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%set(gca,'xticklabel',[])
set(gca,'fontname','courier','fontsize',9)
box on
figure_file_name = [dir_gr,'paper_nh4_hflux_Vfinal'] ;
printeps(figure_file_name)


fig = figure('position',[50 50 500 300])
axp1=subplot(1,1,1);
pcolor(x./1000,y,ttot1sm-ttot2sm) ; shading interp
xlabel('km')
cb=colorbar;
ylabel(cb,'mmol m^{-2} d^{-1}','fontname','courier','fontsize',9)
%sgtitle('NH_4^+ 2013-2017','fontname','avenir','fontsize',14)
colormap(axp1,blu_red)
set(gca,'fontname','avenir','fontsize',9)
box on
%set(gca,'xticklabel',[])
%title('B) mean','fontname','avenir','fontsize',9,'fontweight','normal')
%set(gca,'xticklabel',[])
caxis([-500 500])
ylim([-250 0])
xlim([0 max(x(:)./1000)])
title('\Delta total F_{NH_4^+}' ,'fontweight','normal','fontsize',12 ,'fontname','courier')

text(50 , +5 ,'PL' ,'rotation',45 ,'fontweight','bold','fontsize',8 ,'fontname','courier')
text(320 , +5 ,'SP' ,'rotation',45 ,'fontweight','bold','fontsize',8 ,'fontname','courier')
text(370 , +5 ,'PD' ,'rotation',45 ,'fontweight','bold','fontsize',8 ,'fontname','courier')
text(550 , +5 ,'SB' ,'rotation',45 ,'fontweight','bold','fontsize',8 ,'fontname','courier')
figure_file_name = [dir_gr,'paper_nh4_hflux_pcolor'] ;
printeps(figure_file_name)

return

