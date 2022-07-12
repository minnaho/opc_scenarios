addpath(genpath('/data/project3/kesf/tools_matlab/matlab_paths/'))

new=1;

opcsc = 'l1617'

vname3 = 'O2';  % name of the variable, O2, Alk, DIC , DON, ORG for sum of particulate matter
%anth =0 ;
%vname = 'C' ;
%vname = 'PON' ;
vname1 = 'NO3' ;
vname2 = 'NH4' ;
%vname = 'temp' ;
%vname = 'salt' ;

dir_gr = '/data/project3/kesf/tools_matlab/applications/budget/anth_2021/graph/';
s2d = 86400 ;

list1 = 0:20:40


for vari=1:3
cpt=1
for dd = 1:length(list1)

depthmin = list1(dd) ; % shallower limit
depthmax = depthmin+20 ; % deeper limit

if vari==1
vname=vname1;
elseif vari==2
vname=vname2;
elseif vari==3
vname=vname3;
end

if new==1
%% NEW PERIOD
rep_out = '/data/project3/kesf/tools_matlab/applications/budget/anth_2021/';
fout =  [rep_out,'budget_L2_mask10_',vname,'_',num2str(depthmin),'_to_',num2str(depthmax),'_',opcsc,'.nc'];
fout0 =  [rep_out,'budget_L2_mask10_',vname,'_',num2str(depthmin),'_to_',num2str(depthmax),'_natur.nc'];
file0 = fout0 ;
time0 = double(ncread([fout0],'time')) ;
dt0 = ncread([fout0],'dt') ;
vol0 = ncread([fout0],'volume') ;
area0 = ncread([fout0],'area') ;
bgc0(:,cpt,vari) = ncread([fout0],'bgc').*s2d ./dt0 ./area0 ;
inv0(:,cpt,vari) = ncread([fout0],'invm') ./vol0 ;

file = fout ;
time = double(ncread([fout],'time')) ;
dt = ncread([fout],'dt') ;
vol = ncread([fout],'volume') ;
area = ncread([fout],'area') ;
bgc(:,cpt,vari) = ncread([fout],'bgc').*s2d ./dt ./area ;
inv(:,cpt,vari) = ncread([fout],'invm') ./vol ;

end

%% OLD PERIOD
rep_out = '/data/project3/kesf/tools_matlab/applications/budget/anth_2021/old/';
fouto =  [rep_out,'budget_L2_mask10_',vname,'_',num2str(depthmin),'_to_',num2str(depthmax),'_anthr.nc'];
fout0o =  [rep_out,'budget_L2_mask10_',vname,'_',num2str(depthmin),'_to_',num2str(depthmax),'_natur.nc'];
file0o = fout0o ;
time0o = double(ncread([fout0o],'time')) ;
dt0o = ncread([fout0o],'dt') ;
vol0o = ncread([fout0o],'volume') ;
area0o = ncread([fout0o],'area') ;
bgc0o(:,cpt,vari) = ncread([fout0o],'bgc').*s2d ./dt0o ./area0o ;
inv0o(:,cpt,vari) = ncread([fout0o],'invm') ./vol0o ;

fileo = fouto ;
timeo = double(ncread([fouto],'time')) ;
dto = ncread([fouto],'dt') ;
volo = ncread([fouto],'volume') ;
areao = ncread([fouto],'area') ;
bgco(:,cpt,vari) = ncread([fouto],'bgc').*s2d ./dto ./areao ;
invo(:,cpt,vari) = ncread([fouto],'invm') ./volo ;

cpt=cpt+1;
end % dd
end % vari

if new==1
month = str2num( datestr(time,'mm')  );
year = str2num( datestr(time,'yyyy')  );
month0 = str2num( datestr(time0,'mm')  );
year0 = str2num( datestr(time0,'yyyy')  );
time2d = repmat(time,1,8); time2d = time2d(1:46,:);
depth = -repmat([10:20:150]',1,46)';
end

montho = str2num( datestr(timeo,'mm')  );
yearo = str2num( datestr(timeo,'yyyy')  );
month0o = str2num( datestr(time0o,'mm')  );
year0o = str2num( datestr(time0o,'yyyy')  );
time2do = repmat(timeo,1,8); time2do = time2do(1:46,:);
deptho = -repmat([10:20:150]',1,46)';

disp('loading is done')

%% old
%% compare change in oxygen concentration at depth and change in biological uptake of nutrients NH4+NO3
ts1o = bgco(1:46,:,1)-bgc0o(1:46,:,1) ; % CHANGE IN BIO NO3
ts2o = bgco(1:46,:,2)-bgc0o(1:46,:,2) ; % CHANGE IN BIO NH4
%ts3o = bgco(1:140,:,3) ;
ts3o = invo(1:46,:,3)-inv0o(1:46,:,3) ; % OXYGEN CONC
ts12o = ts1o+ts2o ; % CHANGE IN BIO total DIN
vts1o = nanmean(ts1o,2); % VERTICALLY AVERAGED CHANGE BIO NO3
vts2o = nanmean(ts2o,2); % VERTICALLY AVERAGED CHANGE BIO NH4
ts12o = ts1o+ts2o ; % VERTICALLY AVERAGED CHANGE BIO TOTAL DIN
vts3o = nanmean(ts3o,2); % VERTICALLY AVERAGED CHANGE CONC O2

%ts12o(ts12o>0)=NaN; % KEEP ONLY UPTAKE
vts12o = nanmean(ts12o,2); % VERTICALLY AVERAGED UPTAKE

ts3_2o = ts3o ;
ts3_2o(ts3_2o>0)=NaN; % KEEP ONLY REDUCTION IN OXYGEN CONC. --> GENERALLY AT DEPTH ONLY, AS SURFACE SHOWS INCREASE
vts3_2o = nanmean(ts3_2o,2); % % VERTICALLY AVERAGED REDUCTION IN OXUGEN CONC negavtive only.
%vts3_2o = nanmean(ts3o,2); % % VERTICALLY AVERAGED REDUCTION IN OXUGEN CONC.

clear list12o list3o 
i=1
for year = 1997:2000
list12o(i,1) = nanmean(vts12o(month0o>0 & month0o<13 & year0o==year))  ;
list3o(i,1) = nanmean(vts3_2o(month0o>0 & month0o<13 & year0o==year))  ;
i=i+1
end

%% new
%% compare change in oxygen concentration at depth and change in biological uptake of nutrients NH4+NO3
ts1 = bgc(1:61,:,1)-bgc0(1:61,:,1) ; % CHANGE IN BIO NO3
ts2 = bgc(1:61,:,2)-bgc0(1:61,:,2) ; % CHANGE IN BIO NH4
%ts3 = bgc(1:61,:,3) ; % OXYGEN CONC
ts3 = inv(1:61,:,3)-inv0(1:61,:,3) ; % OXYGEN CONC
ts12 = ts1+ts2 ; % CHANGE IN BIO total DIN
vts1 = nanmean(ts1,2); % VERTICALLY AVERAGED CHANGE BIO NO3
vts2 = nanmean(ts2,2); % VERTICALLY AVERAGED CHANGE BIO NH4
ts12 = ts1+ts2 ; % VERTICALLY AVERAGED CHANGE BIO TOTAL DIN
vts3 = nanmean(ts3,2); % VERTICALLY AVERAGED CHANGE CONC O2

%ts12(ts12>0)=NaN; % KEEP ONLY UPTAKE
vts12 = nanmean(ts12,2); % VERTICALLY AVERAGED UPTAKE

ts3_2 = ts3 ;
ts3_2(ts3_2>0)=NaN; % KEEP ONLY REDUCTION IN OXYGEN CONC. --> GENERALLY AT DEPTH ONLY, AS SURFACE SHOWS INCREASE
vts3_2 = nanmean(ts3_2,2); % % VERTICALLY AVERAGED REDUCTION IN OXUGEN CONC negative only.
%vts3_2 = nanmean(ts3,2); % % VERTICALLY AVERAGED REDUCTION IN OXUGEN CONC.

clear list12 list3
i=1
for year = 2013:2017
list12(i,1) = nanmean(vts12(month0>0 & month0<13 & year0==year))  ;
list3(i,1) = nanmean(vts3_2(month0>0 & month0<13 & year0==year))  ;
i=i+1
end

%% merge both periods

li12 =  [list12o ; list12] ;
li3 = [list3o ; list3] ;

name = [1997:2000 2013:2017];

fig = figure('position',[100 100 550 500]) ;
hold on ; box on
tbl = table( li12(:) , li3(:) ) ;
mdl = fitlm(tbl,'linear') ;
plot(mdl)
for i=1:length(name)
text(double(li12(i)) , double(li3(i)) ,num2str(name(i)) ,'fontsize',13,'fontname','courier')
end
title(['R$^{2}$= ',num2str(mdl.Rsquared.Ordinary ,'%0.2f')] ,'fontweight','normal','Interpreter', 'latex')
legend off
ylabel('$\Delta$ Oxygen conc.','Interpreter', 'latex' ,'fontweight','bold')
xlabel('$\Delta$ NPP','Interpreter', 'latex' ,'fontweight','bold')
set(gca, 'fontsize',14,'fontname','courier')

figure_file_name = [dir_gr,'DINUPTAKE_vs_O2CONC'] ;
printeps(figure_file_name)

return

figure
hold on
tbl = table( li12(:) , li3(:) ) ;
mdl = fitlm(tbl,'linear') ;
plot(mdl)
for i=1:length(name)
text(double(li12(i)) , double(li3(i)) ,num2str(name(i)))
end

title([num2str(mdl.Rsquared.Ordinary ,'%0.2f')] ,'fontweight','normal')
legend off



return

figure ; plot(vts1+vts2 , vts3 ,'.')
figure ; plot(vts12 , vts3_2 ,'.')

i=1
for year = 1997:2000
list12(i) = nanmean(vts12(month0o>2 & month0o<6 & year0o==year))  ;
list3(i) = nanmean(vts3_2(month0o>5 & month0o<9 & year0o==year))  ;
i=i+1
end

figure
tbl = table( list12(:) , list3(:) ) ;
mdl = fitlm(tbl,'linear') ;
plot(mdl)
title([num2str(mdl.Rsquared.Ordinary ,'%0.2f')] ,'fontweight','normal')
legend off






return

profile1 = nanmean(bgco(1:46,:,1)-bgc0o(1:46,:,1),1) ;
profile2 = nanmean(bgco(1:46,:,2)-bgc0o(1:46,:,2),1) ;
profile3 = nanmean(invo(1:46,:,3)-inv0o(1:46,:,3),1) ;
%profile3 = nanmean(bgco(1:46,:,3)-bgc0o(1:46,:,3),1) ;

figure
hold on
plot([0 0] , [deptho(1) deptho(end)]  ,'--k')
pr1=plot(nanmean(profile1+profile2,1) , deptho(1,:)  ,'.-k') ;
pr2=plot(nanmean(profile1,1) , deptho(1,:)  ,'.-r') ;
pr3=plot(nanmean(profile2,1) , deptho(1,:)  ,'.-g') ;
lg = legend([pr1 pr2 pr3],{'DIN','NO3','NH4'} ,'location','southwest') ;
ylim([-150 -10])
title('ANTH-CTRL')
xlabel('mmol d^{-1}')
ylabel('z (m)')


figure
hold on
plot([0 0] , [deptho(1) deptho(end)]  ,'--k')
pr1=plot(nanmean(profile3,1) , deptho(1,:)  ,'.-k') ;
lg = legend([pr1],{'O_2'} ,'location','southwest') ;
ylim([-150 -10])
title('ANTH-CTRL')
xlabel('mmol d^{-1}')
ylabel('z (m)')


figure
plot(nanmean(profile1+profile2,1) , nanmean(profile3,1) ,'.' ,'markersize',30)




ts1 = bgco(1:46,:,1)-bgc0o(1:46,:,1) ;
ts2 = bgco(1:46,:,2)-bgc0o(1:46,:,2) ;
%ts3 = bgco(1:46,:,3)-bgc0o(1:46,:,3) ;
ts3 = invo(1:46,:,3)-inv0o(1:46,:,3) ;

vts1 = nanmean(ts1,2);
vts2 = nanmean(ts2,2);
ts12 = ts1+ts2 ;
vts3 = nanmean(ts3,2);

ts12(ts12>0)=NaN;
vts12 = nanmean(ts12,2);

ts3_2 = ts3 ;
ts3_2(ts3_2>0)=NaN;
vts3_2 = nanmean(ts3_2,2);



figure ; plot(vts1+vts2 , vts3 ,'.')
figure ; plot(vts12 , vts3_2 ,'.')


i=1
for year = 1997:2000
list12(i) = nanmean(vts12(month0o>2 & month0o<6 & year0o==year))  ;
list3(i) = nanmean(vts3_2(month0o>5 & month0o<9 & year0o==year))  ;
i=i+1
end

figure
tbl = table( list12(:) , list3(:) ) ;
mdl = fitlm(tbl,'linear') ;
plot(mdl)
title([num2str(mdl.Rsquared.Ordinary ,'%0.2f')] ,'fontweight','normal')
legend off


%%%
% new period %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%

%% biology flux

time2d = repmat(time,1,8);
depth = -repmat([0:20:140]',1,61)';

fig = figure('position' ,[50 50 1400 400])
ax1 = subplot(1,4,1);
contourf(time2d , depth , bgc0 ,30 , 'edgecolor','none')
caxis([-10 60])
datetick('x','yy')
title('CTRL')
xlim([time(1) time(end)])
cb = colorbar ;
ylim([-100 0])
ylabel(cb,'mmol d^{-1}')
colormap(ax1,parula)

ax2 = subplot(1,4,2);
contourf(time2d , depth , bgc ,30 , 'edgecolor','none')
caxis([-10 60])
datetick('x','yy')
title('ANTH')
xlim([time(1) time(end)])
cb = colorbar ;
ylim([-100 0])
ylabel(cb,'mmol d^{-1}')
colormap(ax2,parula)

ax3 = subplot(1,4,3);
contourf(time2d , depth , bgc-bgc0 ,30 , 'edgecolor','none')
caxis([-10 10])
datetick('x','yy')
xlim([time(1) time(end)])
cb = colorbar;
ylim([-100 0])
title('ANTH-CTRL')
ylabel(cb,'mmol d^{-1}')
colormap(ax3,cmocean('balance'))

ax4 = subplot(1,4,4);
hold on
plot([0 0] , [depth(1) depth(end)]  ,'--k')
plot(nanmean(bgc-bgc0,1) , depth(1,:)  ,'-k')
ylim([-100 0])
title('ANTH-CTRL')
xlabel('mmol d^{-1}')

sgtitle('Biology')

figure_file_name = [dir_gr,'NH4_budget_bio'] ;
printeps(figure_file_name)

%% inventory

fig = figure('position' ,[50 50 1400 400])
ax1 = subplot(1,4,1);
contourf(time2d , depth , inv0 ,30 , 'edgecolor','none')
caxis([0 30])
datetick('x','yy')
title('CTRL')
xlim([time(1) time(end)])
cb = colorbar ;
ylim([-100 0])
ylabel(cb,'mmol m^{-3}')
colormap(ax1,parula)

ax2 = subplot(1,4,2);
contourf(time2d , depth , inv ,30 , 'edgecolor','none')
caxis([0 30])
datetick('x','yy')
title('ANTH')
xlim([time(1) time(end)])
cb = colorbar ;
ylim([-100 0])
ylabel(cb,'mmol m^{-3}')
colormap(ax2,parula)

ax3 = subplot(1,4,3);
contourf(time2d , depth , inv-inv0 ,30 , 'edgecolor','none')
caxis([-3 3])
datetick('x','yy')
xlim([time(1) time(end)])
cb = colorbar;
ylim([-100 0])
title('ANTH-CTRL')
ylabel(cb,'mmol m^{-3}')
colormap(ax3,cmocean('balance'))

ax4 = subplot(1,4,4);
hold on
plot([0 0] , [depth(1) depth(end)]  ,'--k')
plot(nanmean(inv-inv0,1) , depth(1,:)  ,'-k')
ylim([-100 0])
title('ANTH-CTRL')
xlabel('mmol m^{-3}')

sgtitle('Concentration')

figure_file_name = [dir_gr,'NH4_budget_conc'] ;
printeps(figure_file_name)

%%%
% old period %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%

%% biology flux

time2do = repmat(timeo,1,8); time2do = time2do(1:46,:);
deptho = -repmat([10:20:150]',1,46)';

fig = figure('position' ,[50 50 1400 400])
ax1 = subplot(1,4,1);
contourf(time2do , deptho , bgc0o(1:46,1:8) ,30 , 'edgecolor','none')
caxis([-2 2])
datetick('x','yy')
title('CTRL')
xlim([timeo(1) timeo(46)])
cb = colorbar ;
ylim([-150 -10])
ylabel(cb,'mmol d^{-1}')
colormap(ax1,parula)

ax2 = subplot(1,4,2);
contourf(time2do , deptho , bgco(1:46,1:8) ,30 , 'edgecolor','none')
caxis([-2 2])
datetick('x','yy')
title('ANTH')
xlim([timeo(1) timeo(46)])
cb = colorbar ;
ylim([-150 -10])
ylabel(cb,'mmol d^{-1}')
colormap(ax2,parula)

ax3 = subplot(1,4,3);
contourf(time2do , deptho , bgco(1:46,1:8)-bgc0o(1:46,1:8) ,30 , 'edgecolor','none')
caxis([-2 2])
datetick('x','yy')
xlim([timeo(1) timeo(46)])
cb = colorbar;
ylim([-150 -10])
title('ANTH-CTRL')
ylabel(cb,'mmol d^{-1}')
colormap(ax3,cmocean('balance'))

ax4 = subplot(1,4,4);
hold on
plot([0 0] , [deptho(1) deptho(end)]  ,'--k')
plot(nanmean(bgco(1:46,1:8)-bgc0o(1:46,1:8),1) , deptho(1,:)  ,'.-k')
ylim([-150 -10])
title('ANTH-CTRL')
xlabel('mmol d^{-1}')

sgtitle('Biology')

figure_file_name = [dir_gr,'NH4_budget_bio0'] ;
printeps(figure_file_name)

%% inventory

fig = figure('position' ,[50 50 1400 400])
ax1 = subplot(1,4,1);
contourf(time2do , deptho , inv0o ,30 , 'edgecolor','none')
caxis([0 2])
datetick('x','yy')
title('CTRL')
xlim([timeo(1) timeo(46)])
cb = colorbar ;
ylim([-150 -10])
ylabel(cb,'mmol m^{-3}')
colormap(ax1,parula)

ax2 = subplot(1,4,2);
contourf(time2do , deptho , invo(1:46,:) ,30 , 'edgecolor','none')
caxis([0 2])
datetick('x','yy')
title('ANTH')
xlim([timeo(1) timeo(46)])
cb = colorbar ;
ylim([-150 -10])
ylabel(cb,'mmol m^{-3}')
colormap(ax2,parula)

ax3 = subplot(1,4,3);
contourf(time2do , deptho , invo(1:46,:)-inv0o ,30 , 'edgecolor','none')
caxis([-2 2])
datetick('x','yy')
xlim([timeo(1) timeo(46)])
cb = colorbar;
ylim([-150 -10])
title('ANTH-CTRL')
ylabel(cb,'mmol m^{-3}')
colormap(ax3,cmocean('balance'))

ax4 = subplot(1,4,4);
hold on
plot([0 0] , [deptho(1) deptho(end)] ,'--k')
plot( nanmean(invo(1:46,:)-inv0o,1) ,  deptho(1,:) ,'.-k' )
ylim([-150 -10])
title('ANTH-CTRL')
xlabel('mmol m^{-3}')

sgtitle('Concentration')

figure_file_name = [dir_gr,'NH4_budget_conc0'] ;
printeps(figure_file_name)





%%
%%
%%


return


fig = figure('position' ,[50 50 1200 400])
ax1 = subplot(1,3,1);
contourf(time2d , depth , phys0 ,30 , 'edgecolor','none')
caxis([-20 20])
datetick('x','yy')
xlim([time(1) time(end)])
cb = colorbar ;
ylim([-100 0])
ylabel(cb,'mmol d^{-1}')
title('CTRL')
colormap(ax1,parula)

ax2 = subplot(1,3,2);
contourf(time2d , depth , phys ,30 , 'edgecolor','none')
caxis([-20 20])
datetick('x','yy')
xlim([time(1) time(end)])
cb = colorbar ;
ylim([-100 0])
ylabel(cb,'mmol d^{-1}')
title('ANTH')
colormap(ax2,parula)

ax3 = subplot(1,3,3);
contourf(time2d , depth , phys-phys0 ,30 , 'edgecolor','none')
caxis([-10 10])
datetick('x','yy')
title('ANTH-CTRL')
xlim([time(1) time(end)])
cb = colorbar;
ylim([-100 0])
ylabel(cb,'mmol d^{-1}')
colormap(ax3,cmocean('balance'))
sgtitle('physics')



return





fig = figure('position' ,[50 50 1400 400])
ax1 = subplot(1,3,1);
contourf(time2d , depth , phys0+bgc0 ,30 , 'edgecolor','none')
caxis([-20 20])
datetick('x','yy')
title('CTRL')
xlim([time(1) time(end)])
cb = colorbar ;
ylim([-100 0])
ylabel(cb,'mmol d^{-1}')
colormap(ax1,parula)

ax2 = subplot(1,3,2);
contourf(time2d , depth , phys+bgc ,30 , 'edgecolor','none')
caxis([-20 20])
datetick('x','yy')
title('ANTH')
xlim([time(1) time(end)])
cb = colorbar ;
ylim([-100 0])
ylabel(cb,'mmol d^{-1}')
colormap(ax2,parula)

ax3 = subplot(1,3,3);
contourf(time2d , depth , (phys+bgc)-(phys0+bgc0) ,30 , 'edgecolor','none')
caxis([-3 3])
datetick('x','yy')
title('ANTH')
xlim([time(1) time(end)])
cb = colorbar ;
ylim([-100 0])
ylabel(cb,'mmol d^{-1}')
colormap(ax3,cmocean('balance'))
sgtitle('total balance')




%%



return



fig = figure('position',[100 100 1000 400]) ;
hold on
plot(time,inv-inv0,'k' , 'linewidth',1.5)
datetick('x','yyyy')
title([vname,' ',num2str(depthmin),' to ',num2str(depthmax),' anthr'])
grid on
grid minor

%end
return

fig = figure('position',[100 100 1000 400]) ;
hold on
plot(time,inv0,'k' , 'linewidth',1.5)
plot(time,inv,'r' , 'linewidth',1.5)
legend('ctrl','anth')
datetick('x','yyyy')

fig = figure('position',[100 100 1000 400]) ;
hold on
plot(time,inv-inv0,'k' , 'linewidth',1.5)
datetick('x','yyyy')
title([vname,' ',num2str(depthmin),' to ',num2str(depthmax),' anthr'])


fig = figure('position',[100 100 1000 400]) ;
hold on
plot(time,bgc-bgc0,'g' , 'linewidth',1.5)
plot(time,phys-phys0,'r' , 'linewidth',1.5)
plot(time,dndt-dndt0,'c' , 'linewidth',1.5)
plot(time,(phys+bgc-dndt)-(phys0+bgc0-dndt0) ,'k' , 'linewidth',1.5)
datetick('x','yyyy')
title([vname,' ',num2str(depthmin),' to ',num2str(depthmax),' anthr'])

figure
subplot(1,3,[1 2])
bar([nanmean(bgc-bgc0) nanmean(phys-phys0) nanmean(bgc-bgc0)+nanmean(phys-phys0)] )
set(gca, 'XTickLabel', {'bgc' 'phys' 'dO_{2}/dt'})
subplot(1,3,[1 2])
bar([nanmean(bgc-bgc0) nanmean(phys-phys0) nanmean(bgc-bgc0)+nanmean(phys-phys0)] )
set(gca, 'XTickLabel', {'bgc' 'phys' 'dO_{2}/dt'})

%end
return

mean(JO2(1:34))

return
nitrif = ncread([fout],'nitrification').*s2d ./dt ./area ;
denitrif = ncread([fout],'denitrification').*s2d ./dt ./area ;
uptake_diat = ncread([fout],'uptake_diat').*s2d ./dt ./area ;
uptake_sp = ncread([fout],'uptake_sp').*s2d ./dt ./area ;
biofix = uptake_diat+uptake_sp ;
merid = ncread([fout],'meridional').*s2d./dt ./area ;
zonal = ncread([fout],'zonal').*s2d./dt ./area ;
vert = ncread([fout],'vertical_flux').*s2d./dt ./area ;
dndt = ncread([fout],'dNdT').*s2d./dt ./area ;
inv1 = ncread([fout],'inv1') ;
tot = bgc + merid + zonal + vert - dndt ;


%% winter : 11-01
list1 = find(month>=11 | month==1) ;
%% upwelling: 02-05
list2 = find(month>=2 & month<=5) ;
%% summer: 06-10
list3 = find(month>=6 & month<=10) ;

%% inventory
var = inv1./1e3 ;
[nanmean(var(list1)) nanmean(var(list2)) nanmean(var(list3))]/1e6
%% bio fixation
var = biofix ;
[nanmean(var(list1)) nanmean(var(list2)) nanmean(var(list3))]
%% nitrification
var =  nitrif;
[nanmean(var(list1)) nanmean(var(list2)) nanmean(var(list3))]
%% sediment
var = denitrif ;
[nanmean(var(list1)) nanmean(var(list2)) nanmean(var(list3))]
%% bgc
var = bgc ;
[nanmean(var(list1)) nanmean(var(list2)) nanmean(var(list3))]
%% physical flux
phys = zonal+merid+vert ;
var = phys ;
[nanmean(var(list1)) nanmean(var(list2)) nanmean(var(list3))]
%% dndt
var = dndt ;
[nanmean(var(list1)) nanmean(var(list2)) nanmean(var(list3))]


figure
hold on
plot(time,bgc,'g')
plot(time,phys,'b')
plot(time,dndt,'c')
plot(time,phys+bgc-dndt ,'k')
datetick('x','mm')


figure
hold on
plot(time,nitrif,'-og')
plot(time,-denitrif,'--g')
plot(time,-biofix,'-xg')
plot(time,bgc,'-g')
datetick('x','mm')

return

cpt = 1
for i=1:12
list = find(months==i | year~=1999) ;
diag_tot(cpt) = nanmean(tot(list),1) ;
diag_vert(cpt) = nanmean(vert(list),1) ;
diag_zonal(cpt) = nanmean(zonal(list),1) ;
diag_merid(cpt) = nanmean(merid(list),1) ;
diag_dndt(cpt) = nanmean(dndt(list),1) ;
diag_bgc(cpt) = nanmean(bgc(list),1) ;
cpt=cpt+1;
end % i
clear cpt
timec = 1:12;

%% CLIMATOLOGIES
ind = 1 ;
clim_diag_tot = diag_tot ;
clim_vert = diag_vert ;
clim_zonal = diag_zonal ;
clim_merid = diag_merid ;
clim_dndt = diag_dndt ;
clim_bgc = diag_bgc ;
%dt = [31*86400 28*86400 31*86400 30*86400 31*86400 30*86400 31*86400 31*86400 30*86400 31*86400 30*86400 31*86400];

return

s2d = 86400 ;
%% total
domask_LTER
fig  = figure('visible','on','position',[0 0 800 800]);
subplot 111
                hold on
                m_proj('mercator','long',[-120.6 -119],'lat',[33.8 34.8]);
%                m_proj('mercator','long',[minlon maxlon],'lat',[minlat maxlat]);
                m_pcolor(lon,lat,mask) ; shading flat
[c m] =         m_contour(lon,lat,-abs(h) ,[-50 -500 -1000] ,'k' ) ;
                clabel(c,m,'FontSize',9,'Color','red','LabelSpacing',400, 'FontName', 'Courier')
                m_gshhs_h('patch',[.9 .9 .9]);
                m_grid('linewi',1,'tickdir','out','FontSize',12,'xtick',3, 'FontName', 'Courier');
                m_ruler([.5 .9],.9,1,'fontsize',10, 'FontName', 'Courier');
		m_text( -120.45,34.52 , 'Conception. Pt.','fontsize',10, 'FontName', 'Courier')
		m_text( -119.47,34.42 , 'Rincon Pt.','fontsize',10, 'FontName', 'Courier')
		title('LTER (Santa Barbara)')

figure_file_name = [dir_gr,'Location'] ;
printpng(figure_file_name)


fig  = figure('visible','on','position',[0 0 800 800]);
subplot 111
hold on
grid on
box on
plot(time, bgc ,'.-g','linewidth',2)
plot(time, (zonal+merid+vert) ,'.-b','linewidth',2)
%plot(timec, s2d.*clim_vert./dt ,'.-r','linewidth',2)
plot(time, dndt ,'.-c','linewidth',2)
plot(time, tot ,'.-k','linewidth',2)
legend('Biology','Physics','dN/dt','0') %,'Total')
title(['Convergence of NO_3'])
ylabel('mmol N m^{-3} d^{-1}')
datetick('x')
%xlim([0 15])
%set(gca,'xtick',([1 2 3 4 5 6 7 8 9 10 11 12]),'xticklabel',['01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12'],'tickdir','out');
set(gca,'fontname','Courier','fontsize',12)
xlabel('Month')

figure_file_name = [dir_gr,'NH4_convergence_L1'] ;
printpng(figure_file_name)


bgc_list = bgc(24:35) ; 
nanmean(bgc_list(2:5))
nanmean(bgc_list(6:10))
nanmean(bgc_list(11:12))





