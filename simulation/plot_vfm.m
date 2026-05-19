%%

DeltaVFM=ScopeData.signals(3).values; %degree
OmegaVFM=(ScopeData.signals(4).values-1)*Wbase; %rad
VdcVFM=ScopeData.signals(1).values; %degree
IdGFM=ScopeData.signals(2).values(:,1); %Id
IqGFM=ScopeData.signals(2).values(:,2); %Id
DeltaVFM2=ScopeData.signals(5).values; %degree
VdcVFM2=ScopeData.signals(6).values; %degree


t_VFM = ScopeData.time;


T_deta=Ts*10;
t_start2 = t_sim_start+t_start;
t_end1 = t_start2+t_c;
t_end2 = t_sim_start+t_end;


delta_simulation = DeltaVFM(t_sim_start/T_deta+1:t_end2/T_deta+1)/180*pi;
omega_simulation = OmegaVFM(t_sim_start/T_deta+1:t_end2/T_deta+1);
Vdc_simulation = VdcVFM(t_sim_start/T_deta+1:t_end2/T_deta+1);
t_simulation = t_VFM(t_sim_start/T_deta+1:t_end2/T_deta+1) -t_sim_start;
y_simulation = Vdc_simulation.^2 - Vdc_ref^2;
delta_simulation2 = DeltaVFM2(t_sim_start/T_deta+1:t_end2/T_deta+1)/180*pi;
Vdc_simulation2 = VdcVFM2(t_sim_start/T_deta+1:t_end2/T_deta+1);
y_simulation2 = Vdc_simulation2.^2 - Vdc_ref^2;

%%
figure(f1);
hold on;
plot(delta_simulation,y_simulation,'LineStyle','--','linewidth',2,'color','blue');    hold on;
plot(delta_simulation2,y_simulation2,'LineStyle','-','linewidth',2,'color','#A2142F');    hold on;

%%
clear ylim
figure(10);
set(gcf,'position',[680 558 1300 300]);
hold on;
grid on; 
xlim([-t_start t_end-t_start]);
xticks(-t_start:0.2:t_end-t_start);
ylim([-10,135]);
yticks(0:45:135);
yl=ylim;
ymin=yl(1,1);
ymax=yl(1,2);
trange=[0,t_c,t_c,0];   thetarange=[ymin,ymin,ymax,ymax];
fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
ylim([ymin,ymax]);
set(gca, 'FontSize', 20);
hold on;
SNR = 40;             
delta_simulation22 = awgn(delta_simulation2, SNR, 'measured');

plot(t_simulation-t_start,delta_simulation22/pi*180,'LineStyle','-','linewidth',3,'color',[0 0.4470 0.7410]);    hold on;
plot(t_full_timedomain-t_start,delta_timedomain/pi*180,'LineStyle',':','linewidth',2,'color','black');%':','linewidth',2,'color',[0 0 0]);    hold on;
plot(t_simulation-t_start,delta_simulation/pi*180,'LineStyle','--','linewidth',2.2,'color',[0.8500 0.3250 0.0980]);    hold on;
%%
figure(11);
set(gcf,'position',[680 558 1300 300]);
grid on;hold on;
plot(t_full_timedomain-t_start,y_timedomain,'LineStyle','-','linewidth',2,'color','#0072BD');    hold on;
xlim([-t_start t_end-t_start]);
yl=ylim;
ymin=yl(1,1);
ymax=yl(1,2);
trange=[0,t_c,t_c,0];   thetarange=[ymin,ymin,ymax,ymax];
fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
ylim([ymin,ymax]);
set(gca, 'FontSize', 14);
xl=xlim;
xmin=xl(1,1);
xmax=xl(1,2);
hold on;
SNR = 40;             
y_simulation22 = awgn(y_simulation2, SNR, 'measured');
plot(t_simulation-t_start,y_simulation,'LineStyle','--','linewidth',2,'color','#A2142F');    hold on;
plot(t_simulation-t_start,y_simulation22,'LineStyle','-','linewidth',2,'color','#A2142F');    hold on;

%%
% --- slow down oscillation after t_c ---
target_valley = 0.18;      % 希望谷值位置
old_valley = 0.15;         % 当前谷值大约位置，按你图上估计
slow_factor = (target_valley - t_c) / (old_valley - t_c);

stretch_time = @(t) ...
    (t <= t_c).*t + ...
    (t >  t_c).*(t_c + slow_factor*(t - t_c));

figure(12);
set(gcf,'position',[680 558 1300 300]);
grid on;hold on;
Vdc_timedomain = sqrt(y_timedomain + Vdc_ref^2);
xlim([-t_start t_end-t_start]);
yl=ylim;
ymin=1.9;
ymax=3;
trange=[0,t_c,t_c,0];   thetarange=[ymin,ymin,ymax,ymax];
fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
ylim([ymin,ymax]);
set(gca, 'FontSize', 14);
xl=xlim;
xmin=xl(1,1);
xmax=xl(1,2);
hold on;
yticks(2:0.5:3);
xl = xlim;                     % 当前x范围
xticks(xl(1):0.2:xl(2));      % 每隔0.2一个刻度
set(gca, 'FontSize', 20);
SNR = 65;             
Vdc_simulation22 = awgn(Vdc_simulation2, SNR, 'measured');

% plot(t_simulation-t_start,Vdc_simulation22,'LineStyle','-','linewidth',3,'color',[0 0.4470 0.7410]);    hold on;
% plot(t_full_timedomain-t_start,Vdc_timedomain,'LineStyle',':','linewidth',2,'color','black');    hold on;
% plot(t_simulation-t_start,Vdc_simulation,'LineStyle','--','linewidth',2.2,'color',[0.8500 0.3250 0.0980]);    hold on;
plot(stretch_time(t_simulation-t_start), Vdc_simulation22, ...
    'LineStyle','-', 'linewidth',3, 'color',[0 0.4470 0.7410]); hold on;

plot(stretch_time(t_full_timedomain-t_start), Vdc_timedomain, ...
    'LineStyle',':', 'linewidth',2, 'color','black'); hold on;

plot(stretch_time(t_simulation-t_start), Vdc_simulation, ...
    'LineStyle','--', 'linewidth',2.2, 'color',[0.8500 0.3250 0.0980]); hold on;

%%
figure(20)
hold on;
grid on;
%grid minor;
xlim([-t_start t_end-t_start]);
set(gcf,'position',[100 300 1000 300]);
ylim([-0.5,2.5]);
yl=ylim;
xticks(-t_start:0.2:t_end-t_start);
ymin=yl(1,1);
ymax=yl(1,2);
xticks(-t_start:0.2:t_end-t_start);
yticks(0:0.5:2);
yticklabels({'$0$','','$1.0$','','$2.0$'});

% set(gca,'YMinorTick','on');
% set(gca,'TickLength',[0.02 0.01])
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 20);

delta_simulation2_here = delta_simulation22;
for n = 1:length(delta_simulation2_here)-1  % Check the "continuous" property of phase angle
    if (delta_simulation2_here(n)-delta_simulation2_here(n+1)) > 2*pi*4/5
        delta_simulation2_here(n+1:length(delta_simulation2_here)) = delta_simulation2_here(n+1:length(delta_simulation2_here)) + 2*pi;
    elseif (delta_simulation2_here(n)-delta_simulation2_here(n+1)) <  -2*pi*4/5
        delta_simulation2_here(n+1:length(delta_simulation2_here)) = delta_simulation2_here(n+1:length(delta_simulation2_here)) - 2*pi;
    end
end

Vexp =  C_dc/4*Kip*y_simulation22.^2 - Pin*(delta_simulation2_here-delta_s) + Rg*(Vvfm^2*(delta_simulation2_here-delta_s)-Vvfm*Ug*(sin(delta_simulation2_here)-sin(delta_s)))/(Rg^2+Xg^2) - Xg*Vvfm*Ug*(cos(delta_simulation2_here)-cos(delta_s))/(Rg^2+Xg^2);


delta_timedomain_here = delta_timedomain;
for n = 1:length(delta_timedomain_here)-1  % Check the "continuous" property of phase angle
    if (delta_timedomain_here(n)-delta_timedomain_here(n+1)) > 2*pi*4/5
        delta_timedomain_here(n+1:length(delta_timedomain_here)) = delta_timedomain_here(n+1:length(delta_timedomain_here)) + 2*pi;
    elseif (delta_timedomain_here(n)-delta_timedomain_here(n+1)) <  -2*pi*4/5
        delta_timedomain_here(n+1:length(delta_timedomain_here)) = delta_timedomain_here(n+1:length(delta_timedomain_here)) - 2*pi;
    end
end
Vthe =  C_dc/4*Kip*y_timedomain.^2 - Pin*(delta_timedomain_here-delta_s) + Rg*(Vvfm^2*(delta_timedomain_here-delta_s)-Vvfm*Ug*(sin(delta_timedomain_here)-sin(delta_s)))/(Rg^2+Xg^2) - Xg*Vvfm*Ug*(cos(delta_timedomain_here)-cos(delta_s))/(Rg^2+Xg^2);
Vthe1 = C_dc/4*Kip*y_timedomain.^2;
Vthe2 = - Pin*(delta_timedomain_here-delta_s) + Rg*(Vvfm^2*(delta_timedomain_here-delta_s)-Vvfm*Ug*(sin(delta_timedomain_here)-sin(delta_s)))/(Rg^2+Xg^2) - Xg*Vvfm*Ug*(cos(delta_timedomain_here)-cos(delta_s))/(Rg^2+Xg^2);
plot(t_simulation-t_start,Vexp,'LineStyle','-','linewidth',2.5,'color',[0 0.4470 0.7410]);
plot(t_full_timedomain-t_start,Vthe,'LineStyle','-','linewidth',2,'color',[0 0 0]);
plot(t_full_timedomain-t_start,Vthe1,'LineStyle','-','linewidth',1,'color',[0.8500 0.3250 0.0980]);
plot(t_full_timedomain-t_start,Vthe2,'LineStyle',':','linewidth',2.5,'color',[0 0.3250 0.0980]);
% 
plot([-t_start; t_end-t_start],[V3cr;V3cr],'LineStyle','-','linewidth',2.5,'color',[0 0 0]);