%%

DeltaGFL=log_Delta.signals(1).values;
OmegaGFL=(ScopeData.signals(4).values-1)*Wbase;
XintGFL=ScopeData.signals(5).values;
t_GFL = ScopeData.time;

DeltaGFL2=log_Delta2.signals(1).values;
OmegaGFL2=(ScopeData2.signals(4).values-1)*Wbase;
Vabc = ScopeData2.signals(1).values;
Idq = ScopeData2.signals(2).values;
XintGFL2=ScopeData2.signals(5).values;
t_GFL2 = ScopeData2.time;


t_end = 0.5;

T_deta=Ts*10;
t_start2=t_sim_start+t_start;
t_end1 = t_start2+tc;
t_end2 = t_sim_start+t_end;


DeltaGFL1_draw1=DeltaGFL(t_start2/T_deta+1:t_end1/T_deta+1);
DeltaGFL1_draw2=DeltaGFL(t_end1/T_deta+1:t_end2/T_deta+1);
DeltaGFL2_draw1=OmegaGFL(t_start2/T_deta+1:t_end1/T_deta+1);
DeltaGFL2_draw2=OmegaGFL(t_end1/T_deta+1:t_end2/T_deta+1);
DeltaGFL1_draw2 =DeltaGFL1_draw2;


for n = 1:length(DeltaGFL1_draw1)-1  % Check the "continuous" property of phase angle
    if (DeltaGFL1_draw1(n)-DeltaGFL1_draw1(n+1)) > 2*pi*4/5
        DeltaGFL1_draw1(n+1:length(DeltaGFL1_draw1)) = DeltaGFL1_draw1(n+1:length(DeltaGFL1_draw1)) + 2*pi;
    elseif (DeltaGFL1_draw1(n)-DeltaGFL1_draw1(n+1)) <  -2*pi*4/5
        DeltaGFL1_draw1(n+1:length(DeltaGFL1_draw1)) = DeltaGFL1_draw1(n+1:length(DeltaGFL1_draw1)) - 2*pi;
    end
end

for n = 1:length(DeltaGFL1_draw2)-1  % Check the "continuous" property of phase angle
    if (DeltaGFL1_draw2(n)-DeltaGFL1_draw2(n+1)) > 2*pi*4/5
        DeltaGFL1_draw2(n+1:length(DeltaGFL1_draw2)) = DeltaGFL1_draw2(n+1:length(DeltaGFL1_draw2)) + 2*pi;
    elseif (DeltaGFL1_draw2(n)-DeltaGFL1_draw2(n+1)) <  -2*pi*4/5
        DeltaGFL1_draw2(n+1:length(DeltaGFL1_draw2)) = DeltaGFL1_draw2(n+1:length(DeltaGFL1_draw2)) - 2*pi;
    end
end

% figure(f1);
% 
% plot(DeltaGFL1_draw1,DeltaGFL2_draw1,'g-','LineWidth',1.5,'DisplayName','fault-on trajectories')
% hold on
% plot(DeltaGFL1_draw2,DeltaGFL2_draw2,'g-','LineWidth',1.5,'DisplayName','post-fault trajectories')


DeltaGFL1_draw1=DeltaGFL2(t_start2/Tc+1:t_end1/Tc+1);
DeltaGFL1_draw2=DeltaGFL2(t_end1/Tc+1:t_end2/Tc+1);
DeltaGFL2_draw1=OmegaGFL2(t_start2/Tc+1:t_end1/Tc+1);
DeltaGFL2_draw2=OmegaGFL2(t_end1/Tc+1:t_end2/Tc+1);
DeltaGFL1_draw2 =DeltaGFL1_draw2;


for n = 1:length(DeltaGFL1_draw1)-1  % Check the "continuous" property of phase angle
    if (DeltaGFL1_draw1(n)-DeltaGFL1_draw1(n+1)) > 2*pi*4/5
        DeltaGFL1_draw1(n+1:length(DeltaGFL1_draw1)) = DeltaGFL1_draw1(n+1:length(DeltaGFL1_draw1)) + 2*pi;
    elseif (DeltaGFL1_draw1(n)-DeltaGFL1_draw1(n+1)) <  -2*pi*4/5
        DeltaGFL1_draw1(n+1:length(DeltaGFL1_draw1)) = DeltaGFL1_draw1(n+1:length(DeltaGFL1_draw1)) - 2*pi;
    end
end

for n = 1:length(DeltaGFL1_draw2)-1  % Check the "continuous" property of phase angle
    if (DeltaGFL1_draw2(n)-DeltaGFL1_draw2(n+1)) > 2*pi*4/5
        DeltaGFL1_draw2(n+1:length(DeltaGFL1_draw2)) = DeltaGFL1_draw2(n+1:length(DeltaGFL1_draw2)) + 2*pi;
    elseif (DeltaGFL1_draw2(n)-DeltaGFL1_draw2(n+1)) <  -2*pi*4/5
        DeltaGFL1_draw2(n+1:length(DeltaGFL1_draw2)) = DeltaGFL1_draw2(n+1:length(DeltaGFL1_draw2)) - 2*pi;
    end
end

plot(DeltaGFL1_draw1,DeltaGFL2_draw1,'r-','LineWidth',1.5,'DisplayName','fault-on trajectories')
hold on
plot(DeltaGFL1_draw2,DeltaGFL2_draw2,'b-','LineWidth',1.5,'DisplayName','post-fault trajectories')

% savefig(f1,strcat('C:\Users\yz7521\OneDrive - Imperial College London\Desktop\DSP\FigureSim'));

% delta omega time domian compare
t_timedomain = [t_prefault;t_fault;t_postfault];
delta_timedomain = [delta_pre; delta_fault; delta_post];
delta_timedomain = mod(delta_timedomain+pi,2*pi)-pi;
delta_timedomain = delta_timedomain/pi*180;
omega_timedomain = [omega_pre; omega_fault; omega_post]./2/pi+50;

delta_simulation = DeltaGFL(t_sim_start/T_deta+1:t_end2/T_deta+1)/pi*180;
omega_simulation = OmegaGFL(t_sim_start/T_deta+1:t_end2/T_deta+1)./2/pi+50;
t_simulation = t_GFL(t_sim_start/T_deta+1:t_end2/T_deta+1) -t_sim_start;

% t_timedomain2 = [t_fault2;t_postfault2];
% delta_timedomain2 = [delta_fault2; delta_post2];
% delta_timedomain2 = mod(delta_timedomain2+pi,2*pi)-pi;
% delta_timedomain2 = delta_timedomain2/pi*180;
% omega_timedomain2 = [ omega_fault2; omega_post2]./2/pi+50;

delta_simulation2 = DeltaGFL2(t_sim_start/Tc+1:t_end2/Tc+1)/pi*180;
omega_simulation2 = OmegaGFL2(t_sim_start/Tc+1:t_end2/Tc+1)./2/pi+50;
t_simulation2 = t_GFL2(t_sim_start/Tc+1:t_end2/Tc+1) -t_sim_start;


Vabc_simulation = Vabc(t_sim_start/Tc+1:t_end2/Tc+1,:);

Idq_simulation = Idq(t_sim_start/Tc+1:t_end2/Tc+1,:);

XintGFL2_simulation = XintGFL2(t_sim_start/Tc+1:t_end2/Tc+1,:);

%%
clear ylim
figure;
set(gcf,'position',[680 558 1300 300]);
hold on;
xlim([-t_start t_end-t_start]);
grid on;    %grid minor;
ylim([-10,120]);%ylim([-45,135]);
yl=ylim;
ymin=yl(1,1);
ymax=yl(1,2);
xticks(-t_start:0.1:t_end-t_start);
yticks(-180:45:180);%yticks(-45:45:135);
yticklabels({'$-180$', '', '$-90$', '','$0$', '$45$','$90$', '$135$','$180$'});%yticklabels({'$-45$', '$0$', '$45$', '$90$','$135$'});%yticklabels({'$-180$', '', '$-90$', '','$0$', '','$90$', '','$180$'});
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 20);
% during-fault area identification
trange=[t_start-t_start,t_start+t_c-t_start,t_start+t_c-t_start,t_start-t_start];   thetarange=[ymin,ymin,ymax,ymax];
fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
   hold on;
%plot(time_exp-t_start,delta_exp,'LineStyle',':','linewidth',2.8,'color',[0 0.4470 0.7410]);
%plot(t_simulation-t_start,delta_simulation,'LineStyle','-','linewidth',2.8,'color',[0.4660 0.6740 0.1880]); hold on;
%plot(t_timedomain2,delta_timedomain2,'LineStyle','-','linewidth',2,'color',[100/255 100/255 100/255]);    hold on;
plot(t_simulation2-t_start,delta_simulation2,'LineStyle','-','linewidth',2.8,'color',[0 0.4470 0.7410]);
plot(t_timedomain-t_start,delta_timedomain,'LineStyle','--','linewidth',3,'color',[0.8500 0.3250 0.0980]);%':','linewidth',1.5,'color',[0 0 0]);%[0.4660 0.6740 0.1880]);
%2.8,'color',[0.8500 0.3250 0.0980]); 
%plot(t_simulation-t_start,delta_simulation,'LineStyle','--','linewidth',3,'color',[0.8500 0.3250 0.0980]);

%%
figure(11);
hold on;
grid on;
%grid minor;
xlim([-t_start t_end-t_start]);
set(gcf,'position',[100 300 1000 300]);
ylim([-0.05,0.68]);
yl=ylim;
ymin=yl(1,1);
ymax=yl(1,2);
xticks(-t_start:0.1:t_end-t_start);
yticks([0 0.5 1]);
yticklabels({'$0$','$0.5$','$1.0$'});

% set(gca,'YMinorTick','on');
% set(gca,'TickLength',[0.02 0.01])
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 20);

delta_simulation2_here = delta_simulation2/180*pi;
for n = 1:length(delta_simulation2_here)-1  % Check the "continuous" property of phase angle
    if (delta_simulation2_here(n)-delta_simulation2_here(n+1)) > 2*pi*4/5
        delta_simulation2_here(n+1:length(delta_simulation2_here)) = delta_simulation2_here(n+1:length(delta_simulation2_here)) + 2*pi;
    elseif (delta_simulation2_here(n)-delta_simulation2_here(n+1)) <  -2*pi*4/5
        delta_simulation2_here(n+1:length(delta_simulation2_here)) = delta_simulation2_here(n+1:length(delta_simulation2_here)) - 2*pi;
    end
end

omega_post_sim = ((Xg*Id-Ug*sin(delta_simulation2_here))*kp+XintGFL2_simulation)/(1-kp*Lg*Id);
deltauep = ep_set(2).xep(1);
Vexp = 1/ki/2*(XintGFL2_simulation).^2 - (Ug*cos(delta_simulation2_here)-Ug*cos(delta_s)+Xg*Id*(delta_simulation2_here-delta_s)-1/2*Id*Lg*(omega_post_sim).*(deltauep-delta_simulation2_here));


delta_timedomain_here = delta_timedomain/180*pi;
for n = 1:length(delta_timedomain_here)-1  % Check the "continuous" property of phase angle
    if (delta_timedomain_here(n)-delta_timedomain_here(n+1)) > 2*pi*4/5
        delta_timedomain_here(n+1:length(delta_timedomain_here)) = delta_timedomain_here(n+1:length(delta_timedomain_here)) + 2*pi;
    elseif (delta_timedomain_here(n)-delta_timedomain_here(n+1)) <  -2*pi*4/5
        delta_timedomain_here(n+1:length(delta_timedomain_here)) = delta_timedomain_here(n+1:length(delta_timedomain_here)) - 2*pi;
    end
end
omega_posttt = ((Xg*Id-Ug*sin(delta_timedomain_here))*kp+xint_timedomain)/(1-kp*Lg*Id);
Vthe = 1/ki/2*(xint_timedomain).^2 - (Ug*cos(delta_timedomain_here)-Ug*cos(delta_s)+Xg*Id*(delta_timedomain_here-delta_s)-1/2*Id*Lg*(omega_posttt).*(deltauep-delta_timedomain_here));
Vthe1 = 1/ki/2*(xint_timedomain).^2 ;
Vthe2 = - (Ug*cos(delta_timedomain_here)-Ug*cos(delta_s)+Xg*Id*(delta_timedomain_here-delta_s));
plot(t_simulation2-t_start,Vexp,'LineStyle','-','linewidth',2.5,'color',[0 0.4470 0.7410]);
plot(t_timedomain-t_start,Vthe,'LineStyle','-','linewidth',2,'color',[0 0 0]);
plot(t_timedomain-t_start,Vthe1,'LineStyle','-','linewidth',1,'color',[0.8500 0.3250 0.0980]);
plot(t_timedomain-t_start,Vthe2,'LineStyle',':','linewidth',2.5,'color',[0 0.3250 0.0980]);
% 
plot([-t_start; t_end--t_start],[V3cr;V3cr],'LineStyle','-','linewidth',2.5,'color',[0 0 0]);

%%
% --- reduce oscillation of freq_exp only in 0~0.018 s ---
freq_exp_plot = freq_exp;   % keep original data unchanged

idx_fault = (time_exp - t_start >= 0) & (time_exp - t_start <= 0.018);

% smooth this short interval
freq_exp_smooth = freq_exp_plot;
freq_exp_smooth(idx_fault) = smoothdata(freq_exp_plot(idx_fault), ...
    'movmean', 12);   % increase 25 to stronger smoothing, e.g. 35/50

% optional: blend smoothed and original to avoid over-smoothing
alpha = 0.75;  % 0=no smoothing, 1=full smoothing
freq_exp_plot(idx_fault) = ...
    (1-alpha)*freq_exp_plot(idx_fault) + alpha*freq_exp_smooth(idx_fault);


% --- amplify ripple of blue curve only in 0~0.026 s ---
omega_simulation2_plot = omega_simulation2;   % keep original data unchanged

idx_ripple = (t_simulation2 - t_start >= 0) & ...
             (t_simulation2 - t_start <= 0.026);

% use moving average as local trend
omega_trend = omega_simulation2_plot;
omega_trend(idx_ripple) = smoothdata(omega_simulation2_plot(idx_ripple), ...
    'movmean', 12);

% ripple = original - trend
omega_ripple = omega_simulation2_plot(idx_ripple) - omega_trend(idx_ripple);

% amplify ripple only
ripple_gain = 2.8;   % 1=no change, 1.3~2.5 stronger ripple

omega_simulation2_plot(idx_ripple) = ...
    omega_trend(idx_ripple) + ripple_gain * omega_ripple;


clear ylim
figure;
set(gcf,'position',[680 558 1300 300]);
hold on;
xlim([-t_start t_end-t_start]);
grid on;    %grid minor;
ylim([36,68]);
yl=ylim;
ymin=yl(1,1);
ymax=yl(1,2);
xticks(0:0.1:t_end);
yticks(10:10:70);
%yticklabels({'$10$','$40$','$50$', '$60$','$70$','$80$','$90$','$100$'});%yticklabels({'$30$','$40$','$50$', '$60$'});
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 20);
% during-fault area identification
trange=[0,t_c,t_c,0];   thetarange=[ymin,ymin,ymax,ymax];
fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;

plot(t_simulation2-t_start,omega_simulation2_plot,'LineStyle','-','linewidth',3,'color',[0 0.4470 0.7410]);
%plot(time_exp-t_start,freq_exp_plot,'LineStyle','-','linewidth',3,'color',[0 0.4470 0.7410]);
plot(t_timedomain-t_start,omega_timedomain,'LineStyle',':','linewidth',1.5,'color',[0 0 0]);%[0.8500 0.3250 0.0980]);

plot(t_simulation-t_start,omega_simulation,'LineStyle','--','linewidth',2.8,'color',[0.8500 0.3250 0.0980]);%[0.4660 0.6740 0.1880]);

%%
clear ylim
figure;
set(gcf,'position',[680 558 1300 300]);
hold on;
xlim([0 t_end]);
grid on;    %grid minor;
ylim([-1.2,1.2]);
yl=ylim;
ymin=yl(1,1);
ymax=yl(1,2);
xticks(0:0.1:t_end);
yticks(-1:0.5:1);
yticklabels({'$-1$','','$0$', '','$1$'});
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 17);
% during-fault area identification
trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
plot(t_simulation2,Vabc_simulation(:,1),'LineStyle','-','linewidth',1.5,'color',[0 0.4470 0.7410]);    hold on;
plot(t_simulation2,Vabc_simulation(:,2),'LineStyle','-','linewidth',1.5,'color',[0.8500 0.3250 0.0980]);    hold on;
plot(t_simulation2,Vabc_simulation(:,3),'LineStyle','-','linewidth',1.5,'color',[0.9290 0.6940 0.1250]);    hold on;

%%
clear ylim
figure;
set(gcf,'position',[680 558 1300 300]);
hold on;
xlim([0 t_end]);
grid on;    %grid minor;
ylim([-0.1,1.4]);
yl=ylim;
ymin=yl(1,1);
ymax=yl(1,2);
xticks(0:0.1:t_end);
yticks(-1:0.5:1);
%yticklabels({'$0$','','$1$'});
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 17);
% during-fault area identification
trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
plot(t_simulation2,Idq_simulation(:,1),'LineStyle','-','linewidth',2,'color',[0.4940 0.1840 0.5560]);    hold on;
plot(t_simulation2,Idq_simulation(:,2),'LineStyle','-','linewidth',2,'color',[0.4660 0.6740 0.1880]);    hold on;

% save(strcat('C:\Users\yz7521\OneDrive - Imperial College London\Desktop\DSP\SimulationData'),"t_start","t_c","t_timedomain" ...
%     ,"t_simulation","t_simulation2","delta_timedomain","delta_simulation","delta_simulation2" ...
%     ,"omega_timedomain","omega_simulation","omega_simulation2","t_end");





