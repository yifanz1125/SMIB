%%

DeltaVFM=ScopeData.signals(3).values; %degree
OmegaVFM=(ScopeData.signals(4).values-1)*Wbase; %rad
VdcVFM=ScopeData.signals(1).values; %degree
IdGFM=ScopeData.signals(2).values(:,1); %Id
IqGFM=ScopeData.signals(2).values(:,2); %Id


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


%%
figure(f1);
hold on;
plot(delta_simulation,y_simulation,'LineStyle','--','linewidth',2,'color','#A2142F');    hold on;

figure(f2);
xl=xlim;
xmin=xl(1,1);
xmax=xl(1,2);
hold on;
plot(t_simulation,delta_simulation,'LineStyle','--','linewidth',2,'color','#A2142F');    hold on;
xlim([xmin,xmax]);



figure(f3);
xl=xlim;
xmin=xl(1,1);
xmax=xl(1,2);
hold on;
plot(t_simulation,y_simulation,'LineStyle','--','linewidth',2,'color','#A2142F');    hold on;
xlim([xmin,xmax]);


