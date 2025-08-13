%%

DeltaGFM=ScopeData.signals(3).values; %degree
OmegaGFM=(ScopeData.signals(4).values-1)*Wbase; %rad
IdGFM=ScopeData.signals(2).values(:,1); %Id
IqGFM=ScopeData.signals(2).values(:,2); %Id


t_GFM = ScopeData.time;


T_deta=Ts*10;
t_start2 = t_sim_start+t_start;
t_end1 = t_start2+t_c;
t_end2 = t_sim_start+t_end;


delta_simulation = DeltaGFM(t_sim_start/T_deta+1:t_end2/T_deta+1)/180*pi;
omega_simulation = OmegaGFM(t_sim_start/T_deta+1:t_end2/T_deta+1);
t_simulation = t_GFM(t_sim_start/T_deta+1:t_end2/T_deta+1) -t_sim_start;
Iq_simulation = IqGFM(t_sim_start/T_deta+1:t_end2/T_deta+1);


%%
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
plot(t_simulation,omega_simulation,'LineStyle','--','linewidth',2,'color','#A2142F');    hold on;
xlim([xmin,xmax]);


