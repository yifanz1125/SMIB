%%
load("mydata2.mat");
t_fault = Fault_time.signals(1).values(end);
delta_ori = GFL_Test.signals(9).values;
freq_ori = GFL_Test.signals(5).values;
time_ori = GFL_Test.time;
n_start = find(time_ori==t_fault);
Tss= 1e-4;

t_start = 0.1;
t_end = 0.5;
t_c = 18e-3;

for n = 1:length(delta_ori)-1  % Check the "continuous" property of phase angle
    if (delta_ori(n)-delta_ori(n+1)) > 2*pi*4/5
        delta_ori(n+1:length(delta_ori)) = delta_ori(n+1:length(delta_ori)) + 2*pi;
    elseif (delta_ori(n)-delta_ori(n+1)) <  -2*pi*4/5
        delta_ori(n+1:length(delta_ori)) = delta_ori(n+1:length(delta_ori)) - 2*pi;
    end
end
%delta_ori =delta_ori-2*pi;

delta_exp = delta_ori(((n_start-t_start/Tss)):(n_start-t_start/Tss)+t_end/Tss)*180/pi;
freq_exp = freq_ori(((n_start-t_start/Tss)):(n_start-t_start/Tss)+t_end/Tss);
time_exp = time_ori(((n_start-t_start/Tss)):(n_start-t_start/Tss)+t_end/Tss) - t_fault +t_start;

omega_exp_fault = (freq_ori((n_start):(n_start+t_c/Tss))-50)*2*pi;
delta_exp_fault = delta_ori((n_start):(n_start+t_c/Tss));
omega_exp_post = (freq_ori((n_start+t_c/Tss):(n_start+t_end/Tss))-50)*2*pi;
delta_exp_post = delta_ori((n_start+t_c/Tss):(n_start+t_end/Tss));




%%

figure
hold on;
plot(delta_exp_fault,omega_exp_fault,'LineStyle','-','linewidth',2,'color',[0.9290 0.6940 0.1250]);
plot(delta_exp_post,omega_exp_post,'LineStyle','-','linewidth',2,'color',[0.8500 0.3250 0.0980]);


figure
hold on;
plot(time_exp,delta_exp/180*pi,'LineStyle','-','linewidth',2,'color',[0.6350 0.0780 0.1840]);


figure
hold on;
plot(time_exp,freq_exp,'LineStyle','-','linewidth',2,'color',[0.6350 0.0780 0.1840]);