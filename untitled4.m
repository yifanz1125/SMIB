p11=relay.signals(1).values;
p12=relay.signals(2).values;
p21=relay.signals(3).values;
p22=relay.signals(4).values;
tt=relay.time;


T_deta= Ts*10;
t_start= 1.5;
t_end1 = t_start+tc2;
t_end2 = 2.5;

p11=p11(t_start/T_deta+1:t_end2/T_deta+1);
p12=p12(t_start/T_deta+1:t_end2/T_deta+1);
p21=p21(t_start/T_deta+1:t_end2/T_deta+1);
p22=p22(t_start/T_deta+1:t_end2/T_deta+1);
tt = tt(t_start/T_deta+1:t_end2/T_deta+1);

figure;
subplot(2, 1, 1);
yyaxis left
plot(tt,p11,'-','LineWidth',2);
ylim([0 20]);
yticklabels({'0', '', '10', '','20'});
grid on; 
yyaxis right
plot(tt,p12,'-','LineWidth',2);
ylim([-180 180]);
yticks(-180:90:180);
yticklabels({'-180', '', '0', '','180'});
xticks(1.5:0.1:2.5);
xticklabels([]);
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=12;


subplot(2, 1, 2);
yyaxis left
plot(tt,p21,'-','LineWidth',2);
ylim([0 2]);
yticks(0:0.5:2);
yticklabels({'0', '', '1', '','2'});
yyaxis right
plot(tt,p22,'-','LineWidth',2);
ylim([-180 180]);
yticks(-180:90:180);
yticklabels({'-180', '', '0', '','180'});
grid on; 
xticklabels({'1.5','','','','','2.0','','','','','2.5'});
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=12;
 
