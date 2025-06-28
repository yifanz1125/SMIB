%%
VSC.LCL.Lf = Lf/Wbase;
VSC.LCL.rlf = 1e-8;
VSC.LCL.Cf = Cf2/Wbase;
VSC.LCL.rcf = 1e-8;
VSC.Net.Lg = Lg;
VSC.Net.rg = 1e-3;

f_res1 = sqrt(1/VSC.Net.Lg/VSC.LCL.Cf)/2/pi
f_res2 = sqrt((VSC.Net.Lg+VSC.LCL.Lf)/VSC.LCL.Cf/VSC.Net.Lg/VSC.LCL.Lf)/2/pi

f_pos=logspace(1,4,1e6); %Hz
f_neg=-flip(f_pos);
f_tt=[f_neg,f_pos];
s=1i*f_tt*2*pi;  % omega
n_tt=size(f_tt,2);
fbd_L = min(f_pos);
fbd_H = max(f_pos);
VSC.Ctrl.Ts = Tc;
%% Cf
Tf.LCL.Zcf=1./(s+1i*Wbase)/VSC.LCL.Cf+VSC.LCL.rcf;
Tf.LCL.Zlg=(s+1i*Wbase)*VSC.Net.Lg+VSC.Net.rg;
Tf.LCL.Zlf=(s+1i*Wbase)*VSC.LCL.Lf+VSC.LCL.rlf;
Tf.LCL.Zpara=(Tf.LCL.Zcf.*Tf.LCL.Zlg)./(Tf.LCL.Zcf+Tf.LCL.Zlg);
Tf.LCL.Zseries=Tf.LCL.Zpara+Tf.LCL.Zlf;
Tf.LCL.Yseries=1./Tf.LCL.Zseries;
[Tf.LCL.Yseries_mag,Tf.LCL.Yseries_ang]=Fcn_Cal_BodeMagAng(Tf.LCL.Yseries);


VSC.Ctrl.CCL.kpi = w_i_GFL1*Lf/Wbase;
VSC.Ctrl.CCL.kii = w_i_GFL1*w_i_GFL1/4*Lf/Wbase;
Tf.CCL.PIc=VSC.Ctrl.CCL.kpi+VSC.Ctrl.CCL.kii./s;
Tf.PWM.Gdel=exp(-s*1.5*VSC.Ctrl.Ts);
Tf.CCL.Gcol=Tf.PWM.Gdel.*Tf.CCL.PIc.*Tf.LCL.Yseries;   
[Tf.CCL.Gcol_mag,Tf.CCL.Gcol_ang]=Fcn_Cal_BodeMagAng(Tf.CCL.Gcol);

Tf.LCL.vfilter = Tf.LCL.Zpara./Tf.LCL.Zseries;
Tf.LCL.ifilter = Tf.LCL.Zcf./(Tf.LCL.Zcf+Tf.LCL.Zlg);
[Tf.LCL.vfilter_mag,Tf.LCL.vfilter_ang]=Fcn_Cal_BodeMagAng(Tf.LCL.vfilter);
[Tf.LCL.ifilter_mag,Tf.LCL.ifilter_ang]=Fcn_Cal_BodeMagAng(Tf.LCL.ifilter);
%% paralle inductance
VSC.LCL.r1cf = 2;
VSC.LCL.l1cf = 0.05/Wbase;
Tf.LCL.Zcf1=1./(s+1i*Wbase)/VSC.LCL.Cf+VSC.LCL.rcf+VSC.LCL.r1cf.*VSC.LCL.l1cf.*(s+1i*Wbase)./(VSC.LCL.r1cf+VSC.LCL.l1cf.*(s+1i*Wbase));
Tf.LCL.Zpara1=(Tf.LCL.Zcf1.*Tf.LCL.Zlg)./(Tf.LCL.Zcf1+Tf.LCL.Zlg);
Tf.LCL.Zseries1=Tf.LCL.Zpara1+Tf.LCL.Zlf;
Tf.LCL.Yseries1=1./Tf.LCL.Zseries1;
[Tf.LCL.Yseries_mag1,Tf.LCL.Yseries_ang1]=Fcn_Cal_BodeMagAng(Tf.LCL.Yseries1);
Tf.CCL.Gcol1=Tf.PWM.Gdel.*Tf.CCL.PIc.*Tf.LCL.Yseries1;   
[Tf.CCL.Gcol_mag1,Tf.CCL.Gcol_ang1]=Fcn_Cal_BodeMagAng(Tf.CCL.Gcol1);

Tf.LCL.vfilter1 = Tf.LCL.Zpara1./Tf.LCL.Zseries1;
Tf.LCL.ifilter1 = Tf.LCL.Zcf1./(Tf.LCL.Zcf1+Tf.LCL.Zlg);
[Tf.LCL.vfilter_mag1,Tf.LCL.vfilter_ang1]=Fcn_Cal_BodeMagAng(Tf.LCL.vfilter1);
[Tf.LCL.ifilter_mag1,Tf.LCL.ifilter_ang1]=Fcn_Cal_BodeMagAng(Tf.LCL.ifilter1);

%% paralle resistor
VSC.LCL.rpcf = 10;
Tf.LCL.Zcftem = 1./(s+1i*Wbase)/VSC.LCL.Cf+VSC.LCL.rcf;
Tf.LCL.Zcf2 = Tf.LCL.Zcftem.*VSC.LCL.rpcf./(VSC.LCL.rpcf+Tf.LCL.Zcftem);
Tf.LCL.Zpara2=(Tf.LCL.Zcf2.*Tf.LCL.Zlg)./(Tf.LCL.Zcf2+Tf.LCL.Zlg);
Tf.LCL.Zseries2=Tf.LCL.Zpara2+Tf.LCL.Zlf;
Tf.LCL.Yseries2=1./Tf.LCL.Zseries2;
[Tf.LCL.Yseries_mag2,Tf.LCL.Yseries_ang2]=Fcn_Cal_BodeMagAng(Tf.LCL.Yseries2);
Tf.CCL.Gcol2=Tf.PWM.Gdel.*Tf.CCL.PIc.*Tf.LCL.Yseries2;   
[Tf.CCL.Gcol_mag2,Tf.CCL.Gcol_ang2]=Fcn_Cal_BodeMagAng(Tf.CCL.Gcol2);
Tf.LCL.vfilter = Tf.LCL.Zpara./Tf.LCL.Zseries;
Tf.LCL.ifilter = Tf.LCL.Zcf./Tf.LCL.Zpara;

Tf.LCL.vfilter2 = Tf.LCL.Zpara2./Tf.LCL.Zseries2;
Tf.LCL.ifilter2 = Tf.LCL.Zcf2./(Tf.LCL.Zcf2+Tf.LCL.Zlg);
[Tf.LCL.vfilter_mag2,Tf.LCL.vfilter_ang2]=Fcn_Cal_BodeMagAng(Tf.LCL.vfilter2);
[Tf.LCL.ifilter_mag2,Tf.LCL.ifilter_ang2]=Fcn_Cal_BodeMagAng(Tf.LCL.ifilter2);

%% paralle inductance & capacitor
VSC.LCL.r1cf = 2;
VSC.LCL.c1cf = 0;%0.001/Wbase;
VSC.LCL.l1cf = 0.05/Wbase;
Tf.LCL.Zcf3=1./(s+1i*Wbase)/VSC.LCL.Cf+VSC.LCL.rcf+1./(1./VSC.LCL.r1cf + VSC.LCL.c1cf*(s+1i*Wbase));%1./(s+1i*Wbase)/VSC.LCL.Cf+VSC.LCL.rcf + 1./(1./VSC.LCL.r1cf+1./VSC.LCL.l1cf*(s+1i*Wbase)+VSC.LCL.c1cf*(s+1i*Wbase));
Tf.LCL.Zpara3=(Tf.LCL.Zcf3.*Tf.LCL.Zlg)./(Tf.LCL.Zcf3+Tf.LCL.Zlg);
Tf.LCL.Zseries3=Tf.LCL.Zpara3+Tf.LCL.Zlf;
Tf.LCL.Yseries3=1./Tf.LCL.Zseries3;
[Tf.LCL.Yseries_mag3,Tf.LCL.Yseries_ang3]=Fcn_Cal_BodeMagAng(Tf.LCL.Yseries3);
Tf.CCL.Gcol3=Tf.PWM.Gdel.*Tf.CCL.PIc.*Tf.LCL.Yseries3;   
[Tf.CCL.Gcol_mag3,Tf.CCL.Gcol_ang3]=Fcn_Cal_BodeMagAng(Tf.CCL.Gcol3);

Tf.LCL.vfilter3 = Tf.LCL.Zpara3./Tf.LCL.Zseries3;
Tf.LCL.ifilter3 = Tf.LCL.Zcf3./(Tf.LCL.Zcf3+Tf.LCL.Zlg);
[Tf.LCL.vfilter_mag3,Tf.LCL.vfilter_ang3]=Fcn_Cal_BodeMagAng(Tf.LCL.vfilter3);
[Tf.LCL.ifilter_mag3,Tf.LCL.ifilter_ang3]=Fcn_Cal_BodeMagAng(Tf.LCL.ifilter3);
%%
figure;
set(gcf,'position',[500 100 1000 500]);
% Positive frequency
subplot(2,2,2)
semilogx(f_pos,Tf.LCL.Yseries_mag(n_tt/2+1:end),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-.'); grid on; hold on;
semilogx(f_pos,Tf.CCL.Gcol_mag(n_tt/2+1:end),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-'); hold on;
semilogx(f_pos,Tf.LCL.Yseries_mag1(n_tt/2+1:end),'linewidth',1.5,'Color',[0.8500 0.3250 0.0980],'LineStyle','-.'); hold on;
semilogx(f_pos,Tf.CCL.Gcol_mag1(n_tt/2+1:end),'linewidth',1.5,'Color',[0.8500 0.3250 0.0980],'LineStyle','-'); hold on;
semilogx(f_pos,Tf.LCL.Yseries_mag2(n_tt/2+1:end),'linewidth',1.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','-.'); hold on;
semilogx(f_pos,Tf.CCL.Gcol_mag2(n_tt/2+1:end),'linewidth',1.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','-'); hold on;

semilogx(f_pos,Tf.LCL.Yseries_mag3(n_tt/2+1:end),'linewidth',1.5,'Color',[0.4940 0.1840 0.5560],'LineStyle','-.'); hold on;
semilogx(f_pos,Tf.CCL.Gcol_mag3(n_tt/2+1:end),'linewidth',1.5,'Color',[0.4940 0.1840 0.5560],'LineStyle','-'); hold on;
set(gca,'XLim',[fbd_L fbd_H]);

subplot(2,2,4)
semilogx(f_pos,Tf.LCL.Yseries_ang(n_tt/2+1:end),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-.'); grid on;hold on;
semilogx(f_pos,Tf.CCL.Gcol_ang(n_tt/2+1:end),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-'); hold on;
semilogx(f_pos,Tf.LCL.Yseries_ang1(n_tt/2+1:end),'linewidth',1.5,'Color',[0.8500 0.3250 0.0980],'LineStyle','-.'); hold on;
semilogx(f_pos,Tf.CCL.Gcol_ang1(n_tt/2+1:end),'linewidth',1.5,'Color',[0.8500 0.3250 0.0980],'LineStyle','-'); hold on;
semilogx(f_pos,Tf.LCL.Yseries_ang2(n_tt/2+1:end),'linewidth',1.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','-.'); hold on;
semilogx(f_pos,Tf.CCL.Gcol_ang2(n_tt/2+1:end),'linewidth',1.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','-'); hold on;

semilogx(f_pos,Tf.LCL.Yseries_ang3(n_tt/2+1:end),'linewidth',1.5,'Color',[0.4940 0.1840 0.5560],'LineStyle','-.'); hold on;
semilogx(f_pos,Tf.CCL.Gcol_ang3(n_tt/2+1:end),'linewidth',1.5,'Color',[0.4940 0.1840 0.5560],'LineStyle','-'); hold on;
set(gca,'YLim',[-180 180]);
set(gca,'XLim',[fbd_L fbd_H]);
xlabel('Positive Frequency (Hz)','interpreter','latex','FontSize',12)

% Negative frequency
subplot(2,2,1)
semilogx(f_neg,Tf.LCL.Yseries_mag(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-.'); grid on;hold on;
semilogx(f_neg,Tf.CCL.Gcol_mag(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-');hold on;
semilogx(f_neg,Tf.LCL.Yseries_mag1(1:n_tt/2),'linewidth',1.5,'Color',[0.8500 0.3250 0.0980],'LineStyle','-.'); hold on;
semilogx(f_neg,Tf.CCL.Gcol_mag1(1:n_tt/2),'linewidth',1.5,'Color',[0.8500 0.3250 0.0980],'LineStyle','-');hold on;
semilogx(f_neg,Tf.LCL.Yseries_mag2(1:n_tt/2),'linewidth',1.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','-.'); hold on;
semilogx(f_neg,Tf.CCL.Gcol_mag2(1:n_tt/2),'linewidth',1.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','-');hold on;

semilogx(f_neg,Tf.LCL.Yseries_mag3(1:n_tt/2),'linewidth',1.5,'Color',[0.4940 0.1840 0.5560],'LineStyle','-.'); hold on;
semilogx(f_neg,Tf.CCL.Gcol_mag3(1:n_tt/2),'linewidth',1.5,'Color',[0.4940 0.1840 0.5560],'LineStyle','-');hold on;
set(gca,'XLim',[-fbd_H -fbd_L]);
ylabel('Magnitude (dB)','interpreter','latex','FontSize',12)

subplot(2,2,3)
semilogx(f_neg,Tf.LCL.Yseries_ang(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-.'); grid on;hold on;
semilogx(f_neg,Tf.CCL.Gcol_ang(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-');hold on;
semilogx(f_neg,Tf.LCL.Yseries_ang1(1:n_tt/2),'linewidth',1.5,'Color',[0.8500 0.3250 0.0980],'LineStyle','-'); hold on;
semilogx(f_neg,Tf.CCL.Gcol_ang1(1:n_tt/2),'linewidth',1.5,'Color',[0.8500 0.3250 0.0980],'LineStyle','-');hold on;
semilogx(f_neg,Tf.LCL.Yseries_ang2(1:n_tt/2),'linewidth',1.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','-.'); hold on;
semilogx(f_neg,Tf.CCL.Gcol_ang2(1:n_tt/2),'linewidth',1.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','-');hold on;

semilogx(f_neg,Tf.LCL.Yseries_ang3(1:n_tt/2),'linewidth',1.5,'Color',[0.4940 0.1840 0.5560],'LineStyle','-.'); hold on;
semilogx(f_neg,Tf.CCL.Gcol_ang3(1:n_tt/2),'linewidth',1.5,'Color',[0.4940 0.1840 0.5560],'LineStyle','-');hold on;
set(gca,'XLim',[-fbd_H -fbd_L]);
set(gca,'YLim',[-180 180]);
ylabel('Phase (degree)','interpreter','latex','FontSize',12)
xlabel('Negative Frequency (Hz)','interpreter','latex','FontSize',12)
%%
figure;
set(gcf,'position',[500 100 1000 500]);
% Positive frequency
subplot(2,2,2)
semilogx(f_pos,Tf.LCL.vfilter_mag(n_tt/2+1:end),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-'); grid on; hold on;
semilogx(f_pos,Tf.LCL.ifilter_mag(n_tt/2+1:end),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle',':'); hold on;
semilogx(f_pos,Tf.LCL.vfilter_mag1(n_tt/2+1:end),'linewidth',1.5,'Color',[0.8500 0.3250 0.0980],'LineStyle','-'); hold on;
semilogx(f_pos,Tf.LCL.ifilter_mag1(n_tt/2+1:end),'linewidth',1.5,'Color',[0.8500 0.3250 0.0980],'LineStyle',':'); hold on;
semilogx(f_pos,Tf.LCL.vfilter_mag2(n_tt/2+1:end),'linewidth',1.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','-'); hold on;
semilogx(f_pos,Tf.LCL.ifilter_mag2(n_tt/2+1:end),'linewidth',1.5,'Color',[0.9290 0.6940 0.1250],'LineStyle',':'); hold on;

semilogx(f_pos,Tf.LCL.vfilter_mag3(n_tt/2+1:end),'linewidth',1.5,'Color',[0.4940 0.1840 0.5560],'LineStyle','-'); hold on;
semilogx(f_pos,Tf.LCL.ifilter_mag3(n_tt/2+1:end),'linewidth',1.5,'Color',[0.4940 0.1840 0.5560],'LineStyle',':'); hold on;
set(gca,'XLim',[fbd_L fbd_H]);

subplot(2,2,4)
semilogx(f_pos,Tf.LCL.vfilter_ang(n_tt/2+1:end),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-'); grid on;hold on;
semilogx(f_pos,Tf.LCL.ifilter_ang(n_tt/2+1:end),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle',':'); hold on;
semilogx(f_pos,Tf.LCL.vfilter_ang1(n_tt/2+1:end),'linewidth',1.5,'Color',[0.8500 0.3250 0.0980],'LineStyle','-'); hold on;
semilogx(f_pos,Tf.LCL.ifilter_ang1(n_tt/2+1:end),'linewidth',1.5,'Color',[0.8500 0.3250 0.0980],'LineStyle',':'); hold on;
semilogx(f_pos,Tf.LCL.vfilter_ang2(n_tt/2+1:end),'linewidth',1.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','-'); hold on;
semilogx(f_pos,Tf.LCL.ifilter_ang2(n_tt/2+1:end),'linewidth',1.5,'Color',[0.9290 0.6940 0.1250],'LineStyle',':'); hold on;

semilogx(f_pos,Tf.LCL.vfilter_ang3(n_tt/2+1:end),'linewidth',1.5,'Color',[0.4940 0.1840 0.5560],'LineStyle','-'); hold on;
semilogx(f_pos,Tf.LCL.ifilter_ang3(n_tt/2+1:end),'linewidth',1.5,'Color',[0.4940 0.1840 0.5560],'LineStyle',':'); hold on;
set(gca,'YLim',[-180 180]);
set(gca,'XLim',[fbd_L fbd_H]);
xlabel('Positive Frequency (Hz)','interpreter','latex','FontSize',12)

% Negative frequency
subplot(2,2,1)
semilogx(f_neg,Tf.LCL.vfilter_mag(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-'); grid on;hold on;
semilogx(f_neg,Tf.LCL.ifilter_mag(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle',':');hold on;
semilogx(f_neg,Tf.LCL.vfilter_mag1(1:n_tt/2),'linewidth',1.5,'Color',[0.8500 0.3250 0.0980],'LineStyle','-'); hold on;
semilogx(f_neg,Tf.LCL.ifilter_mag1(1:n_tt/2),'linewidth',1.5,'Color',[0.8500 0.3250 0.0980],'LineStyle',':');hold on;
semilogx(f_neg,Tf.LCL.vfilter_mag2(1:n_tt/2),'linewidth',1.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','-'); hold on;
semilogx(f_neg,Tf.LCL.ifilter_mag2(1:n_tt/2),'linewidth',1.5,'Color',[0.9290 0.6940 0.1250],'LineStyle',':');hold on;

semilogx(f_neg,Tf.LCL.vfilter_mag3(1:n_tt/2),'linewidth',1.5,'Color',[0.4940 0.1840 0.5560],'LineStyle','-'); hold on;
semilogx(f_neg,Tf.LCL.ifilter_mag3(1:n_tt/2),'linewidth',1.5,'Color',[0.4940 0.1840 0.5560],'LineStyle',':');hold on;
set(gca,'XLim',[-fbd_H -fbd_L]);
ylabel('Magnitude (dB)','interpreter','latex','FontSize',12)

subplot(2,2,3)
semilogx(f_neg,Tf.LCL.vfilter_ang(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-'); grid on;hold on;
semilogx(f_neg,Tf.LCL.ifilter_ang(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle',':');hold on;
semilogx(f_neg,Tf.LCL.vfilter_ang1(1:n_tt/2),'linewidth',1.5,'Color',[0.8500 0.3250 0.0980],'LineStyle','-'); hold on;
semilogx(f_neg,Tf.LCL.ifilter_ang1(1:n_tt/2),'linewidth',1.5,'Color',[0.8500 0.3250 0.0980],'LineStyle',':');hold on;
semilogx(f_neg,Tf.LCL.vfilter_ang2(1:n_tt/2),'linewidth',1.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','-'); hold on;
semilogx(f_neg,Tf.LCL.ifilter_ang2(1:n_tt/2),'linewidth',1.5,'Color',[0.9290 0.6940 0.1250],'LineStyle',':');hold on;

semilogx(f_neg,Tf.LCL.vfilter_ang3(1:n_tt/2),'linewidth',1.5,'Color',[0.4940 0.1840 0.5560],'LineStyle','-'); hold on;
semilogx(f_neg,Tf.LCL.ifilter_ang3(1:n_tt/2),'linewidth',1.5,'Color',[0.4940 0.1840 0.5560],'LineStyle',':');hold on;
set(gca,'XLim',[-fbd_H -fbd_L]);
set(gca,'YLim',[-180 180]);
ylabel('Phase (degree)','interpreter','latex','FontSize',12)
xlabel('Negative Frequency (Hz)','interpreter','latex','FontSize',12)