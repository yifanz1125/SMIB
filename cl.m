%%
VSC.LCL.Lf = Lf/Wbase;
VSC.LCL.rlf = Lf/5;
VSC.LCL.Cf = Cf/Wbase*3;
VSC.LCL.rcf = 1e-2;
VSC.Net.Lg = Xg/Wbase;
VSC.Net.rg = Rg;
Tc = 1/20e3;

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
w_i_GFL1 = 1200*2*pi;

%% Cf
Tf.LCL.Zcf=1./(s+1i*Wbase)/VSC.LCL.Cf+VSC.LCL.rcf;
Tf.LCL.Zlg=(s+1i*Wbase)*VSC.Net.Lg+VSC.Net.rg;
Tf.LCL.Zlf=(s+1i*Wbase)*VSC.LCL.Lf+VSC.LCL.rlf;
Tf.LCL.Zpara=(Tf.LCL.Zcf.*Tf.LCL.Zlg)./(Tf.LCL.Zcf+Tf.LCL.Zlg);
Tf.LCL.Zseries=Tf.LCL.Zpara+Tf.LCL.Zlf;
Tf.LCL.Yseries=1./Tf.LCL.Zseries;
[Tf.LCL.Yseries_mag,Tf.LCL.Yseries_ang]=Fcn_Cal_BodeMagAng(Tf.LCL.Yseries);

VSC.Ctrl.CCL.kpi = w_i_GFL1*Lf/Wbase*1.5;
VSC.Ctrl.CCL.kii = w_i_GFL1*w_i_GFL1/4*Lf/Wbase;
Tf.CCL.PIc=VSC.Ctrl.CCL.kpi+VSC.Ctrl.CCL.kii./s;
Tf.PWM.Gdel=exp(-s*1.5*VSC.Ctrl.Ts);
Tf.CCL.Gcol=Tf.PWM.Gdel.*Tf.CCL.PIc.*Tf.LCL.Yseries;   
[Tf.CCL.Gcol_mag,Tf.CCL.Gcol_ang]=Fcn_Cal_BodeMagAng(Tf.CCL.Gcol);

% Closed-loop transfer function from current reference to converter-side
% inductor current. The 1.5*Tc PWM delay has already been included in Gcol.
Tf.CCL.Gcl=Tf.CCL.Gcol./(1+Tf.CCL.Gcol);
[Tf.CCL.Gcl_mag,Tf.CCL.Gcl_ang]=Fcn_Cal_BodeMagAng(Tf.CCL.Gcl);

Tf.LCL.vfilter = Tf.LCL.Zpara./Tf.LCL.Zseries;
Tf.LCL.ifilter = Tf.LCL.Zcf./(Tf.LCL.Zcf+Tf.LCL.Zlg);
[Tf.LCL.vfilter_mag,Tf.LCL.vfilter_ang]=Fcn_Cal_BodeMagAng(Tf.LCL.vfilter);
[Tf.LCL.ifilter_mag,Tf.LCL.ifilter_ang]=Fcn_Cal_BodeMagAng(Tf.LCL.ifilter);

%% Quasi-stationary EVA voltage-feedback loop
% User settings: -3-dB bandwidths in Hz. For the second-order LPF,
% zeta=1/sqrt(2) gives a Butterworth response, so BW_LPF2_Hz is also its
% -3-dB bandwidth.
VSC.Ctrl.EVA.BW_LPF1_Hz = 30;
VSC.Ctrl.EVA.BW_LPF2_Hz = 30;
VSC.Ctrl.EVA.zeta_LPF2 = 1/sqrt(2);

VSC.Ctrl.EVA.w_LPF1 = 2*pi*VSC.Ctrl.EVA.BW_LPF1_Hz;
VSC.Ctrl.EVA.w_LPF2 = 2*pi*VSC.Ctrl.EVA.BW_LPF2_Hz;

% One first-order LPF is applied to each measured dq-voltage component.
Tf.EVA.G_LPF1 = VSC.Ctrl.EVA.w_LPF1./ ...
    (s+VSC.Ctrl.EVA.w_LPF1);

% Second-order low-pass filter with independently adjustable bandwidth.
Tf.EVA.G_LPF2 = VSC.Ctrl.EVA.w_LPF2^2./ ...
    (s.^2+2*VSC.Ctrl.EVA.zeta_LPF2*VSC.Ctrl.EVA.w_LPF2.*s ...
    +VSC.Ctrl.EVA.w_LPF2^2);

% The algebraic virtual-impedance block is an admittance because it maps
% voltage error to current reference. The sign of capacitor-voltage
% feedback is negative, so the characteristic equation is 1+Gloop=0.
if hypot(Rv0,Xv0)==0
    error('Rv0 and Xv0 cannot both be zero for the algebraic EVA model.');
end
Tf.EVA.Yv = 1/(Rv0+1i*Xv0);

% Current reference -> capacitor voltage:
% closed current loop cascaded with Zcf || Zlg.
Tf.EVA.Giv = Tf.CCL.Gcl.*Tf.LCL.Zpara;
[Tf.EVA.Giv_mag,Tf.EVA.Giv_ang]=Fcn_Cal_BodeMagAng(Tf.EVA.Giv);

% Open-loop gains of the capacitor-voltage feedback loop.
Tf.EVA.Gloop_noLPF = Tf.EVA.Yv.*Tf.EVA.Giv;
Tf.EVA.Gloop_LPF1 = Tf.EVA.G_LPF1.*Tf.EVA.Gloop_noLPF;
Tf.EVA.Gloop_LPF2 = Tf.EVA.G_LPF2.*Tf.EVA.Gloop_noLPF;

[Tf.EVA.Gloop_noLPF_mag,Tf.EVA.Gloop_noLPF_ang]= ...
    Fcn_Cal_BodeMagAng(Tf.EVA.Gloop_noLPF);
[Tf.EVA.Gloop_LPF1_mag,Tf.EVA.Gloop_LPF1_ang]= ...
    Fcn_Cal_BodeMagAng(Tf.EVA.Gloop_LPF1);
[Tf.EVA.Gloop_LPF2_mag,Tf.EVA.Gloop_LPF2_ang]= ...
    Fcn_Cal_BodeMagAng(Tf.EVA.Gloop_LPF2);







%%
figure;
set(gcf,'position',[500 100 1000 500]);
% Positive frequency
subplot(2,2,2)
semilogx(f_pos,Tf.LCL.Yseries_mag(n_tt/2+1:end),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-.'); grid on; hold on;
semilogx(f_pos,Tf.CCL.Gcol_mag(n_tt/2+1:end),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-'); hold on;

set(gca,'XLim',[fbd_L fbd_H]);

subplot(2,2,4)
semilogx(f_pos,Tf.LCL.Yseries_ang(n_tt/2+1:end),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-.'); grid on;hold on;
semilogx(f_pos,Tf.CCL.Gcol_ang(n_tt/2+1:end),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-'); hold on;

set(gca,'YLim',[-180 180]);
set(gca,'XLim',[fbd_L fbd_H]);
xlabel('Positive Frequency (Hz)','interpreter','latex','FontSize',12)

% Negative frequency
subplot(2,2,1)
semilogx(f_neg,Tf.LCL.Yseries_mag(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-.'); grid on;hold on;
semilogx(f_neg,Tf.CCL.Gcol_mag(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-');hold on;
set(gca,'XLim',[-fbd_H -fbd_L]);
ylabel('Magnitude (dB)','interpreter','latex','FontSize',12)

subplot(2,2,3)
semilogx(f_neg,Tf.LCL.Yseries_ang(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-.'); grid on;hold on;
semilogx(f_neg,Tf.CCL.Gcol_ang(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-');hold on;
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

set(gca,'XLim',[fbd_L fbd_H]);

subplot(2,2,4)
semilogx(f_pos,Tf.LCL.vfilter_ang(n_tt/2+1:end),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-'); grid on;hold on;
semilogx(f_pos,Tf.LCL.ifilter_ang(n_tt/2+1:end),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle',':'); hold on;

set(gca,'YLim',[-180 180]);
set(gca,'XLim',[fbd_L fbd_H]);
xlabel('Positive Frequency (Hz)','interpreter','latex','FontSize',12)

% Negative frequency
subplot(2,2,1)
semilogx(f_neg,Tf.LCL.vfilter_mag(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-'); grid on;hold on;
semilogx(f_neg,Tf.LCL.ifilter_mag(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle',':');hold on;

set(gca,'XLim',[-fbd_H -fbd_L]);
ylabel('Magnitude (dB)','interpreter','latex','FontSize',12)

subplot(2,2,3)
semilogx(f_neg,Tf.LCL.vfilter_ang(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-'); grid on;hold on;
semilogx(f_neg,Tf.LCL.ifilter_ang(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle',':');hold on;

set(gca,'XLim',[-fbd_H -fbd_L]);
set(gca,'YLim',[-180 180]);
ylabel('Phase (degree)','interpreter','latex','FontSize',12)
xlabel('Negative Frequency (Hz)','interpreter','latex','FontSize',12)

%% Current closed loop cascaded with Zcf || Zlg: i_ref -> v_c
figure;
set(gcf,'position',[500 100 1000 500]);

% Positive frequency
subplot(2,2,2)
semilogx(f_pos,Tf.EVA.Giv_mag(n_tt/2+1:end), ...
    'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-');
grid on;
set(gca,'XLim',[fbd_L fbd_H]);
title('$G_{i,\mathrm{cl}}(Z_{cf}\parallel Z_{lg})$', ...
    'interpreter','latex','FontSize',12)

subplot(2,2,4)
semilogx(f_pos,Tf.EVA.Giv_ang(n_tt/2+1:end), ...
    'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-');
grid on;
set(gca,'YLim',[-180 180]);
set(gca,'XLim',[fbd_L fbd_H]);
xlabel('Positive Frequency (Hz)','interpreter','latex','FontSize',12)

% Negative frequency
subplot(2,2,1)
semilogx(f_neg,Tf.EVA.Giv_mag(1:n_tt/2), ...
    'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-');
grid on;
set(gca,'XLim',[-fbd_H -fbd_L]);
ylabel('Magnitude (dB)','interpreter','latex','FontSize',12)

subplot(2,2,3)
semilogx(f_neg,Tf.EVA.Giv_ang(1:n_tt/2), ...
    'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-');
grid on;
set(gca,'XLim',[-fbd_H -fbd_L]);
set(gca,'YLim',[-180 180]);
ylabel('Phase (degree)','interpreter','latex','FontSize',12)
xlabel('Negative Frequency (Hz)','interpreter','latex','FontSize',12)

%% EVA capacitor-voltage loop: no LPF, first-order LPF, second-order LPF
figure;
set(gcf,'position',[500 100 1000 500]);

clr_noLPF = [0 0 0];
clr_LPF1 = [0 0.4470 0.7410];
clr_LPF2 = [0.8500 0.3250 0.0980];

% Positive frequency
subplot(2,2,2)
semilogx(f_pos,Tf.EVA.Gloop_noLPF_mag(n_tt/2+1:end), ...
    'linewidth',1.5,'Color',clr_noLPF,'LineStyle','-');
grid on; hold on;
semilogx(f_pos,Tf.EVA.Gloop_LPF1_mag(n_tt/2+1:end), ...
    'linewidth',1.5,'Color',clr_LPF1,'LineStyle','--');
semilogx(f_pos,Tf.EVA.Gloop_LPF2_mag(n_tt/2+1:end), ...
    'linewidth',1.5,'Color',clr_LPF2,'LineStyle','-.');
yline(0,':','Color',[0.5 0.5 0.5],'HandleVisibility','off');
set(gca,'XLim',[fbd_L fbd_H]);
legend('No LPF', ...
    sprintf('First-order LPF: %.1f Hz',VSC.Ctrl.EVA.BW_LPF1_Hz), ...
    sprintf('Second-order LPF: %.1f Hz',VSC.Ctrl.EVA.BW_LPF2_Hz), ...
    'Location','best','Interpreter','latex');
title('$Y_vG_{i,\mathrm{cl}}(Z_{cf}\parallel Z_{lg})G_{\mathrm{LPF}}$', ...
    'interpreter','latex','FontSize',12)

subplot(2,2,4)
semilogx(f_pos,Tf.EVA.Gloop_noLPF_ang(n_tt/2+1:end), ...
    'linewidth',1.5,'Color',clr_noLPF,'LineStyle','-');
grid on; hold on;
semilogx(f_pos,Tf.EVA.Gloop_LPF1_ang(n_tt/2+1:end), ...
    'linewidth',1.5,'Color',clr_LPF1,'LineStyle','--');
semilogx(f_pos,Tf.EVA.Gloop_LPF2_ang(n_tt/2+1:end), ...
    'linewidth',1.5,'Color',clr_LPF2,'LineStyle','-.');
yline(-180,':','Color',[0.5 0.5 0.5],'HandleVisibility','off');
set(gca,'YLim',[-180 180]);
set(gca,'XLim',[fbd_L fbd_H]);
xlabel('Positive Frequency (Hz)','interpreter','latex','FontSize',12)

% Negative frequency
subplot(2,2,1)
semilogx(f_neg,Tf.EVA.Gloop_noLPF_mag(1:n_tt/2), ...
    'linewidth',1.5,'Color',clr_noLPF,'LineStyle','-');
grid on; hold on;
semilogx(f_neg,Tf.EVA.Gloop_LPF1_mag(1:n_tt/2), ...
    'linewidth',1.5,'Color',clr_LPF1,'LineStyle','--');
semilogx(f_neg,Tf.EVA.Gloop_LPF2_mag(1:n_tt/2), ...
    'linewidth',1.5,'Color',clr_LPF2,'LineStyle','-.');
yline(0,':','Color',[0.5 0.5 0.5],'HandleVisibility','off');
set(gca,'XLim',[-fbd_H -fbd_L]);
ylabel('Magnitude (dB)','interpreter','latex','FontSize',12)

subplot(2,2,3)
semilogx(f_neg,Tf.EVA.Gloop_noLPF_ang(1:n_tt/2), ...
    'linewidth',1.5,'Color',clr_noLPF,'LineStyle','-');
grid on; hold on;
semilogx(f_neg,Tf.EVA.Gloop_LPF1_ang(1:n_tt/2), ...
    'linewidth',1.5,'Color',clr_LPF1,'LineStyle','--');
semilogx(f_neg,Tf.EVA.Gloop_LPF2_ang(1:n_tt/2), ...
    'linewidth',1.5,'Color',clr_LPF2,'LineStyle','-.');
yline(-180,':','Color',[0.5 0.5 0.5],'HandleVisibility','off');
set(gca,'XLim',[-fbd_H -fbd_L]);
set(gca,'YLim',[-180 180]);
ylabel('Phase (degree)','interpreter','latex','FontSize',12)
xlabel('Negative Frequency (Hz)','interpreter','latex','FontSize',12)








%%
figure;
set(gcf,'position',[500 100 1000 500]);
% Positive frequency
subplot(2,2,2)
semilogx(f_pos,Tf.LCL.Yseries_mag(n_tt/2+1:end),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-.'); grid on; hold on;
semilogx(f_pos,Tf.CCL.Gcol_mag(n_tt/2+1:end),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-'); hold on;

set(gca,'XLim',[fbd_L fbd_H]);

subplot(2,2,4)
semilogx(f_pos,Tf.LCL.Yseries_ang(n_tt/2+1:end),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-.'); grid on;hold on;
semilogx(f_pos,Tf.CCL.Gcol_ang(n_tt/2+1:end),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-'); hold on;

set(gca,'YLim',[-180 180]);
set(gca,'XLim',[fbd_L fbd_H]);
xlabel('Positive Frequency (Hz)','interpreter','latex','FontSize',12)

% Negative frequency
subplot(2,2,1)
semilogx(f_neg,Tf.LCL.Yseries_mag(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-.'); grid on;hold on;
semilogx(f_neg,Tf.CCL.Gcol_mag(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-');hold on;
set(gca,'XLim',[-fbd_H -fbd_L]);
ylabel('Magnitude (dB)','interpreter','latex','FontSize',12)

subplot(2,2,3)
semilogx(f_neg,Tf.LCL.Yseries_ang(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-.'); grid on;hold on;
semilogx(f_neg,Tf.CCL.Gcol_ang(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-');hold on;
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

set(gca,'XLim',[fbd_L fbd_H]);

subplot(2,2,4)
semilogx(f_pos,Tf.LCL.vfilter_ang(n_tt/2+1:end),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-'); grid on;hold on;
semilogx(f_pos,Tf.LCL.ifilter_ang(n_tt/2+1:end),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle',':'); hold on;

set(gca,'YLim',[-180 180]);
set(gca,'XLim',[fbd_L fbd_H]);
xlabel('Positive Frequency (Hz)','interpreter','latex','FontSize',12)

% Negative frequency
subplot(2,2,1)
semilogx(f_neg,Tf.LCL.vfilter_mag(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-'); grid on;hold on;
semilogx(f_neg,Tf.LCL.ifilter_mag(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle',':');hold on;

set(gca,'XLim',[-fbd_H -fbd_L]);
ylabel('Magnitude (dB)','interpreter','latex','FontSize',12)

subplot(2,2,3)
semilogx(f_neg,Tf.LCL.vfilter_ang(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle','-'); grid on;hold on;
semilogx(f_neg,Tf.LCL.ifilter_ang(1:n_tt/2),'linewidth',1.5,'Color',[0 0.4470 0.7410],'LineStyle',':');hold on;

set(gca,'XLim',[-fbd_H -fbd_L]);
set(gca,'YLim',[-180 180]);
ylabel('Phase (degree)','interpreter','latex','FontSize',12)
xlabel('Negative Frequency (Hz)','interpreter','latex','FontSize',12)

