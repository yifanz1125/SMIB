%% GFM circular current limiter + Q-V droop energy surface
% Optimized: precompute P(delta) once, then use cumtrapz/interpolation.

delta = -2*pi:0.01:2*pi;
UN = Vgfm;

Uf_vec = zeros(size(delta));
P_actual = zeros(size(delta));
P_unlimited = zeros(size(delta));
Q_actual = zeros(size(delta));
I_actual = zeros(size(delta));
Re_vec = zeros(size(delta));

for kk = 1:length(delta)
    dd = delta(kk);
    [P_actual(kk),Q_actual(kk),Uf_vec(kk),Re_vec(kk),I_actual(kk),~] = ...
        circle_qv_power(dd,Ug,Rg,Xg,UN,Ilim,kq,Qref1);

    Uf0 = solve_Uf_qv(dd,Ug,Rg,Xg,UN,kq,Qref1);
    P_unlimited(kk) = Rg/(Rg^2+Xg^2)*(Uf0^2-Uf0*Ug*cos(dd)) ...
                    + Xg/(Rg^2+Xg^2)*Uf0*Ug*sin(dd);
end

%% Current-limit boundary
delta_scan = linspace(0,pi,3000);
I_scan = zeros(size(delta_scan));
for kk = 1:length(delta_scan)
    dd = delta_scan(kk);
    Uf0 = solve_Uf_qv(dd,Ug,Rg,Xg,UN,kq,Qref1);
    I_scan(kk) = sqrt((Uf0^2+Ug^2-2*Uf0*Ug*cos(dd))/(Rg^2+Xg^2));
end

idx = find(I_scan>=Ilim,1);
if isempty(idx)
    deltacc = pi;
elseif idx == 1
    deltacc = delta_scan(1);
else
    current_boundary = @(dd) current_unlimited_qv_local(dd,Ug,Rg,Xg,UN,kq,Qref1)-Ilim;
    deltacc = fzero(current_boundary,[delta_scan(idx-1),delta_scan(idx)]);
end

delta_wrap = mod(delta+pi,2*pi)-pi;
idx_in = abs(delta_wrap)<=deltacc;
idx_out = ~idx_in;

%% Current plot
figure; hold on
plot(delta,I_actual,'r-','LineWidth',1.5);
plot(delta,Ilim*ones(size(delta)),'k-','LineWidth',1.2);
plot(deltacc,Ilim,'k.','MarkerSize',18);
plot(-deltacc,Ilim,'k.','MarkerSize',18);
xlabel('\\delta'); ylabel('I'); grid on; box on;

%% Uf plot
figure; hold on
plot(delta,Uf_vec,'LineWidth',1.5);
plot(delta,UN*ones(size(delta)),'k--','LineWidth',1);
xlabel('\\delta'); ylabel('U_f'); grid on; box on;

%% Re plot
figure;
plot(delta,Re_vec,'LineWidth',1.5);
xlabel('\\delta'); ylabel('R_e'); grid on; box on;

%% Piecewise power plot
figure; hold on
P_in = P_unlimited; P_in(~idx_in) = NaN;
P_ext = P_unlimited; P_ext(~idx_out) = NaN;
P_lim = P_actual; P_lim(~idx_out) = NaN;
plot(delta,P_in,'b-','LineWidth',2);
plot(delta,P_ext,'b--','LineWidth',1.5);
plot(delta,P_lim,'r-','LineWidth',2);
plot(delta,Pm*ones(size(delta)),'k-','LineWidth',1.2);
xlabel('\\delta'); ylabel('P'); grid on; box on;

%% Precompute potential energy once
delta_s = prefault_SEP(1);
J_ori = J/Ws;
lamda = 1;

dphi = 0.0025*pi;
delta_energy = unique(sort([-deltacc:dphi:pi,-deltacc,delta_s,deltacc,pi]));
P_energy = zeros(size(delta_energy));

for kk = 1:length(delta_energy)
    [P_energy(kk),~,~,~,~,~] = circle_qv_power(...
        delta_energy(kk),Ug,Rg,Xg,UN,Ilim,kq,Qref1);
end

Phi_raw = cumtrapz(delta_energy,P_energy-Pm);
Phi_at_sep = interp1(delta_energy,Phi_raw,delta_s,'pchip');
Phi_energy = Phi_raw-Phi_at_sep;
Phi_fun = @(dd) interp1(delta_energy,Phi_energy,dd,'pchip','extrap');

VV = @(delta_val,omega_val) ...
    J_ori/2*(omega_val*Ws).^2 ...
    +Phi_fun(delta_val) ...
    +lamda*D*omega_val.*(delta_val-delta_s) ...
    +lamda/2*D^2/J/Ws*(delta_val-delta_s).^2;

%% 3D energy surface
dx = 0.01*pi;
x1 = unique(sort([-deltacc:dx:pi,deltacc,delta_s]));
x2 = -0.2:0.001:0.2;
[y1,y2] = meshgrid(x1,x2);

zz_total = ...
    J_ori/2*(y2*Ws).^2 ...
    +Phi_fun(y1) ...
    +lamda*D.*y2.*(y1-delta_s) ...
    +lamda/2*D^2/J/Ws.*(y1-delta_s).^2;

figure;
surf(y1,y2,zz_total,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.82);
hold on; colormap turbo; shading interp; colorbar;
xlim([-deltacc,pi]); ylim([-0.1,0.1]); zlim([-0.2,2]); clim([0,2]);
view(3); grid on; box on;
xlabel('\\delta'); ylabel('\\omega'); zlabel('V');

%% Shadow
zmin = -0.2;
shadow = zmin*ones(size(y1));
shadow(~(y1<=deltacc)) = NaN;
surf(y1,y2,shadow,'FaceColor',[1 0 0],'FaceAlpha',0.15,'EdgeColor','none');

%% SEP / UEP
z_sep = VV(delta_s,0);
if ~exist('delta_uep_cir_qv','var') || isempty(delta_uep_cir_qv)
    delta_uep_cir_qv = NaN;
end
if ~isnan(delta_uep_cir_qv)
    z_uep = VV(delta_uep_cir_qv,0);
else
    z_uep = NaN;
end

plot3(delta_s,0,z_sep,'ko','MarkerSize',8,'MarkerFaceColor','g');
text(delta_s,0,z_sep,'  SEP','FontSize',12,'FontWeight','bold');

if ~isnan(z_uep)
    plot3(delta_uep_cir_qv,0,z_uep,'ko','MarkerSize',8,'MarkerFaceColor','r');
    text(delta_uep_cir_qv,0,z_uep,'  UEP','FontSize',12,'FontWeight','bold');
    Vcut = z_uep;
    [Xp,Yp] = meshgrid([min(x1),max(x1)],[min(x2),max(x2)]);
    Zp = Vcut*ones(size(Xp));
    surf(Xp,Yp,Zp,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.28,'EdgeColor','none');
    contour3(y1,y2,zz_total,[Vcut Vcut],'k-','LineWidth',1.5);
else
    Vcut = NaN;
end

%% Level set in f1
if exist('f1','var') && ~isnan(Vcut)
    figure(f1); hold on
    x1_ls = unique(sort([-deltacc:0.01*pi:pi,deltacc,delta_s,delta_uep_cir_qv]));
    yl_tmp = ylim;
    x2_ls = linspace(yl_tmp(1),yl_tmp(2),400);
    [DL,YL] = meshgrid(x1_ls,x2_ls);

    Vlevel = ...
        J_ori/2*(YL*Ws).^2 ...
        +Phi_fun(DL) ...
        +lamda*D.*YL.*(DL-delta_s) ...
        +lamda/2*D^2/J/Ws.*(DL-delta_s).^2;

    contour(DL,YL,Vlevel,[Vcut Vcut],'Color',[1 0.25 0.25],'LineWidth',2.2);
end

function I = current_unlimited_qv_local(delta,Ug,Rg,Xg,UN,kq,Qref1)
    Uf = solve_Uf_qv(delta,Ug,Rg,Xg,UN,kq,Qref1);
    Den = max(Uf^2+Ug^2-2*Uf*Ug*cos(delta),1e-12);
    I = sqrt(Den/(Rg^2+Xg^2));
end
