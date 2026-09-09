%% =========================================================
%  GFM + VA + QV droop energy surface and level set
%
%  Required variables in workspace:
%  Rg, Xg, Ug, Vgfm, Pm, J, Ws, D, Ilim, kq, Qref1,
%  prefault_SEP, delta_uep_va_qv
%
%  Required external function:
%  solve_Uf_qv(delta,Ug,Rg,Xg,UN,kq,Qref)
%
%  Note:
%  Because Uf = Uf(delta), the potential part is evaluated fully
%  by numerical integration. No analytical integral is used.
%% =========================================================

%% ===================== 基本量 =====================
delta = -2*pi:0.01:2*pi;

UN = Vgfm;

if ~exist('kq','var')
    kq = 0.1;
end

if ~exist('Qref1','var')
    fsep = @(dd) Rg*(Vgfm^2 - Vgfm*Ug*cos(dd))/(Rg^2 + Xg^2) ...
               + Xg*Vgfm*Ug*sin(dd)/(Rg^2 + Xg^2) ...
               - Pm;

    deltas = fsolve(fsep,0);

    Qref1 = Xg*(Vgfm^2 - Vgfm*Ug*cos(deltas))/(Rg^2+Xg^2) ...
          - Rg*Vgfm*Ug*sin(deltas)/(Rg^2+Xg^2);
end

Uf_vec = zeros(size(delta));
Pv = zeros(size(delta));
Pi = zeros(size(delta));
Iv = zeros(size(delta));
Xva = zeros(size(delta));

for kk = 1:length(delta)
    dd = delta(kk);

    Uf = solve_Uf_qv(dd,Ug,Rg,Xg,UN,kq,Qref1);
    Uf_vec(kk) = Uf;

    Den = Uf^2 + Ug^2 - 2*Uf*Ug*cos(dd);
    Den = max(Den,1e-12);

    Iv(kk) = sqrt(Den/(Rg^2+Xg^2));

    % no-limit power with QV voltage Uf
    Pv(kk) = Rg*(Uf^2 - Uf*Ug*cos(dd))/(Rg^2+Xg^2) ...
           + Xg*Uf*Ug*sin(dd)/(Rg^2+Xg^2);

    % VA-limited equivalent total reactance
    if Iv(kk) <= Ilim
        Pi(kk) = Pv(kk);
        Xva(kk) = Xg;
    else
        Xvar = sqrt(max(Den/Ilim^2 - Rg^2,0));
        Xva(kk) = Xvar;

        Pi(kk) = Rg/Den*Ilim^2*(Uf^2 - Uf*Ug*cos(dd)) ...
               + Xvar/Den*Ilim^2*Uf*Ug*sin(dd);
    end
end

%% ===================== 限流边界：数值求 I(delta)=Ilim =====================
delta_scan = linspace(0,pi,3000);
I_scan = zeros(size(delta_scan));

for kk = 1:length(delta_scan)
    dd = delta_scan(kk);
    Uf = solve_Uf_qv(dd,Ug,Rg,Xg,UN,kq,Qref1);
    I_scan(kk) = sqrt((Uf^2 + Ug^2 - 2*Uf*Ug*cos(dd))/(Rg^2+Xg^2));
end

idx = find(I_scan >= Ilim,1);

if isempty(idx)
    deltacc = pi;
    warning('No current-limit boundary found in [0,pi]. Use deltacc = pi.');
elseif idx == 1
    deltacc = delta_scan(idx);
else
    fun_I = @(dd) current_qv(dd,Ug,Rg,Xg,UN,kq,Qref1) - Ilim;
    deltacc = fzero(fun_I,[delta_scan(idx-1),delta_scan(idx)]);
end

%% ===================== 周期映射 =====================
delta_wrap = mod(delta + pi, 2*pi) - pi;

idx_in  = abs(delta_wrap) <= deltacc;
idx_out = ~idx_in;

%% ===================== 电流图 =====================
figure;
hold on;
plot(delta, Iv, 'r-', 'LineWidth', 1.5);
plot(delta, Ilim*ones(size(delta)), 'k-', 'LineWidth', 1.2);
plot(deltacc, Ilim, 'k.', 'MarkerSize', 18);
plot(-deltacc, Ilim, 'k.', 'MarkerSize', 18);

xlabel('\delta');
ylabel('I');
grid on;
box on;

%% ===================== Uf 图 =====================
figure;
hold on;
plot(delta, Uf_vec, 'LineWidth', 1.5);
plot(delta, UN*ones(size(delta)), 'k--', 'LineWidth', 1.0);
xlabel('\delta');
ylabel('U_f');
grid on;
box on;

%% ===================== 功率图（分段） =====================
figure;
hold on;

% --- Pv ---
Pv_in  = Pv;  Pv_in(~idx_in)   = NaN;
Pv_out = Pv;  Pv_out(~idx_out) = NaN;

plot(delta, Pv_in,  'b-',  'LineWidth', 2);
plot(delta, Pv_out, 'b--', 'LineWidth', 1.5);

% --- Pi ---
Pi_out = Pi;  Pi_out(~idx_out) = NaN;

plot(delta, Pi_out, 'r-', 'LineWidth', 2);

% --- Pm ---
plot(delta, Pm*ones(size(delta)), 'k-', 'LineWidth', 1.2);

xlabel('\delta');
ylabel('P');
legend('P no limit','P no-limit extension','P VA+QV limit','P_m');
grid on;
box on;

%% ===================== 构造统一数值能量函数 =====================
delta_s = prefault_SEP(1);

J_ori = J/Ws;
lamda = 1;

% P_fun is the actual piecewise post-fault power under VA+QV
P_fun = @(dd) arrayfun(@(xx) P_va_qv_scalar(xx,Ug,Rg,Xg,UN,Ilim,kq,Qref1), dd);

% Entire potential is numerical:
% V = kinetic - Pm*(delta-delta_s) + integral_{delta_s}^{delta} P_fun(x) dx + damping terms
VV = @(delta_val, omega_val) ...
    J_ori/2*(omega_val*Ws)^2 ...
    - Pm*(delta_val-delta_s) ...
    + integral(P_fun, delta_s, delta_val, 'ArrayValued', true) ...
    + lamda*D*(omega_val)*(delta_val-delta_s) ...
    + lamda/2*D^2/J/Ws*(delta_val-delta_s)^2;

%% =================================
dx = 0.01*pi;
x1 = unique(sort([-deltacc:dx:pi, deltacc]));
x2 = -0.2:0.001:0.2;

[y1, y2] = meshgrid(x1, x2);
zz_total = NaN(size(y1));

for a = 1:length(x1)
    for b = 1:length(x2)
        delta_now = y1(b,a);
        y_now     = y2(b,a);

        zz_total(b,a) = VV(delta_now, y_now);
    end
end

figure;
surf(y1, y2, zz_total, ...
    'EdgeColor', 'none', ...
    'FaceColor', 'interp', ...
    'FaceAlpha', 0.82);
hold on
colormap turbo;
shading interp;
colorbar;

xlim([-deltacc, pi]);
ylim([-0.1, 0.1]);
zlim([-0.2, 2]);
clim([0 2]);

view(3);
grid on;
box on;
xlabel('\delta');
ylabel('\omega');
zlabel('V');

%% ===================== 阴影投影（只给内区） =====================
zmin = -0.2;
shadow = zmin * ones(size(y1));

idx_shadow = ~(y1 <= deltacc);
shadow(idx_shadow) = NaN;

surf(y1, y2, shadow, ...
    'FaceColor', [1 0 0], ...
    'FaceAlpha', 0.15, ...
    'EdgeColor', 'none');

%% ===================== SEP / UEP 坐标 =====================
y_sep = 0;
y_uep = 0;

z_sep = VV(delta_s, y_sep);

if ~exist('delta_uep_va_qv','var') || isempty(delta_uep_va_qv)
    warning('delta_uep_va_qv does not exist. UEP marker and level set will be skipped.');
    delta_uep_va_qv = NaN;
end

if ~isnan(delta_uep_va_qv)
    z_uep_va_qv = VV(delta_uep_va_qv, y_uep);
else
    z_uep_va_qv = NaN;
end

% ===================== 标出 SEP / UEP =====================
plot3(delta_s, y_sep, z_sep, 'ko', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'g');

text(delta_s, y_sep, z_sep, '  SEP', ...
    'FontSize', 12, ...
    'Color', 'k', ...
    'FontWeight', 'bold');

if ~isnan(delta_uep_va_qv) && ~isnan(z_uep_va_qv)
    plot3(delta_uep_va_qv, y_uep, z_uep_va_qv, 'ko', ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', 'r');

    text(delta_uep_va_qv, y_uep, z_uep_va_qv, '  UEP', ...
        'FontSize', 12, ...
        'Color', 'k', ...
        'FontWeight', 'bold');
end

%% ===================== 经过 [delta_uep_va_qv, 0] 的能量切面 =====================
y_cut = 0;

if ~isnan(delta_uep_va_qv)
    Vcut = VV(delta_uep_va_qv, y_cut);
else
    Vcut = NaN;
end

if ~isnan(Vcut)
    x_plane = [min(x1), max(x1)];
    y_plane = [min(x2), max(x2)];
    [Xp, Yp] = meshgrid(x_plane, y_plane);
    Zp = Vcut * ones(size(Xp));

    surf(Xp, Yp, Zp, ...
        'FaceColor', [0.5 0.5 0.5], ...
        'FaceAlpha', 0.28, ...
        'EdgeColor', 'none');

    % 交线
    contour3(y1, y2, zz_total, [Vcut Vcut], ...
        'k-', 'LineWidth', 1.5);
end

%% ===================== 在 f1 里画 level set =====================
if exist('f1','var') && ~isnan(Vcut)
    figure(f1);
    hold on;

    x1_ls = -deltacc:0.01*pi:pi;
    yl_tmp = ylim;
    x2_ls = linspace(yl_tmp(1), yl_tmp(2), 400);

    [DL, YL] = meshgrid(x1_ls, x2_ls);
    Vlevel = NaN(size(DL));

    for a = 1:length(x1_ls)
        for b = 1:length(x2_ls)
            delta_now = DL(b,a);
            y_now     = YL(b,a);

            Vlevel(b,a) = VV(delta_now, y_now);
        end
    end

    contour(DL, YL, Vlevel, [Vcut Vcut], ...
        'Color', [1 0.25 0.25], ...
        'LineWidth', 2.2);
end

%% ===================== local functions =====================
function I = current_qv(delta,Ug,Rg,Xg,UN,kq,Qref1)
    Uf = solve_Uf_qv(delta,Ug,Rg,Xg,UN,kq,Qref1);
    Den = Uf^2 + Ug^2 - 2*Uf*Ug*cos(delta);
    Den = max(Den,1e-12);
    I = sqrt(Den/(Rg^2 + Xg^2));
end

function P = P_va_qv_scalar(delta,Ug,Rg,Xg,UN,Ilim,kq,Qref1)

    Uf = solve_Uf_qv(delta,Ug,Rg,Xg,UN,kq,Qref1);

    Den = Uf^2 + Ug^2 - 2*Uf*Ug*cos(delta);
    Den = max(Den,1e-12);

    I0 = sqrt(Den/(Rg^2 + Xg^2));

    if I0 <= Ilim
        P = Rg*(Uf^2 - Uf*Ug*cos(delta))/(Rg^2+Xg^2) ...
          + Xg*Uf*Ug*sin(delta)/(Rg^2+Xg^2);
    else
        Xvar = sqrt(max(Den/Ilim^2 - Rg^2,0));

        P = Rg/Den*Ilim^2*(Uf^2 - Uf*Ug*cos(delta)) ...
          + Xvar/Den*Ilim^2*Uf*Ug*sin(delta);
    end
end
