%% =========================================================
%  GFM virtual impedance current limiter:
%  energy surface + level set
%
%  Required variables in workspace:
%  Rg, Xg, Ug, Vgfm, Pm, J, Ws, D, Ilim, prefault_SEP
%
%  Optional:
%  Kvi          default = 0.2
%  delta_uep_vi  UEP angle from VI equilibrium search
%  f1           existing phase-plane figure handle
%% =========================================================

%% ===================== 基本量 =====================
delta = -2*pi:0.01:2*pi;

if ~exist('Kvi','var')
    Kvi = 0.2;
end

E = Vgfm;

% 原始功率，不限流
Pv = Rg*(E^2 - E*Ug*cos(delta))./(Rg^2+Xg^2) ...
   + Xg*E*Ug*sin(delta)./(Rg^2+Xg^2);

% 原始电流，用来判断是否进入虚拟阻抗限流
Den = E^2 + Ug^2 - 2*E*Ug*cos(delta);
Iv0 = sqrt(Den./(Rg^2+Xg^2));

% 进入限流的边界：Xv = 0, Iv0 = Ilim
arg_vi = (E^2 + Ug^2 - Ilim^2*(Xg^2 + Rg^2))/(2*E*Ug);

if abs(arg_vi) <= 1
    deltac = acos(arg_vi);
else
    deltac = NaN;
end

%% ===================== 虚拟阻抗 Xv 与功率 Pi =====================
Xv = zeros(size(delta));
Pi = zeros(size(delta));

for k = 1:length(delta)
    [Pi(k), Xv(k)] = VI_power_scalar(delta(k), E, Ug, Rg, Xg, Ilim, Kvi);
end

%% ===================== 周期映射 =====================
delta_wrap = mod(delta + pi, 2*pi) - pi;

if ~isnan(deltac)
    idx_in  = abs(delta_wrap) <= deltac;
    idx_out = ~idx_in;
else
    idx_in  = Iv0 < Ilim;
    idx_out = ~idx_in;
end

%% ===================== 电流图 =====================
figure;
hold on;

% 实际加虚拟阻抗后的电流
Iv_vi = sqrt(Den./(Rg^2 + (Xg + Xv).^2));

plot(delta, Iv0,  'r--', 'LineWidth', 1.3);
plot(delta, Iv_vi,'r-',  'LineWidth', 1.8);
plot(delta, Ilim*ones(size(delta)), 'k-', 'LineWidth', 1.2);

if ~isnan(deltac)
    plot(deltac, Ilim, 'k.', 'MarkerSize', 18);
    plot(-deltac, Ilim, 'k.', 'MarkerSize', 18);
end

xlabel('\delta');
ylabel('I');
legend('I without VI','I with VI','I_{lim}');
grid on;
box on;

%% ===================== Xv 图 =====================
figure;
hold on;
plot(delta, Xv, 'LineWidth', 1.8);
xlabel('\delta');
ylabel('X_v');
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
Pi_in  = Pi;  Pi_in(~idx_in)   = NaN;
Pi_out = Pi;  Pi_out(~idx_out) = NaN;

plot(delta, Pi_out, 'r-', 'LineWidth', 2);

% --- Pm ---
plot(delta, Pm*ones(size(delta)), 'k-', 'LineWidth', 1.2);

xlabel('\delta');
ylabel('P');
legend('normal P','normal P extension','VI-limited P','P_m');
grid on;
box on;

%% ===================== 构造 V3 / V3_2 =====================
delta_s = prefault_SEP(1);

% 拼接点仍然取进入限流边界 deltac
if isnan(deltac)
    error('No real deltac exists. Check Ilim or system parameters.');
end

deltacc = deltac;
J_ori = J/Ws;
lamda = 0;

% ------- 内区 V3：正常功率解析积分 -------
VV3 = @(delta_val, omega_val) ...
    J_ori/2*(omega_val*Ws)^2 ...
    - Pm*(delta_val-delta_s) ...
    + Rg*(E^2*(delta_val-delta_s) - E*Ug*(sin(delta_val)-sin(delta_s))) / (Rg^2+Xg^2) ...
    - Xg*E*Ug*(cos(delta_val)-cos(delta_s)) / (Rg^2+Xg^2) ...
    + lamda*D*(omega_val)*(delta_val-delta_s) ...
    + lamda/2*D^2/J/Ws*(delta_val-delta_s)^2;

% ------- VI 限流功率函数句柄（供数值积分） -------
Pi_fun = @(x) arrayfun(@(xx) VI_power_scalar(xx, E, Ug, Rg, Xg, Ilim, Kvi), x);

% ------- 拼接常数 -------
% 要求：V3_2(deltacc,0) = V3(deltacc,0)
V3_dacc_0 = VV3(deltacc, 0);
Cmatch = V3_dacc_0 ...
       + Pm*(deltacc - delta_s) ...
       - lamda/2*D^2/J/Ws*(deltacc-delta_s)^2;

% ------- 外区 V3_2：对 VI 限流功率数值积分 -------
VV3_2 = @(delta_val, omega_val) ...
    J_ori/2*(omega_val*Ws)^2 ...
    - Pm*(delta_val-delta_s) ...
    + integral(Pi_fun, deltacc, delta_val, 'ArrayValued', true) ...
    + lamda*D*(omega_val)*(delta_val-delta_s) ...
    + lamda/2*D^2/J/Ws*(delta_val-delta_s)^2 ...
    + Cmatch;

%% =================================
dx = 0.01*pi;
x1 = unique(sort([-deltacc:dx:pi, deltacc]));
x2 = -0.2:0.001:0.2;

[y1, y2] = meshgrid(x1, x2);
zz_total = NaN(size(y1));

tol = 1e-10;

for a = 1:length(x1)
    for b = 1:length(x2)
        delta_now = y1(b,a);
        y_now     = y2(b,a);

        if delta_now <= deltacc + tol
            zz_total(b,a) = VV3(delta_now, y_now);
        else
            zz_total(b,a) = VV3_2(delta_now, y_now);
        end
    end
end

figure;
surf(y1, y2, zz_total, ...
    'EdgeColor', 'none', ...
    'FaceColor', 'interp', ...
    'FaceAlpha', 0.82);
hold on;

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

% --- SEP ---
if delta_s <= deltacc
    z_sep = VV3(delta_s, y_sep);
elseif delta_s <= pi
    z_sep = VV3_2(delta_s, y_sep);
else
    z_sep = NaN;
end

% --- VI 模型 UEP ---
if ~exist('delta_uep_vi','var') || isempty(delta_uep_vi)
    warning('delta_uep_vi does not exist. UEP marker and level set will be skipped.');
    delta_uep_vi = NaN;
end

if delta_uep_vi <= deltacc
    z_uep_vi = VV3(delta_uep_vi, y_uep);
elseif delta_uep_vi <= pi
    z_uep_vi = VV3_2(delta_uep_vi, y_uep);
else
    z_uep_vi = NaN;
end

% ===================== 标出 SEP / UEP =====================
plot3(delta_s, y_sep, z_sep, 'ko', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'g');

text(delta_s, y_sep, z_sep, '  SEP', ...
    'FontSize', 12, ...
    'Color', 'k', ...
    'FontWeight', 'bold');

if ~isnan(delta_uep_vi) && ~isnan(z_uep_vi)
    plot3(delta_uep_vi, y_uep, z_uep_vi, 'ko', ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', 'r');

    text(delta_uep_vi, y_uep, z_uep_vi, '  UEP', ...
        'FontSize', 12, ...
        'Color', 'k', ...
        'FontWeight', 'bold');
end

%% ===================== 经过 [delta_uep_vi, 0] 的能量切面 =====================
y_cut = 0;

if ~isnan(delta_uep_vi)
    if delta_uep_vi <= deltacc
        Vcut = VV3(delta_uep_vi, y_cut);
    elseif delta_uep_vi <= pi
        Vcut = VV3_2(delta_uep_vi, y_cut);
    else
        Vcut = NaN;
    end
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

    tol = 1e-10;

    for a = 1:length(x1_ls)
        for b = 1:length(x2_ls)
            delta_now = DL(b,a);
            y_now     = YL(b,a);

            if delta_now <= deltacc + tol
                Vlevel(b,a) = VV3(delta_now, y_now);
            else
                Vlevel(b,a) = VV3_2(delta_now, y_now);
            end
        end
    end

    contour(DL, YL, Vlevel, [Vcut Vcut], ...
        'Color', [1 0.25 0.25], ...
        'LineWidth', 2.2);
end

%% ===================== local functions =====================
function [P, Xv] = VI_power_scalar(delta,E,Ug,Rg,Xg,Ilim,Kvi)

    Den = E^2 + Ug^2 - 2*E*Ug*cos(delta);
    Den = max(Den,1e-12);

    I0 = sqrt(Den/(Rg^2 + Xg^2));

    if I0 < Ilim
        Xv = 0;
    else
        Xv = solve_Xv_vi(Den,Rg,Xg,Ilim,Kvi);
    end

    Xt = Xg + Xv;

    P = Rg/(Rg^2 + Xt^2) * (E^2 - E*Ug*cos(delta)) ...
      + Xt/(Rg^2 + Xt^2) * E*Ug*sin(delta);
end

function Xv = solve_Xv_vi(Den,Rg,Xg,Ilim,Kvi)

    fun = @(z) z - Kvi*(sqrt(Den./(Rg^2 + (Xg+z).^2)) - Ilim);

    lo = 0;
    hi = max(Kvi*(sqrt(Den/(Rg^2 + Xg^2)) - Ilim),1e-6);

    while fun(hi) < 0
        hi = 2*hi;
        if hi > 1e4
            error('solve_Xv_vi:NoBracket','Cannot bracket Xv solution.');
        end
    end

    Xv = fzero(fun,[lo hi]);
    Xv = max(real(Xv),0);
end
