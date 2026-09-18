%% ===================== EVA 基本量 =====================
% Run SMIB_transient_CL first with system = "GFM" and limit_type = "EVA".
% This script uses prefault_SEP, delta_uep_eva and f1 generated there.
% It reproduces all functions of current_limit_va_GFMdamping: current and
% power curves, the piecewise 3-D energy surface, the lambda = 1 practical
% stability boundary, and Moon closed sets for lambda = 0:0.2:1.
%
% Rg+jXg: physical grid impedance.
% Rv0+jXv0: nominal virtual impedance in normal operation.
% During limiting, the virtual-impedance magnitude is enlarged while the
% actual angle defined by Rv0 and Xv0 is retained.
delta = -2*pi:0.01:2*pi;

Zv0_abs = hypot(Rv0, Xv0);
if Zv0_abs <= eps
    error(['EVA requires a nonzero nominal virtual impedance ' ...
           'Rv0+jXv0 to define its limiting direction.']);
end

cos_phi_v = Rv0/Zv0_abs;
sin_phi_v = Xv0/Zv0_abs;
Rsum0 = Rg + Rv0;
Xsum0 = Xg + Xv0;
Zsum20 = Rsum0^2 + Xsum0^2;

% Normal-operation PCC/POC active power.
Pv = (Xsum0*Vgfm*Ug.*sin(delta) ...
    + Rg*(Vgfm^2 - Vgfm*Ug.*cos(delta)) ...
    + Rv0*(Vgfm*Ug.*cos(delta) - Ug^2))./Zsum20;

Den = Vgfm^2 + Ug^2 - 2*Vgfm*Ug*cos(delta);
projection_v = Rg*cos_phi_v + Xg*sin_phi_v;
radicand_v = projection_v^2 + Den/Ilim^2 - (Rg^2 + Xg^2);
lambda_eva = -projection_v + sqrt(max(radicand_v, 0));
lambda_eva = max(lambda_eva, Zv0_abs);

Rv_eva = lambda_eva*cos_phi_v;
Xv_eva = lambda_eva*sin_phi_v;
Rsum_eva = Rg + Rv_eva;
Xsum_eva = Xg + Xv_eva;
Zsum2_eva = Rsum_eva.^2 + Xsum_eva.^2;

% Current-limited PCC/POC active power.
Pi = (Xsum_eva.*Vgfm*Ug.*sin(delta) ...
    + Rg*(Vgfm^2 - Vgfm*Ug.*cos(delta)) ...
    + Rv_eva.*(Vgfm*Ug.*cos(delta) - Ug^2))./Zsum2_eva;

Iv = sqrt(Den/Zsum20);

arg_deltac = (Vgfm^2+Ug^2-Ilim^2*Zsum20)/(2*Vgfm*Ug);
if abs(arg_deltac) > 1
    error('No real EVA current-limiting boundary for the present parameters.');
end
deltac = acos(arg_deltac);

%% ===================== 周期映射 =====================
delta_wrap = mod(delta + pi, 2*pi) - pi;

idx_in  = abs(delta_wrap) <= deltac;
idx_out = ~idx_in;

%% ===================== 电流图 =====================
figure;
hold on;
plot(delta, Iv, 'r-', 'LineWidth', 1.5);
plot(delta, Ilim*ones(size(delta)), 'k-', 'LineWidth', 1.2);
plot(deltac, Ilim, 'k.', 'MarkerSize', 18);
plot(-deltac, Ilim, 'k.', 'MarkerSize', 18);

xlabel('\delta');
ylabel('I');
grid on;
box on;

%% ===================== 功率图（分段） =====================
figure;
hold on;

% --- Pv ---
Pv_in  = Pv;  Pv_in(~idx_in)   = NaN;
Pv_out = Pv;  Pv_out(~idx_out) = NaN;

plot(delta, Pv_in,  'b-',  'LineWidth', 2);
plot(delta, Pv_out, 'k-', 'LineWidth', 1.5);

% --- Pi ---
Pi_in  = Pi;  Pi_in(~idx_in)   = NaN;
Pi_out = Pi;  Pi_out(~idx_out) = NaN;

plot(delta, Pi_out, 'r-', 'LineWidth', 2);

% --- Pm ---
plot(delta, Pm*ones(size(delta)), 'k-', 'LineWidth', 1.2);

axis([-pi 3/2*pi -2 2]);

xlabel('\delta');
ylabel('P');
grid on;
box on;

%% ===================== 构造 V3 / V3_2 =====================
% 平衡点
delta_s   = prefault_SEP(1);

% 这里用右侧边界 deltacc
deltacc = deltac;
J_ori = J/Ws;
lamda = 1;

% ------- 内区 V3 的函数句柄 -------
VV3 = @(delta_val, omega_val) ...
    J_ori/2*(omega_val*Ws)^2 ...
    - Pm*(delta_val-delta_s) ...
    + (Rg*(Vgfm^2*(delta_val-delta_s) ...
          - Vgfm*Ug*(sin(delta_val)-sin(delta_s))) ...
       + Rv0*(Vgfm*Ug*(sin(delta_val)-sin(delta_s)) ...
          - Ug^2*(delta_val-delta_s)) ...
       - Xsum0*Vgfm*Ug*(cos(delta_val)-cos(delta_s))) / Zsum20 ...
    + lamda*D*(omega_val)*(delta_val-delta_s) + lamda/2*D^2/J/Ws*(delta_val-delta_s)^2;

% ------- EVA 限流功率函数句柄（供数值积分） -------
Den_fun = @(x) Vgfm^2 + Ug^2 - 2*Vgfm*Ug.*cos(x);
lambda_fun = @(x) max(Zv0_abs, -projection_v + sqrt(max( ...
    projection_v^2 + Den_fun(x)/Ilim^2 - (Rg^2+Xg^2), 0)));
Rv_fun = @(x) lambda_fun(x).*cos_phi_v;
Xv_fun = @(x) lambda_fun(x).*sin_phi_v;
Zsum2_fun = @(x) (Rg+Rv_fun(x)).^2 + (Xg+Xv_fun(x)).^2;

Pi_fun = @(x) ( ...
      (Xg+Xv_fun(x)).*Vgfm*Ug.*sin(x) ...
    + Rg*(Vgfm^2 - Vgfm*Ug.*cos(x)) ...
    + Rv_fun(x).*(Vgfm*Ug.*cos(x) - Ug^2)) ...
    ./ Zsum2_fun(x);
% ------- 拼接常数 -------
% 要求：V3_2(deltacc,0) = V3(deltacc,0)
% 因为积分上限=下限时积分为0
V3_dacc_0 = VV3(deltacc, 0);
Cmatch = V3_dacc_0 + Pm*(deltacc - delta_s) - lamda/2*D^2/J/Ws*(deltacc-delta_s)^2;

% ------- 外区 V3_2 的函数句柄 -------
VV3_2 = @(delta_val, omega_val) ...
    J_ori/2*(omega_val*Ws)^2 ...
    - Pm*(delta_val-delta_s) ...
    + integral(Pi_fun, deltacc, delta_val, 'ArrayValued', true) ...
    + lamda*D*(omega_val)*(delta_val-delta_s) + lamda/2*D^2/J/Ws*(delta_val-delta_s)^2 + Cmatch;

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
hsurf = surf(y1, y2, zz_total, ...
    'EdgeColor', 'none', ...
    'FaceColor', 'interp', ...
    'FaceAlpha', 0.82);   % 这里调透明度
hold on
colormap turbo;
shading interp;
colorbar;


xlim([-deltacc, pi]);
ylim([-0.1, 0.1]);
zlim([-0.2, 2]);

clim([0 2])

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

% --- EVA 模型 UEP ---
if delta_uep_eva <= deltacc
    z_uep_eva = VV3(delta_uep_eva, y_uep);
elseif delta_uep_eva <= pi
    z_uep_eva = VV3_2(delta_uep_eva, y_uep);
else
    z_uep_eva = NaN;
end

% ===================== 标出 SEP / UEP =====================

plot3(delta_s, y_sep, z_sep, 'ko', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'g');

plot3(delta_uep_eva, y_uep, z_uep_eva, 'ko', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'r');

text(delta_s, y_sep, z_sep, '  SEP', ...
    'FontSize', 12, ...
    'Color', 'k', ...
    'FontWeight', 'bold');

text(delta_uep_eva, y_uep, z_uep_eva, '  UEP', ...
    'FontSize', 12, ...
    'Color', 'k', ...
    'FontWeight', 'bold');


%% ===================== 经过 [delta_uep_eva, 0] 的能量切面 =====================

y_cut = 0;

if delta_uep_eva <= deltacc
    Vcut = VV3(delta_uep_eva, y_cut);
elseif delta_uep_eva <= pi
    Vcut = VV3_2(delta_uep_eva, y_cut);
else
    Vcut = NaN;
end
% 
% %% ===================== 在三维图中画灰色切面 =====================
% 
x_plane = [min(x1), max(x1)];
y_plane = [min(x2), max(x2)];
[Xp, Yp] = meshgrid(x_plane, y_plane);
Zp = Vcut * ones(size(Xp));

surf(Xp, Yp, Zp, ...
    'FaceColor', [0.5 0.5 0.5], ...
    'FaceAlpha', 0.28, ...
    'EdgeColor', 'none');

% 交线（可选）
contour3(y1, y2, zz_total, [Vcut Vcut], ...
    'k-', 'LineWidth', 1.5);
% %% ===================== 在 f1 里画 level set =====================
% 
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

%% ===================== Moon 不同 lambda 对应的闭集（仅作展示） =====================
% Moon 的 S_lambda 由能量函数 E_lambda 自身的鞍点 X_{u,lambda}
% 所对应的等能量线围成，而不是统一经过物理 UEP [delta_uep_eva,0]。
% 以下 6 条细红线仅用于展示，不参与下方原有稳定域边界及稳定性判断。
lambda_family = 0:0.2:1;
M_moon = J*Ws;

% 在足够宽的角度区间内构造无阻尼势能 U_0(delta)。周期判断同时
% 包含负角度侧的限流区域，从而能够得到围住 SEP 的完整闭合边界。
delta_moon = linspace(delta_s - 4*pi, delta_s + 4*pi, 48001);
delta_wrap_moon = mod(delta_moon + pi, 2*pi) - pi;

Pe_moon = (Xsum0*Vgfm*Ug.*sin(delta_moon) ...
    + Rg*(Vgfm^2 - Vgfm*Ug.*cos(delta_moon)) ...
    + Rv0*(Vgfm*Ug.*cos(delta_moon) - Ug^2))./Zsum20;

idx_lim_moon = abs(delta_wrap_moon) > deltacc;
delta_lim_moon = delta_moon(idx_lim_moon);
Den_lim_moon = Vgfm^2 + Ug^2 - 2*Vgfm*Ug*cos(delta_lim_moon);
Rad_lambda_moon = projection_v^2 + Den_lim_moon/Ilim^2 ...
                - (Rg^2 + Xg^2);
Rad_lambda_moon(abs(Rad_lambda_moon) < 1e-12) = 0;
lambda_lim_moon = -projection_v + sqrt(max(Rad_lambda_moon, 0));
lambda_lim_moon = max(lambda_lim_moon, Zv0_abs);
Rv_lim_moon = lambda_lim_moon*cos_phi_v;
Xv_lim_moon = lambda_lim_moon*sin_phi_v;
Rsum_lim_moon = Rg + Rv_lim_moon;
Xsum_lim_moon = Xg + Xv_lim_moon;
Zsum2_lim_moon = Rsum_lim_moon.^2 + Xsum_lim_moon.^2;

Pe_moon(idx_lim_moon) = ( ...
      Xsum_lim_moon.*Vgfm*Ug.*sin(delta_lim_moon) ...
    + Rg*(Vgfm^2 - Vgfm*Ug.*cos(delta_lim_moon)) ...
    + Rv_lim_moon.*(Vgfm*Ug.*cos(delta_lim_moon) - Ug^2)) ...
    ./ Zsum2_lim_moon;

U0_moon = cumtrapz(delta_moon, Pe_moon - Pm);
U0_moon = U0_moon ...
    - interp1(delta_moon, U0_moon, delta_s, 'pchip');

moon_delta_min = inf;
moon_delta_max = -inf;
moon_omega_min = inf;
moon_omega_max = -inf;

for idx_lambda = 1:numel(lambda_family)
    lambda_plot = lambda_family(idx_lambda);

    % 由 Moon 的式 (24)-(25)：
    % omega_{u,lambda} = -lambda*D/M*(delta_{u,lambda}-delta_s)，
    % Pe-Pm + lambda*(1-lambda)*D^2/M*(delta-delta_s) = 0。
    kappa_lambda = lambda_plot*(1-lambda_plot)*D^2/M_moon;
    grad_Ueff = Pe_moon - Pm ...
              + kappa_lambda.*(delta_moon-delta_s);

    idx_right = find( ...
        delta_moon(1:end-1) > delta_s + 1e-8 & ...
        grad_Ueff(1:end-1) >= 0 & grad_Ueff(2:end) <= 0, ...
        1, 'first');

    if isempty(idx_right)
        warning('Moon set skipped for lambda = %.1f: right saddle not found.', ...
            lambda_plot);
        continue;
    end

    d1 = delta_moon(idx_right);
    d2 = delta_moon(idx_right+1);
    g1 = grad_Ueff(idx_right);
    g2 = grad_Ueff(idx_right+1);
    delta_ulambda = d1 - g1*(d2-d1)/(g2-g1);
    omega_ulambda = -lambda_plot*D/M_moon ...
                  * (delta_ulambda-delta_s);

    % E_lambda 的配方形式：
    % E_lambda = M/2*(omega + lambda*D/M*Delta_delta)^2
    %            + U_eff,lambda(delta)。
    Ueff_moon = U0_moon ...
        + 0.5*kappa_lambda.*(delta_moon-delta_s).^2;
    Ecrit_lambda = interp1(delta_moon, Ueff_moon, ...
        delta_ulambda, 'pchip');

    % 左侧交点与右侧鞍点限定包含 SEP 的闭合连通分量。
    level_gap = Ueff_moon - Ecrit_lambda;
    idx_left_all = find( ...
        delta_moon(2:end) < delta_s & ...
        level_gap(1:end-1) >= 0 & level_gap(2:end) <= 0);

    if isempty(idx_left_all)
        warning('Moon set skipped for lambda = %.1f: left closure not found.', ...
            lambda_plot);
        continue;
    end

    idx_left = idx_left_all(end);
    dl1 = delta_moon(idx_left);
    dl2 = delta_moon(idx_left+1);
    h1 = level_gap(idx_left);
    h2 = level_gap(idx_left+1);
    delta_left = dl1 - h1*(dl2-dl1)/(h2-h1);

    delta_curve = linspace(delta_left, delta_ulambda, 1600);
    U0_curve = interp1(delta_moon, U0_moon, delta_curve, 'pchip');
    Ueff_curve = U0_curve ...
        + 0.5*kappa_lambda.*(delta_curve-delta_s).^2;

    radicand = 2*(Ecrit_lambda-Ueff_curve)/M_moon;
    radius_omega = sqrt(max(radicand, 0));
    radius_omega([1 end]) = 0;

    omega_center = -lambda_plot*D/M_moon ...
                 .* (delta_curve-delta_s);
    omega_upper = omega_center + radius_omega;
    omega_lower = omega_center - radius_omega;

    delta_closed = [delta_curve, fliplr(delta_curve)];
    omega_closed = [omega_upper, fliplr(omega_lower)];

    plot(delta_closed, omega_closed, ...
        'Color', [0.85 0.15 0.15], ...
        'LineWidth', 0.7);

    moon_delta_min = min(moon_delta_min, min(delta_closed));
    moon_delta_max = max(moon_delta_max, max(delta_closed));
    moon_omega_min = min(moon_omega_min, min(omega_closed));
    moon_omega_max = max(moon_omega_max, max(omega_closed));
end

% 必要时只扩大显示范围，保证新增曲线完整可见；不改变原曲线数据。
if isfinite(moon_delta_min)
    xl_now = xlim;
    yl_now = ylim;
    x_margin = 0.03*max(moon_delta_max-moon_delta_min, eps);
    y_margin = 0.03*max(moon_omega_max-moon_omega_min, eps);
    xlim([min(xl_now(1), moon_delta_min-x_margin), ...
          max(xl_now(2), moon_delta_max+x_margin)]);
    ylim([min(yl_now(1), moon_omega_min-y_margin), ...
          max(yl_now(2), moon_omega_max+y_margin)]);
end

%% ===================== 原有实际稳定域边界（保持不变） =====================
% 仍采用原程序 lambda = 1 及其 UEP 临界能量进行稳定性判断。
contour(DL, YL, Vlevel, [Vcut Vcut], ...
    'Color', [1 0.25 0.25], ...
    'LineWidth', 2.2);
