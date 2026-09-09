%% ===================== 基本量 =====================
delta = -2*pi:0.01:2*pi;

Pv = Rg*(Vvfm^2 - Vvfm*Ug*cos(delta))./(Rg^2+Xg^2) ...
   + Xg*Vvfm*Ug*sin(delta)./(Rg^2+Xg^2);

Phi = -pi/4;
Ilim = 2;

Den = Vvfm^2 + Ug^2 - 2*Vvfm*Ug*cos(delta);
Xvar = sqrt(Den/Ilim^2 - Rg^2);

Pi = Rg ./ Den .* Ilim^2 .* (Vvfm^2 - Vvfm*Ug*cos(delta)) ...
   + Xvar ./ Den .* Ilim^2 .* Vvfm*Ug.*sin(delta);

Iv = sqrt((Vvfm^2 + Ug^2 - 2*Vvfm*Ug*cos(delta))./(Xg^2+Rg^2));

deltac = acos((Vvfm^2+Ug^2-Ilim^2*(Xg^2+Rg^2))/(2*Vvfm*Ug));
Ilim_sep = -acos((Pin-Ilim^2*Rg)/(Ilim*Ug))-Phi;

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
plot(delta, Pv_out, 'b--', 'LineWidth', 1.5);

% --- Pi ---
Pi_in  = Pi;  Pi_in(~idx_in)   = NaN;
Pi_out = Pi;  Pi_out(~idx_out) = NaN;

plot(delta, Pi_out, 'r-', 'LineWidth', 2);

% --- Pm ---
plot(delta, Pm*ones(size(delta)), 'k-', 'LineWidth', 1.2);

xlabel('\delta');
ylabel('P');
grid on;
box on;

%% ===================== 构造含阻尼补偿的 V3 =====================

delta_s = prefault_SEP(1);

% 右侧限流切换点，也作为功率极大点分段点
deltacc = acos((Vvfm^2+Ug^2-Ilim^2*(Xg^2+Rg^2))/(2*Vvfm*Ug));

deltauep = delta_uep_va;

Den_grid = Rg^2 + Xg^2;

% ---------- 不限流功率 Pv ----------
Pv_fun = @(d) ...
    Rg*(Vvfm^2 - Vvfm*Ug*cos(d))./Den_grid ...
  + Xg*Vvfm*Ug*sin(d)./Den_grid;

% ---------- 限流功率 Pi ----------
Pi_fun = @(d) ...
    Rg ./ (Vvfm^2 + Ug^2 - 2*Vvfm*Ug*cos(d)) .* Ilim^2 .* ...
    (Vvfm^2 - Vvfm*Ug*cos(d)) ...
  + sqrt((Vvfm^2 + Ug^2 - 2*Vvfm*Ug*cos(d))/Ilim^2 - Rg^2) ...
    ./ (Vvfm^2 + Ug^2 - 2*Vvfm*Ug*cos(d)) .* ...
    Ilim^2 .* Vvfm*Ug.*sin(d);

% ---------- 不限流势能积分：int_delta_s^d (Pv - Pin) ddelta ----------
Pint_v = @(d) ...
    - Pin*(d-delta_s) ...
    + Rg*(Vvfm^2*(d-delta_s) ...
    - Vvfm*Ug*(sin(d)-sin(delta_s))) / Den_grid ...
    - Xg*Vvfm*Ug*(cos(d)-cos(delta_s)) / Den_grid;

% ---------- 拼接常数，保证 deltacc 处连续 ----------
V_dacc_0 = Pint_v(deltacc);
Cmatch = V_dacc_0 + Pin*(deltacc - delta_s);

% ---------- 分段势能 Vp(delta) ----------
Vp_fun = @(d) Vp_piecewise(d, deltacc, delta_s, Pin, Pint_v, Pi_fun, Cmatch);

% ---------- 分段积分 int_da^db (P_piecewise - Pin) ddelta ----------
IntP_piece = @(da,db) Vp_fun(db) - Vp_fun(da);

% ---------- UEP处限流功率斜率 ----------
h = 1e-6;
dPi_uep = (Pi_fun(deltauep+h) - Pi_fun(deltauep-h))/(2*h);

% ---------- UEP稳定特征向量斜率 ----------
s_stable = ...
    (Kpp*dPi_uep - sqrt(Kpp^2*dPi_uep^2 - 4*Kip*dPi_uep))/2;

% ---------- 阻尼补偿系数 ----------
coef = Kpp/(2/C_dc*Kip);

% ---------- 总能量函数 ----------
Vtotal = @(d,y) V_total_limited( ...
    d, y, C_dc, Kip, coef, deltacc, deltauep, ...
    Vp_fun, IntP_piece, s_stable);

%% ===================== 画三维能量面 =====================

dx = 0.01*pi;
x1 = unique(sort([-deltacc:dx:pi, deltacc]));
x2 = -4:0.02:6;

[y1, y2] = meshgrid(x1, x2);
zz_total = NaN(size(y1));

for a = 1:length(x1)
    for b = 1:length(x2)

        delta_now = y1(b,a);
        y_now     = y2(b,a);

        zz_total(b,a) = Vtotal(delta_now, y_now);

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
ylim([-4, 4]);
zlim([-0.5, 5]);

view(3);
grid on;
box on;
xlabel('\delta');
ylabel('y');
zlabel('V');

%% ===================== SEP / UEP =====================

y_sep = 0;
y_uep = 0;

z_sep = Vtotal(delta_s, y_sep);
z_uep_va = Vtotal(delta_uep_va, y_uep);

plot3(delta_s, y_sep, z_sep, 'ko', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'g');

plot3(delta_uep_va, y_uep, z_uep_va, 'ko', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'r');

text(delta_s, y_sep, z_sep, '  SEP', ...
    'FontSize', 12, ...
    'Color', 'k', ...
    'FontWeight', 'bold');

text(delta_uep_va, y_uep, z_uep_va, '  UEP', ...
    'FontSize', 12, ...
    'Color', 'k', ...
    'FontWeight', 'bold');

%% ===================== 经过 UEP 的能量切面 =====================

Vcut = Vtotal(delta_uep_va, 0);

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

%% ===================== 在 f1 里画 level set =====================

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

        Vlevel(b,a) = Vtotal(delta_now, y_now);

    end
end

contour(DL, YL, Vlevel, [Vcut Vcut], ...
    'Color', [0.25 0.25 0.25], ...
    'LineWidth', 2.2);

%% ===================== 本段所需局部函数 =====================

function Vp = Vp_piecewise(d, deltacc, delta_s, Pin, Pint_v, Pi_fun, Cmatch)

    if d <= deltacc
        Vp = Pint_v(d);
    else
        Vp = - Pin*(d-delta_s) ...
             + integral(Pi_fun, deltacc, d, 'ArrayValued', true) ...
             + Cmatch;
    end

end

function V = V_total_limited(d, y, C_dc, Kip, coef, deltacc, deltauep, ...
                             Vp_fun, IntP_piece, s_stable)

    Vp = Vp_fun(d);

    xint = Kip*y;

    if d <= deltacc

        I1 = IntP_piece(d, deltacc);
        I2 = IntP_piece(deltacc, deltauep);

        Wd = ...
            - coef * xint/(deltauep-d) * I1 ...
            + coef * (s_stable) * I2;

    else

        I2 = IntP_piece(d, deltauep);

        Wd = ...
            + coef * (s_stable) * I2;

    end

    V = C_dc/4*Kip*y^2 + Vp + Wd;

end