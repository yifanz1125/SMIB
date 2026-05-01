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

%% ===================== 构造 V3 / V3_2 =====================
% 平衡点
delta_s   = prefault_SEP(1);

% 这里用右侧边界 deltacc
deltacc = acos((Vvfm^2+Ug^2-Ilim^2*(Xg^2+Rg^2))/(2*Vvfm*Ug));

% ------- 内区 V3 的函数句柄 -------
VV3 = @(delta_val, y_val) ...
    C_dc/4*Kip*y_val.^2 ...
    - Pin*(delta_val-delta_s) ...
    + Rg*(Vvfm^2*(delta_val-delta_s) - Vvfm*Ug*(sin(delta_val)-sin(delta_s))) / (Rg^2+Xg^2) ...
    - Xg*Vvfm*Ug*(cos(delta_val)-cos(delta_s)) / (Rg^2+Xg^2);

% ------- Pi 的函数句柄（供数值积分） -------
Pi_fun = @(x) ...
    Rg ./ (Vvfm^2 + Ug^2 - 2*Vvfm*Ug*cos(x)) .* Ilim^2 .* (Vvfm^2 - Vvfm*Ug*cos(x)) ...
    + sqrt((Vvfm^2 + Ug^2 - 2*Vvfm*Ug*cos(x))/Ilim^2 - Rg^2) ...
    ./ (Vvfm^2 + Ug^2 - 2*Vvfm*Ug*cos(x)) .* Ilim^2 .* Vvfm*Ug.*sin(x);
% ------- 拼接常数 -------
% 要求：V3_2(deltacc,0) = V3(deltacc,0)
% 因为积分上限=下限时积分为0
V3_dacc_0 = VV3(deltacc, 0);
Cmatch = V3_dacc_0 + Pin*(deltacc - delta_s);

% ------- 外区 V3_2 的函数句柄 -------
VV3_2 = @(delta_val, y_val) ...
    C_dc/4*Kip*y_val.^2 ...
    - Pin*(delta_val-delta_s) ...
    + integral(Pi_fun, deltacc, delta_val, 'ArrayValued', true) ...
    + Cmatch;

%% =================================
dx = 0.01*pi;
x1 = unique(sort([-deltacc:dx:pi, deltacc]));
x2 = -4:0.02:6;

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
ylim([-4, 4]);
zlim([-0.5, 5]);

view(3);
grid on;
box on;
xlabel('\delta');
ylabel('y');
zlabel('V');

%% ===================== 阴影投影（只给内区） =====================
zmin = -0.5;
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

% --- VA 模型 UEP ---
if delta_uep_va <= deltacc
    z_uep_va = VV3(delta_uep_va, y_uep);
elseif delta_uep_va <= pi
    z_uep_va = VV3_2(delta_uep_va, y_uep);
else
    z_uep_va = NaN;
end

% ===================== 标出 SEP / UEP =====================

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


%% ===================== 经过 [delta_uep_cl, 0] 的能量切面 =====================

y_cut = 0;

if delta_uep_va <= deltacc
    Vcut = VV3(delta_uep_va, y_cut);
elseif delta_uep_va <= pi
    Vcut = VV3_2(delta_uep_va, y_cut);
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

contour(DL, YL, Vlevel, [Vcut Vcut], ...
    'Color', [0.25 0.25 0.25], ...
    'LineWidth', 2.2);
