%% ===================== 基本量 =====================
delta = -2*pi:0.01:2*pi;

Pv = Rg*(Vvfm^2 - Vvfm*Ug*cos(delta))./(Rg^2+Xg^2) ...
   + Xg*Vvfm*Ug*sin(delta)./(Rg^2+Xg^2);

Phi = -pi/4;
Ilim = 2;

Pi = sqrt((Vvfm^2+Ug^2-2*Vvfm*Ug*cos(delta))/Ilim^2 - Xg^2) ...
   ./(Vvfm^2+Ug^2-2*Vvfm*Ug*cos(delta)) .* Ilim^2 .* (Vvfm*Ug*cos(delta)-Ug^2) ...
   + Xg./(Vvfm^2+Ug^2-2*Vvfm*Ug*cos(delta)) .* Ilim^2 .* Vvfm*Ug.*sin(delta) ...
   + Ilim^2*Rg;

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
delta_uep = ep_set(2).xep(1);

% 这里用右侧边界 deltacc
deltacc = acos((Vvfm^2+Ug^2-Ilim^2*(Xg^2+Rg^2))/(2*Vvfm*Ug));

% ------- 内区 V3 的函数句柄 -------
VV3 = @(delta_val, y_val) ...
    C_dc/4*Kip*y_val.^2 ...
    - Pin*(delta_val-delta_s) ...
    + Rg*(Vvfm^2*(delta_val-delta_s) - Vvfm*Ug*(sin(delta_val)-sin(delta_s))) / (Rg^2+Xg^2) ...
    - Xg*Vvfm*Ug*(cos(delta_val)-cos(delta_s)) / (Rg^2+Xg^2) ...
    + V32cr - V322cr;

% ------- Pi 的函数句柄（供数值积分） -------
Pi_fun = @(x) ...
    sqrt((Vvfm^2+Ug^2-2*Vvfm*Ug*cos(x))/Ilim^2 - Xg^2) ...
    ./ (Vvfm^2+Ug^2-2*Vvfm*Ug*cos(x)) .* Ilim^2 .* (Vvfm*Ug*cos(x)-Ug^2) ...
    + Xg ./ (Vvfm^2+Ug^2-2*Vvfm*Ug*cos(x)) .* Ilim^2 .* Vvfm*Ug .* sin(x) ...
    + Ilim^2*Rg;

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

%% ===================== 生成网格 =====================
x1 = -deltacc:0.001*pi:pi;
x2 = -4:0.02:6;

[y1, y2] = meshgrid(x1, x2);

zz   = NaN(size(y1));
zz_2 = NaN(size(y1));

%% ===================== 逐点计算两个曲面 =====================
for a = 1:length(x1)
    for b = 1:length(x2)
        delta_now = y1(b,a);
        y_now     = y2(b,a);

        % 内区：delta <= deltacc
        if delta_now <= deltacc
            zz(b,a) = VV3(delta_now, y_now);
        end

        % 外区：deltacc < delta <= pi
        if delta_now > deltacc && delta_now <= pi
            zz_2(b,a) = VV3_2(delta_now, y_now);
        end
    end
end

%% ===================== 画拼接曲面 =====================
figure;
hold on;

colormap turbo;
shading interp;

% 如果担心 NaN 影响 min/max，就只对有效值取范围
zz_all = [zz(~isnan(zz)); zz_2(~isnan(zz_2))];
caxis([min(zz_all), max(zz_all)*0.6]);

colorbar;

h1 = surf(y1, y2, zz,   'EdgeColor', 'none', 'FaceColor', 'interp');
h2 = surf(y1, y2, zz_2, 'EdgeColor', 'none', 'FaceColor', 'interp');

set([h1 h2], 'FaceAlpha', 0.8);

xlim([-deltacc, pi]);
ylim([min(x2), max(x2)]);
axis tight;
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

%% ===================== V3cr 灰色透明平面 =====================
xp = [-pi, pi; -pi, pi];
yp = [min(x2), min(x2); max(x2), max(x2)];
zp = V3cr * ones(2,2);

surf(xp, yp, zp, ...
    'FaceColor', [0.5 0.5 0.5], ...
    'FaceAlpha', 0.3, ...
    'EdgeColor', 'none');

%% ===================== SEP / UEP 坐标 =====================
y_sep = 0;
y_uep = 0;

% SEP
if delta_s <= deltacc
    z_sep = VV3(delta_s, y_sep);
elseif delta_s <= pi
    z_sep = VV3_2(delta_s, y_sep);
else
    z_sep = NaN;
end

% UEP
if delta_uep <= deltacc
    z_uep = VV3(delta_uep, y_uep);
elseif delta_uep <= pi
    z_uep = VV3_2(delta_uep, y_uep);
else
    z_uep = NaN;
end

%% ===================== 标出 SEP / UEP =====================
% plot3(delta_s, y_sep, z_sep, 'ko', ...
%     'MarkerSize', 8, ...
%     'MarkerFaceColor', 'g');
% 
% plot3(delta_uep, y_uep, z_uep, 'ko', ...
%     'MarkerSize', 8, ...
%     'MarkerFaceColor', 'r');
% 
% text(delta_s, y_sep, z_sep, '  SEP', ...
%     'FontSize', 12, ...
%     'Color', 'k', ...
%     'FontWeight', 'bold');
% 
% text(delta_uep, y_uep, z_uep, '  UEP', ...
%     'FontSize', 12, ...
%     'Color', 'k', ...
%     'FontWeight', 'bold');

xlim([-deltacc, pi]);
ylim([-4, 4]);
zlim([-0.5, 5]);

