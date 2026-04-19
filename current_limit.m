delta = -2*pi:0.01:2*pi;

Pv = Rg*(Vvfm^2 - Vvfm*Ug*cos(delta))./(Rg^2+Xg^2) ...
   + Xg*Vvfm*Ug*sin(delta)./(Rg^2+Xg^2);

Phi = -pi/4;
Ilim = 2;

Pi = Ilim*Ug*cos(delta+Phi) + Ilim^2*Rg;

Iv = sqrt((Vvfm^2 + Ug^2 - 2*Vvfm*Ug*cos(delta))./(Xg^2+Rg^2));

deltac = acos((Vvfm^2+Ug^2-Ilim^2*(Xg^2+Rg^2))/(2*Vvfm*Ug));
Ilim_sep = -acos((Pin-Ilim^2*Rg)/(Ilim*Ug))-Phi;

% ===== 周期映射 =====
delta_wrap = mod(delta + pi, 2*pi) - pi;

idx_in  = abs(delta_wrap) <= deltac;
idx_out = ~idx_in;
%% ===== 电流图（分段 + 2π周期）=====
figure;
hold on;
plot(delta, Iv, 'r-');
plot(delta, Ilim*ones(size(delta)), 'k-');
plot(deltac, Ilim, 'k.');

%% ===== 功率图（分段 + 不连接）=====
figure;
hold on;

% --- Pv ---
Pv_in  = Pv;  Pv_in(~idx_in)  = NaN;
Pv_out = Pv;  Pv_out(~idx_out)= NaN;

plot(delta, Pv_in,  'b-','linewidth',2);   % 实线
plot(delta, Pv_out, 'b--');  % 虚线

% --- Pi ---
Pi_in  = Pi;  Pi_in(~idx_in)  = NaN;
Pi_out = Pi;  Pi_out(~idx_out)= NaN;

plot(delta, Pi_in,  'r--');  % 虚线
plot(delta, Pi_out, 'r-','linewidth',2);   % 实线

% --- Pm ---
plot(delta, Pm*ones(size(delta)), 'k-');

%%
syms deltax yx;
    delta_s = prefault_SEP(1);
    deltauep = ep_set(2).xep(1);
    V3 = C_dc/4*Kip*yx^2 - Pin*(deltax-delta_s) + Rg*(Vvfm^2*(deltax-delta_s)-Vvfm*Ug*(sin(deltax)-sin(delta_s)))/(Rg^2+Xg^2)-Xg*Vvfm*Ug*(cos(deltax)-cos(delta_s))/(Rg^2+Xg^2) +V32cr-V322cr;%- C_dc/2*Kpp*yx*(Rg*(Vvfm^2-Vvfm*Ug*cos(deltax))/(Rg^2+Xg^2)+Xg*Vvfm*Ug*sin(deltax)/(Rg^2+Xg^2) - Pin)/2;%
    V3_2 = C_dc/4*Kip*yx^2 - Pin*(deltax-delta_s) + Ilim*Ug*sin(deltax+Phi) - Ilim*Ug*sin(delta_s+Phi) + Ilim^2*Rg*(deltax-delta_s) ;
    V3=vpa(V3);
    VV3=matlabFunction(V3);
    VV3_2=matlabFunction(V3_2);
    V3d = jacobian(V3);
    VV3d = matlabFunction(V3d);
    deltacc = acos((Vvfm^2+Ug^2-Ilim^2*(Xg^2+Rg^2))/(2*Vvfm*Ug));
    x1=-deltacc:0.001*pi:pi;
    x2=-4:0.02:6;%-8:0.1:8;
    [y1,y2]=meshgrid(x1,x2);
    zz = zeros(length(x2),length(x1));
    zz_2 = zeros(length(x2),length(x1));
    dzz = zeros(length(x2),length(x1));
    for a = 1: length(x1)
        for b = 1: length(x2)
            zz(b,a) = VV3(y1(b,a), y2(b,a));
            zz_2(b,a) = VV3_2(y1(b,a), y2(b,a));
        end
    end


idx_in  = abs(y1) <= deltacc;
idx_out = ~idx_in;

zz_in = zz;
zz_in(idx_out) = NaN;

zz2_out = zz_2;
zz2_out(idx_in) = NaN;

figure;
hold on;

colormap turbo;          % 比 parula 更有对比
shading interp;

% 🔥 核心：限制颜色范围（增强起伏）
caxis([min(zz(:)) max(zz(:))*0.6]);   % 或自己调 0.4~0.8

colorbar;



h1 = surf(y1, y2, zz_in,   'EdgeColor', 'none', 'FaceColor', 'interp');
h2 = surf(y1, y2, zz2_out, 'EdgeColor', 'none', 'FaceColor', 'interp');

% % 设置透明度
 set([h1 h2], 'FaceAlpha', 0.8);   % 0~1（推荐 0.6~0.9）


xlim([-deltacc, pi]);
ylim([min(x2), max(x2)]);
axis tight;
view(3);
grid on;
box on;

xlabel('\delta');
ylabel('y');
zlabel('V');

zmin = -0.5;
shadow = zmin * ones(size(y1));
shadow(~idx_in) = NaN;

surf(y1, y2, shadow, ...
    'FaceColor',[1 0 0], ...
    'FaceAlpha',0.15, ...
    'EdgeColor','none');
    %V3cr

    
% ===== V3cr 灰色透明平面 =====
xp = [-pi, pi; -pi, pi];
yp = [min(x2), min(x2); max(x2), max(x2)];
zp = V3cr * ones(2,2);

surf(xp, yp, zp, ...
    'FaceColor', [0.5 0.5 0.5], ...
    'FaceAlpha', 0.3, ...
    'EdgeColor', 'none');

% ===== SEP / UEP 坐标 =====
delta_sep = delta_s;
y_sep = 0;

delta_uep = ep_set(2).xep(1);
y_uep = 0;

% ===== 按分段规则计算 z 坐标 =====
if abs(delta_sep) <= deltacc
    z_sep = VV3(delta_sep, y_sep);
else
    z_sep = VV3_2(delta_sep, y_sep);
end

if abs(delta_uep) <= deltacc
    z_uep = VV3(delta_uep, y_uep);
else
    z_uep = VV3_2(delta_uep, y_uep);
end

% ===== 标出 SEP / UEP =====
plot3(delta_sep, y_sep, z_sep, 'ko', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'g');

plot3(delta_uep, y_uep, z_uep, 'ko', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'r');

% ===== 文字标注 =====
text(delta_sep, y_sep, z_sep, '  SEP', ...
    'FontSize', 12, ...
    'Color', 'k', ...
    'FontWeight', 'bold');

text(delta_uep, y_uep, z_uep, '  UEP', ...
    'FontSize', 12, ...
    'Color', 'k', ...
    'FontWeight', 'bold');

xlim([-deltacc, pi]);
ylim([-4, 4]);
zlim([-0.5, 5]);

%% ===== 扫描范围 =====
Xg_min = 0;
Xg_max = 1;
Phi_min = -pi;
Phi_max = pi;

nx = 600;
nphi = 600;

Xg_vec  = linspace(Xg_min, Xg_max, nx);
Phi_vec = linspace(Phi_min, Phi_max, nphi);

[Xg_grid, Phi_grid] = meshgrid(Xg_vec, Phi_vec);

% ===== acos 参数检查 =====
arg_sep = (Pin - Ilim^2*Rg)/(Ilim*Ug);
if abs(arg_sep) > 1
    error('Ilim_sep 中 acos 的输入超出 [-1,1]，当前参数下无实数解。');
end

theta_sep = acos(arg_sep);

arg_c = (Vvfm^2 + Ug^2 - Ilim^2*(Xg_grid.^2 + Rg^2)) ./ (2*Vvfm*Ug);

valid = abs(arg_c) <= 1;

deltac = nan(size(Xg_grid));
deltac(valid) = acos(arg_c(valid));

% ===== 条件：Ilim_sep > -deltac =====
% Ilim_sep = -theta_sep - Phi
Ilim_sep = -theta_sep - Phi_grid;

cond = false(size(Xg_grid));
cond(valid) = Ilim_sep(valid) > -deltac(valid);

% ===== 求边界曲线 Phi_bd(Xg) =====
arg_c_line = (Vvfm^2 + Ug^2 - Ilim^2*(Xg_vec.^2 + Rg^2)) ./ (2*Vvfm*Ug);

Phi_bd = nan(size(Xg_vec));
idx_line = abs(arg_c_line) <= 1;
Phi_bd(idx_line) = acos(arg_c_line(idx_line)) - theta_sep;

% ===== 画图 =====
figure;
hold on; box on; grid on;

imagesc(Xg_vec, Phi_vec, cond);
set(gca, 'YDir', 'normal');

colormap([1 1 1; 0.2 0.7 1]);
cb = colorbar;
cb.Ticks = [0 1];
cb.TickLabels = {'Not satisfied', 'Satisfied'};

plot(Xg_vec, Phi_bd, 'r-', 'LineWidth', 2);

xlabel('X_g');
ylabel('\Phi (rad)');
title('Region satisfying  I_{lim,sep} > -\delta_c');

xlim([Xg_min, Xg_max]);
ylim([Phi_min, Phi_max]);