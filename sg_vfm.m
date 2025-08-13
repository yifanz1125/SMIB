%%
Kpp = 2.5*2*pi;
Kip = 20;


Iq0 = ((Ug*cos(prefault_SEP(1))-Vvfm)*Xg+Rg*Ug*sin(prefault_SEP(1)))/(Rg^2+Xg^2);
%%
f_pos = logspace(0,4,1e4);  
f_neg = -flip(f_pos);       
f_tt = [f_neg, f_pos];      
n_tt = length(f_tt);
s = 1i*2*pi*f_tt;

Gc = (Kip * (2 ./ (s * C_dc)) + Kpp)./s;
Gplant = (Vvfm^2 * Xg) ./ ((Rg + Lg * s).^2 + Xg^2) + Iq0 * Vvfm;
Gop = Gc .* Gplant;

[Gop_mag, Gop_ang] = Fcn_Cal_BodeMagAng(Gop);
[Gplant_mag, Gplant_ang] = Fcn_Cal_BodeMagAng(Gplant);

fbd_L = min(f_pos);
fbd_H = max(f_pos);

%% 波特图绘图（subplot）
figure;
set(gcf,'position',[500 100 1000 500]);


% 正频率 Magnitude
subplot(2,2,2)
semilogx(f_pos, Gplant_mag(n_tt/2+1:end), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]);
hold on;
semilogx(f_pos, Gop_mag(n_tt/2+1:end), 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
grid on;
title('正频率 |G_{op}|');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
xlim([fbd_L fbd_H]);

% 正频率 Phase
subplot(2,2,4)
semilogx(f_pos, Gplant_ang(n_tt/2+1:end), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]);
hold on;
semilogx(f_pos, Gop_ang(n_tt/2+1:end), 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
grid on;
title('正频率 ∠G_{op}');
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
ylim([-180 180]);
xlim([fbd_L fbd_H]);

% 负频率 Magnitude
subplot(2,2,1)
semilogx(f_neg, Gplant_mag(1:n_tt/2), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]);
hold on;
semilogx(f_neg, Gop_mag(1:n_tt/2), 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
grid on;
title('负频率 |G_{op}|');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
xlim([-fbd_H -fbd_L]);

% 负频率 Phase
subplot(2,2,3)
semilogx(f_neg, Gplant_ang(1:n_tt/2), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]);
hold on;
semilogx(f_neg, Gop_ang(1:n_tt/2), 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
grid on;
title('负频率 ∠G_{op}');
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
ylim([-180 180]);
xlim([-fbd_H -fbd_L]);

%% 根轨迹图
Kpp_list = linspace(0, 2.5*2*pi*10,100);

s = tf('s');
Np = 4;  % 极点数（根据系统阶数改）

poles_all = zeros(Np, length(Kpp_list));

% 扫描 Kp 并存储极点
for k = 1:length(Kpp_list)
    Kpp = Kpp_list(k);
    Gc = (Kip * (2 / (s * C_dc)) + Kpp)/s;
    Gplant = (Vvfm^2 * Xg) / ((Rg + Lg * s)^2 + Xg^2) + Iq0 * Vvfm;
    Gop = Gc * Gplant;
    Gcl = feedback(Gop, 1);  % 闭环系统
    
    p_tmp = pole(Gcl);
    % 按实部排序，便于追踪轨迹
    [~, idx] = sort(real(p_tmp), 'descend');
    poles_all(:, k) = p_tmp(idx);
end

% 轨迹绘图（彩色）
figure; hold on; grid on;

% colormap 与 colorbar 设置（turbo: 蓝 → 品红）
cmap = turbo(length(Kpp_list));  % 或 cool()

% 绘制所有极点轨迹（同一Kpp颜色相同）
for k = 1:length(Kpp_list)
    % 当前所有模态的极点（实数 + 虚数）
    real_parts = real(poles_all(:,k));
    imag_parts = imag(poles_all(:,k));
    
    plot(real_parts, imag_parts, 'x', ...
        'Color', cmap(k,:), ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6);
end

% 添加模态终点标注（λ标签）
for i = 1:Np
    text(real(poles_all(i,end)) + 0.5, imag(poles_all(i,end)), ...
        ['\lambda_' num2str(i)], 'FontSize', 12, 'FontWeight', 'bold');
end

% 坐标轴样式与背景
xlabel('Real (s^{-1})');
ylabel('Imag (s^{-1})');
title('Closed-loop Pole Trajectories vs. K_{pp}');

xlim([-60 20]);
ylim([-400 400]);
xline(0, '--k', 'LineWidth', 1.2);

% 非稳定区背景色
fill([0 20 20 0], [-500 -500 500 500], [1 0.9 0.8], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.3);

% 模态文本标签（可自定义位置）
text(-10, 150, 'LFO mode', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'w');
text(-55, 300, 'SO mode',  'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'w');

% Colorbar 设置
colormap(cmap);
cb = colorbar;
cb.Label.String = 'DVC K_{pp} (p.u.)';
cb.Ticks = [0 1];
cb.TickLabels = {num2str(Kpp_list(1), '%.2f'), num2str(Kpp_list(end), '%.2f')};


%%
Kip_list = linspace(0.1, 100,100);
Kpp = 2.5*2*pi;
s = tf('s');
Np = 4;  % 极点数（根据系统阶数改）

poles_all = zeros(Np, length(Kpp_list));

% 扫描 Kp 并存储极点
for k = 1:length(Kip_list)
    Kip = Kip_list(k);
    Gc = (Kip * (2 / (s * C_dc)) + Kpp)/s;
    Gplant = (Vvfm^2 * Xg) / ((Rg + Lg * s)^2 + Xg^2) + Iq0 * Vvfm;
    Gop = Gc * Gplant;
    Gcl = feedback(Gop, 1);  % 闭环系统
    
    p_tmp = pole(Gcl);
    % 按实部排序，便于追踪轨迹
    [~, idx] = sort(real(p_tmp), 'descend');
    poles_all(:, k) = p_tmp(idx);
end

% 轨迹绘图（彩色）
figure; hold on; grid on;

% colormap 与 colorbar 设置（turbo: 蓝 → 品红）
cmap = turbo(length(Kip_list));  % 或 cool()

% 绘制所有极点轨迹（同一Kpp颜色相同）
for k = 1:length(Kip_list)
    % 当前所有模态的极点（实数 + 虚数）
    real_parts = real(poles_all(:,k));
    imag_parts = imag(poles_all(:,k));
    
    plot(real_parts, imag_parts, 'x', ...
        'Color', cmap(k,:), ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6);
end

% 添加模态终点标注（λ标签）
for i = 1:Np
    text(real(poles_all(i,end)) + 0.5, imag(poles_all(i,end)), ...
        ['\lambda_' num2str(i)], 'FontSize', 12, 'FontWeight', 'bold');
end

% 坐标轴样式与背景
xlabel('Real (s^{-1})');
ylabel('Imag (s^{-1})');
title('Closed-loop Pole Trajectories vs. K_{ip}');

xlim([-60 20]);
ylim([-400 400]);
xline(0, '--k', 'LineWidth', 1.2);

% 非稳定区背景色
fill([0 20 20 0], [-500 -500 500 500], [1 0.9 0.8], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.3);

% 模态文本标签（可自定义位置）
text(-10, 150, 'LFO mode', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'w');
text(-55, 300, 'SO mode',  'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'w');

% Colorbar 设置
colormap(cmap);
cb = colorbar;
cb.Label.String = 'DVC K_{pp} (p.u.)';
cb.Ticks = [0 1];
cb.TickLabels = {num2str(Kip_list(1), '%.2f'), num2str(Kip_list(end), '%.2f')};

%%
function [Mag,Ang]=Fcn_Cal_BodeMagAng(Input)
    Mag=20*log10(abs(Input));
    Ang=angle(Input)*180/pi;
end