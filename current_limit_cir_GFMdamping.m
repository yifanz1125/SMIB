%% GFM circle current limiter: piecewise power and energy surface
% This script follows the same structure as the VA-limiter plotting script,
% but uses the circular current limiter power expression.
%
% Required variables in workspace:
% Rg, Xg, Ug, Vgfm, Ilim, Pm, J, Ws, D, prefault_SEP, delta_uep_cl, f1

%% ===================== basic quantities =====================
delta = -2*pi:0.01:2*pi;

Pv = Rg*(Vgfm^2 - Vgfm*Ug*cos(delta))./(Rg^2+Xg^2) ...
   + Xg*Vgfm*Ug*sin(delta)./(Rg^2+Xg^2);

Den = Vgfm^2 + Ug^2 - 2*Vgfm*Ug*cos(delta);

% circular limiter: active-power expression in the limited region
Rad = Den/Ilim^2 - Xg^2;
Rad(Rad < 0) = NaN;
Sroot = sqrt(Rad);

Pi = Sroot ./ Den .* Ilim^2 .* (Vgfm*Ug*cos(delta) - Ug^2) ...
   + Xg ./ Den .* Ilim^2 .* Vgfm*Ug.*sin(delta) ...
   + Ilim^2*Rg;

Iv = sqrt(Den./(Xg^2+Rg^2));

deltac = acos((Vgfm^2 + Ug^2 - Ilim^2*(Xg^2 + Rg^2))/(2*Vgfm*Ug));

%% ===================== periodic mapping =====================
delta_wrap = mod(delta + pi, 2*pi) - pi;

idx_in  = abs(delta_wrap) <= deltac;
idx_out = ~idx_in;

%% ===================== current plot =====================
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

%% ===================== piecewise power plot =====================
figure;
hold on;

Pv_in  = Pv;  Pv_in(~idx_in)   = NaN;
Pv_out = Pv;  Pv_out(~idx_out) = NaN;

plot(delta, Pv_in,  'b-',  'LineWidth', 2);
plot(delta, Pv_out, 'b--', 'LineWidth', 1.5);

Pi_out = Pi;  Pi_out(~idx_out) = NaN;
plot(delta, Pi_out, 'r-', 'LineWidth', 2);

plot(delta, Pm*ones(size(delta)), 'k-', 'LineWidth', 1.2);

xlabel('\delta');
ylabel('P');
grid on;
box on;

%% ===================== construct V3 / V3_2 =====================
delta_s = prefault_SEP(1);

deltacc = deltac;
J_ori = J/Ws;
lamda = 1;

% inner region energy, i.e., original power-angle relation
VV3 = @(delta_val, omega_val) ...
    J_ori/2*(omega_val*Ws).^2 ...
    - Pm*(delta_val-delta_s) ...
    + Rg*(Vgfm^2*(delta_val-delta_s) - Vgfm*Ug*(sin(delta_val)-sin(delta_s))) / (Rg^2+Xg^2) ...
    - Xg*Vgfm*Ug*(cos(delta_val)-cos(delta_s)) / (Rg^2+Xg^2) ...
    + lamda*D*(omega_val).*(delta_val-delta_s) ...
    + lamda/2*D^2/J/Ws*(delta_val-delta_s).^2;

% limited-region power function for numerical integration
Pi_fun = @(x) ...
    sqrt(max((Vgfm^2 + Ug^2 - 2*Vgfm*Ug*cos(x))/Ilim^2 - Xg^2, 0)) ...
    ./ (Vgfm^2 + Ug^2 - 2*Vgfm*Ug*cos(x)) .* Ilim^2 .* (Vgfm*Ug*cos(x) - Ug^2) ...
    + Xg ./ (Vgfm^2 + Ug^2 - 2*Vgfm*Ug*cos(x)) .* Ilim^2 .* Vgfm*Ug.*sin(x) ...
    + Ilim^2*Rg;

% matching constant: V3_2(deltacc,0) = V3(deltacc,0)
V3_dacc_0 = VV3(deltacc, 0);
Cmatch = V3_dacc_0 + Pm*(deltacc - delta_s) ...
       - lamda/2*D^2/J/Ws*(deltacc-delta_s)^2;

% outer region energy
VV3_2 = @(delta_val, omega_val) ...
    J_ori/2*(omega_val*Ws).^2 ...
    - Pm*(delta_val-delta_s) ...
    + integral(Pi_fun, deltacc, delta_val, 'ArrayValued', true) ...
    + lamda*D*(omega_val).*(delta_val-delta_s) ...
    + lamda/2*D^2/J/Ws*(delta_val-delta_s).^2 ...
    + Cmatch;

%% ===================== 3D energy surface =====================
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

%% ===================== shadow projection for inner region =====================
zmin = -0.2;
shadow = zmin * ones(size(y1));
shadow(~(y1 <= deltacc)) = NaN;

surf(y1, y2, shadow, ...
    'FaceColor', [1 0 0], ...
    'FaceAlpha', 0.15, ...
    'EdgeColor', 'none');

%% ===================== SEP / UEP coordinates =====================
y_sep = 0;
y_uep = 0;

if delta_s <= deltacc
    z_sep = VV3(delta_s, y_sep);
elseif delta_s <= pi
    z_sep = VV3_2(delta_s, y_sep);
else
    z_sep = NaN;
end

% circle-limiter UEP
if delta_uep_cl <= deltacc
    z_uep_cl = VV3(delta_uep_cl, y_uep);
elseif delta_uep_cl <= pi
    z_uep_cl = VV3_2(delta_uep_cl, y_uep);
else
    z_uep_cl = NaN;
end

plot3(delta_s, y_sep, z_sep, 'ko', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'g');

plot3(delta_uep_cl, y_uep, z_uep_cl, 'ko', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'r');

text(delta_s, y_sep, z_sep, '  SEP', ...
    'FontSize', 12, ...
    'Color', 'k', ...
    'FontWeight', 'bold');

text(delta_uep_cl, y_uep, z_uep_cl, '  UEP', ...
    'FontSize', 12, ...
    'Color', 'k', ...
    'FontWeight', 'bold');

%% ===================== energy cut through [delta_uep_cl, 0] =====================
y_cut = 0;

if delta_uep_cl <= deltacc
    Vcut = VV3(delta_uep_cl, y_cut);
elseif delta_uep_cl <= pi
    Vcut = VV3_2(delta_uep_cl, y_cut);
else
    Vcut = NaN;
end

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

%% ===================== level set in f1 =====================
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
