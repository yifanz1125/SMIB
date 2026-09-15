%% GFM circle current limiter: piecewise power and energy surface
% This script follows the same structure as the VA-limiter plotting script,
% but uses the circular current limiter power expression.
%
% Required variables in workspace:
% Rg, Xg, Ug, Vgfm, Ilim, Pm, J, Ws, D, prefault_SEP, delta_uep_cl, f1

%% ===================== basic quantities =====================
arg_c = (Vgfm^2 + Ug^2 ...
       - Ilim^2*(Xg^2 + Rg^2))/(2*Vgfm*Ug);

deltac = acos(arg_c);

% Explicitly include all switching points in [-2*pi, 2*pi]
delta0 = -2*pi:0.01:2*pi;
k = -1:1;

delta_sw = [2*pi*k - deltac, 2*pi*k + deltac];
delta_sw = delta_sw(delta_sw >= -2*pi & delta_sw <= 2*pi);

delta = unique(sort([delta0, delta_sw]));

Pv = Rg*(Vgfm^2 - Vgfm*Ug*cos(delta))./(Rg^2+Xg^2) ...
   + Xg*Vgfm*Ug*sin(delta)./(Rg^2+Xg^2);

Den = Vgfm^2 + Ug^2 - 2*Vgfm*Ug*cos(delta);

Rad = Den/Ilim^2 - Xg^2;
Rad(abs(Rad) < 1e-12) = 0;
Rad(Rad < 0) = NaN;

Sroot = sqrt(Rad);

Pi = Sroot ./ Den .* Ilim^2 ...
     .* (Vgfm*Ug*cos(delta) - Ug^2) ...
   + Xg ./ Den .* Ilim^2 .* Vgfm*Ug.*sin(delta) ...
   + Ilim^2*Rg;

Iv = sqrt(Den./(Xg^2+Rg^2));

%% ===================== periodic mapping =====================
delta_wrap = mod(delta + pi, 2*pi) - pi;

tol = 1e-10;

% Include the switching point in both curves
idx_in  = abs(delta_wrap) <= deltac + tol;
idx_out = abs(delta_wrap) >= deltac - tol;

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
plot(delta, Pv_out, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);

Pi_out = Pi;  Pi_out(~idx_out) = NaN;
plot(delta, Pi_out, 'r-', 'LineWidth', 2);

plot(delta, Pm*ones(size(delta)), 'k-', 'LineWidth', 1.2);

axis([-pi pi -2 2]);

xlabel('\delta');
ylabel('P');
grid on;
box on;

%% ===================== construct V3 / V3_2 =====================
delta_s = prefault_SEP(1);

deltacc = deltac;
J_ori = J/Ws;
lamda = 0;

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

%% ===================== Moon lambda-family closed sets (display only) =====================
% Moon's S_lambda is bounded by the equipotential through the saddle point
% X_{u,lambda} of E_lambda, rather than by the physical UEP [delta_uep_cl,0].
% The six thin red curves below are display-only and do not enter the actual
% stability test based on the original red boundary.
lambda_family = 0:0.2:1;
M_moon = J*Ws;

% Construct the undamped potential U_0(delta) over a sufficiently wide
% interval. The periodic mode logic also includes the negative-angle
% current-limited region, which is needed to close each S_lambda boundary.
delta_moon = linspace(delta_s - 4*pi, delta_s + 4*pi, 48001);
delta_wrap_moon = mod(delta_moon + pi, 2*pi) - pi;

Pe_moon = Rg*(Vgfm^2 - Vgfm*Ug*cos(delta_moon))./(Rg^2+Xg^2) ...
        + Xg*Vgfm*Ug*sin(delta_moon)./(Rg^2+Xg^2);

idx_lim_moon = abs(delta_wrap_moon) > deltacc;
delta_lim_moon = delta_moon(idx_lim_moon);
Den_lim_moon = Vgfm^2 + Ug^2 - 2*Vgfm*Ug*cos(delta_lim_moon);
Rad_lim_moon = Den_lim_moon/Ilim^2 - Xg^2;
Rad_lim_moon(abs(Rad_lim_moon) < 1e-12) = 0;
Rad_lim_moon = max(Rad_lim_moon, 0);

Pe_moon(idx_lim_moon) = sqrt(Rad_lim_moon)./Den_lim_moon .* Ilim^2 ...
    .* (Vgfm*Ug*cos(delta_lim_moon) - Ug^2) ...
    + Xg./Den_lim_moon .* Ilim^2 .* Vgfm*Ug.*sin(delta_lim_moon) ...
    + Ilim^2*Rg;

U0_moon = cumtrapz(delta_moon, Pe_moon - Pm);
U0_moon = U0_moon ...
    - interp1(delta_moon, U0_moon, delta_s, 'pchip');

moon_delta_min = inf;
moon_delta_max = -inf;
moon_omega_min = inf;
moon_omega_max = -inf;

for idx_lambda = 1:numel(lambda_family)
    lambda_plot = lambda_family(idx_lambda);

    % From Moon's (24)-(25):
    % omega_{u,lambda} = -lambda*D/M*(delta_{u,lambda}-delta_s),
    % Pe-Pm + lambda*(1-lambda)*D^2/M*(delta-delta_s) = 0.
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

    % Complete-square form of E_lambda:
    % E_lambda = M/2*(omega + lambda*D/M*Delta_delta)^2
    %            + U_eff,lambda(delta).
    Ueff_moon = U0_moon ...
        + 0.5*kappa_lambda.*(delta_moon-delta_s).^2;
    Ecrit_lambda = interp1(delta_moon, Ueff_moon, ...
        delta_ulambda, 'pchip');

    % The left intersection of U_eff,lambda = Ecrit_lambda and the right
    % saddle delimit the connected closed component containing the SEP.
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

% Expand the displayed window only when necessary so the added curves are
% visibly closed; this does not alter any original curve data.
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

%% ===================== original stability boundary (unchanged) =====================
% The actual stability test still uses lambda = 1 and its UEP energy level.
contour(DL, YL, Vlevel, [Vcut Vcut], ...
    'Color', [1 0.25 0.25], ...
    'LineWidth', 2.2);
