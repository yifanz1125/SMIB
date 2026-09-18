%% Overlay experimental results on the existing GFM figures
% Run the analytical/Simulink plotting scripts first so that f1, f2 and f3
% already contain the existing curves. The experimental Scope block must
% use the "Structure With Time" format with decimation equal to 1:
%   signal 1: delta2deg (degree)
%   signal 2: omega2IBR (the measured value contains the nominal 1 pu)

required_vars_exp = {'ScopeData_exp','Tc','t_sim_start','t_start','t_end', ...
                     'f1','f2','f3'};
for k_required_exp = 1:numel(required_vars_exp)
    if ~exist(required_vars_exp{k_required_exp},'var')
        error('plot_GFM_cir_exp:MissingVariable', ...
            'Required workspace variable "%s" is missing.', ...
            required_vars_exp{k_required_exp});
    end
end

if ~isfield(ScopeData_exp,'signals') || numel(ScopeData_exp.signals) < 2
    error('plot_GFM_cir_exp:InvalidScopeData', ...
        'ScopeData_exp must contain at least two signals.');
end
if ~isfield(ScopeData_exp,'time')
    error('plot_GFM_cir_exp:MissingTime', ...
        'ScopeData_exp must contain the logged time vector.');
end
if ~isscalar(Tc) || ~isfinite(Tc) || Tc <= 0
    error('plot_GFM_cir_exp:InvalidSampleTime', ...
        'Tc must be a positive finite scalar.');
end

%% ===================== Experimental scope channels =====================
DeltaGFM_exp = ScopeData_exp.signals(1).values(:);       % delta (degree)
OmegaGFM_exp = ScopeData_exp.signals(2).values(:) - 1;   % Delta omega (pu)
t_GFM_exp = ScopeData_exp.time(:);

n_scope_exp = numel(t_GFM_exp);
signal_lengths_exp = [numel(DeltaGFM_exp), numel(OmegaGFM_exp)];
if any(signal_lengths_exp ~= n_scope_exp)
    error('plot_GFM_cir_exp:LengthMismatch', ...
        ['ScopeData_exp.time and experimental channels 1 and 2 ' ...
         'must have equal lengths.']);
end

%% ===================== Common time window =====================
% Decimation is 1, hence the logged interval is exactly Tc.
T_deta_exp = Tc;

% Account for a Scope time vector that does not start exactly at zero.
idx_start_exp = round((t_sim_start - t_GFM_exp(1))/T_deta_exp) + 1;
n_window_exp  = round(t_end/T_deta_exp);
idx_end_exp   = idx_start_exp + n_window_exp;

if idx_start_exp < 1 || idx_end_exp > n_scope_exp || ...
        idx_end_exp < idx_start_exp
    error('plot_GFM_cir_exp:InvalidTimeWindow', ...
        ['Requested samples [%d,%d] exceed the available ' ...
         'ScopeData_exp range [1,%d].'], ...
        idx_start_exp, idx_end_exp, n_scope_exp);
end

idx_experiment = idx_start_exp:idx_end_exp;

% Use Tc to reconstruct the uniform experimental time base. Subtracting
% t_start places the disturbance instant at t = 0, as in the existing plots.
t_experiment = (0:n_window_exp)'*T_deta_exp;
t_experiment_plot = t_experiment - t_start;

% Phase portrait uses radians; the phase-angle time trace remains in degrees.
delta_experiment = DeltaGFM_exp(idx_experiment)*pi/180;
omega_experiment = OmegaGFM_exp(idx_experiment);

experimental_color = '#0072BD';

%% ===================== f1: experimental phase portrait =====================
figure(f1);
hold on;
h_exp_f1 = plot(delta_experiment, omega_experiment, ...
    'LineStyle','-','LineWidth',2,'Color',experimental_color, ...
    'Tag','plot_GFM_cir_experiment');

%% ===================== f2: experimental phase-angle response =====================
figure(f2);
hold on;
h_exp_f2 = plot(t_experiment_plot, DeltaGFM_exp(idx_experiment), ...
    'LineStyle','-','LineWidth',2,'Color',experimental_color, ...
    'Tag','plot_GFM_cir_experiment');
xlim([-t_start, t_end-t_start]);
xticks(-t_start:0.2:t_end-t_start);
refresh_experimental_time_axis(gca);

%% ===================== f3: experimental frequency response =====================
figure(f3);
hold on;
h_exp_f3 = plot(t_experiment_plot, omega_experiment, ...
    'LineStyle','-','LineWidth',2,'Color',experimental_color, ...
    'Tag','plot_GFM_cir_experiment');
xlim([-t_start, t_end-t_start]);
xticks(-t_start:0.2:t_end-t_start);
refresh_experimental_time_axis(gca);


function refresh_experimental_time_axis(ax)
% Re-evaluate the y range after overlaying the experimental curve, then
% extend any existing fault-on shading patch to the updated limits.
    ylim(ax,'auto');
    drawnow;
    yl = ylim(ax);

    h_patch = findobj(ax,'Type','Patch');
    for k_patch = 1:numel(h_patch)
        x_patch = get(h_patch(k_patch),'XData');
        if numel(x_patch) == 4
            set(h_patch(k_patch),'YData',[yl(1), yl(1), yl(2), yl(2)]);
        end
    end
    ylim(ax,yl);
end
