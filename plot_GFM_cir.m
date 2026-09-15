%% Plot Simulink results for the GFM circular-current-limit case
% Run the analytical script first so that f1, f2, f3, t_start and t_end
% already exist. ScopeData must use the "Structure With Time" format:
%   signal 1: delta1deg,  signal 2: omega1IBR,
%   signal 3: delta2deg,  signal 4: omega2IBR.

required_vars = {'ScopeData','Ts','t_sim_start','t_start','t_end', ...
                 'f1','f2','f3'};
for k_required = 1:numel(required_vars)
    if ~exist(required_vars{k_required},'var')
        error('plot_GFM_cir:MissingVariable', ...
            'Required workspace variable "%s" is missing.', ...
            required_vars{k_required});
    end
end

if ~isfield(ScopeData,'signals') || numel(ScopeData.signals) < 4
    error('plot_GFM_cir:InvalidScopeData', ...
        'ScopeData must contain at least four signals.');
end

%% ===================== Scope channels =====================
DeltaGFM  = ScopeData.signals(1).values(:);       % delta1 (degree)
OmegaGFM  = ScopeData.signals(2).values(:) - 1;   % Delta omega1 (pu)
t_GFM = ScopeData.time(:);

%% ===================== Common time window =====================
% Scope recording resolution/decimation interval.
T_deta = Ts*10;

t_end2 = t_sim_start + t_end;

idx_start = round(t_sim_start/T_deta) + 1;
idx_end   = round(t_end2/T_deta) + 1;

n_scope = numel(t_GFM);
signal_lengths = [numel(DeltaGFM), numel(OmegaGFM)];
if any(signal_lengths ~= n_scope)
    error('plot_GFM_cir:LengthMismatch', ...
        'ScopeData.time and the four signal channels must have equal lengths.');
end
if idx_start < 1 || idx_end > n_scope || idx_end < idx_start
    error('plot_GFM_cir:InvalidTimeWindow', ...
        ['Requested samples [%d,%d] exceed the available ScopeData ' ...
         'range [1,%d].'], idx_start, idx_end, n_scope);
end

idx_simulation = idx_start:idx_end;

% The two simulations share the same time vector. The additional shift by
% t_start places the disturbance instant at t = 0.
t_simulation = t_GFM(idx_simulation) - t_sim_start;
t_simulation_plot = t_simulation - t_start;

% Phase portraits use radians; f2 retains phase angle in degrees.
delta_simulation  = DeltaGFM(idx_simulation)*pi/180;
omega_simulation  = OmegaGFM(idx_simulation);

simulation_color = '#A2142F';

%% ===================== f1: phase portraits =====================
figure(f1);
hold on;
delete(findobj(f1,'Tag','plot_GFM_cir_simulation'));
h_sim_f1(1) = plot(delta_simulation, omega_simulation, ...
    'LineStyle','-','LineWidth',2,'Color',simulation_color, ...
    'Tag','plot_GFM_cir_simulation');


%% ===================== f2: phase-angle responses =====================
figure(f2);
hold on;
delete(findobj(f2,'Tag','plot_GFM_cir_simulation'));
h_sim_f2(1) = plot(t_simulation_plot, DeltaGFM(idx_simulation), ...
    'LineStyle','-','LineWidth',2,'Color',simulation_color, ...
    'Tag','plot_GFM_cir_simulation');
xlim([-t_start, t_end-t_start]);
xticks(-t_start:0.2:t_end-t_start);
refresh_time_axis(gca);

%% ===================== f3: frequency responses =====================
figure(f3);
hold on;
delete(findobj(f3,'Tag','plot_GFM_cir_simulation'));
h_sim_f3(1) = plot(t_simulation_plot, omega_simulation, ...
    'LineStyle','-','LineWidth',2,'Color',simulation_color, ...
    'Tag','plot_GFM_cir_simulation');
xlim([-t_start, t_end-t_start]);
xticks(-t_start:0.2:t_end-t_start);
refresh_time_axis(gca);


function refresh_time_axis(ax)
% Re-evaluate the automatic y range after adding Simulink curves and extend
% an existing fault-on shading patch to the new limits.
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
