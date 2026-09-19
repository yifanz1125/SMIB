%% Plot experimental phase angle in a new figure
% Required workspace variables:
%   delta      : N-by-2 array, [absolute timestamp, phase angle in rad]
%   fault_time : M-by-2 array, [absolute timestamp, fault flag]
%   t_start    : displayed prefault duration
%   t_end      : total displayed duration
%
% The last rising edge of fault_time(:,2) is defined as t = 0. The phase
% angle is converted from radians to degrees and wrapped into [0, 360).
t_start = 0.4;
t_c = 0;
angle_wrap_mode = '-180_180';      % '0_360' 或 '-180_180'
%% ===================== Input checks =====================
if ~isnumeric(delta) || size(delta,2) < 2 || isempty(delta)
    error('delta must be a nonempty numeric N-by-2 array.');
end

if ~isnumeric(fault_time) || size(fault_time,2) < 2 || isempty(fault_time)
    error('fault_time must be a nonempty numeric M-by-2 array.');
end

if ~exist('t_start','var') || ~exist('t_end','var')
    error('t_start and t_end must already exist in the workspace.');
end

if t_start < 0 || t_end <= t_start
    error('Require t_start >= 0 and t_end > t_start.');
end

if any(~isfinite(delta(:,1))) || any(diff(delta(:,1)) < 0)
    error('delta(:,1) must contain finite, nondecreasing timestamps.');
end

if any(~isfinite(fault_time(:,1))) || any(diff(fault_time(:,1)) < 0)
    error('fault_time(:,1) must contain finite, nondecreasing timestamps.');
end

%% ===================== Locate the first rising edge =====================
fault_signal = fault_time(:,2);
valid_fault_signal = isfinite(fault_signal);

if ~any(valid_fault_signal)
    error('fault_time(:,2) contains no finite samples.');
end

fault_low = min(fault_signal(valid_fault_signal));
fault_high = max(fault_signal(valid_fault_signal));

if fault_high <= fault_low
    error('No rising edge can be detected because fault_time(:,2) is constant.');
end

% This threshold works for both 0/1 and other two-level fault signals.
fault_threshold = (fault_low + fault_high)/2;
fault_logic = fault_signal > fault_threshold;
idx_fault_rise = find(~fault_logic(1:end-1) & fault_logic(2:end), ...
                      1, 'last') + 1;

if isempty(idx_fault_rise)
    error('No rising edge was found in fault_time(:,2).');
end

fault_start_time = fault_time(idx_fault_rise,1);

%% ===================== Extract and process delta =====================
window_start_time = fault_start_time - t_start;
window_end_time = fault_start_time + (t_end - t_start);

if delta(1,1) > window_start_time || delta(end,1) < window_end_time
    error(['The delta recording does not fully cover the requested interval ' ...
           '[-t_start, t_end-t_start].']);
end

% Select the samples nearest to the requested interval boundaries.
[~, idx_delta_start] = min(abs(delta(:,1) - window_start_time));
[~, idx_delta_end] = min(abs(delta(:,1) - window_end_time));

if idx_delta_end < idx_delta_start
    error('The extracted delta interval is invalid. Check the timestamps.');
end

idx_delta_window = idx_delta_start:idx_delta_end;
delta_extracted = delta(idx_delta_window,1:2);

% Unified format: [time relative to fault start, phase angle in degrees].
delta_extracted(:,1) = delta_extracted(:,1) - fault_start_time;
delta_deg = delta_extracted(:,2)*180/pi;

switch angle_wrap_mode
    case '0_360'
        delta_extracted(:,2) = mod(delta_deg,360);
        angle_lower_limit = 0;
        angle_upper_limit = 360;

    case '-180_180'
        delta_extracted(:,2) = mod(delta_deg+180,360)-180;
        angle_lower_limit = -180;
        angle_upper_limit = 180;

    otherwise
        error(['angle_wrap_mode must be either ''0_360'' ' ...
               'or ''-180_180''.']);
end

valid_delta = isfinite(delta_extracted(:,2));
if ~any(valid_delta)
    error('The extracted phase-angle data contain no finite samples.');
end

%% ===================== Plot in a new figure =====================
f10 = figure(10);
clf(f10);
set(f10, 'Position', [680 558 1010 200]);

gem_colors = orderedcolors("gem");
gem_blue = gem_colors(1,:);

plot(delta_extracted(:,1), delta_extracted(:,2), ...
    'LineStyle','-', ...
    'LineWidth',2.5, ...
    'Color',gem_blue);

grid on;
xlim([-t_start, t_end-t_start]);
xticks(-t_start:0.2:t_end-t_start);

% Adaptive y-axis with a 5% margin, while remaining within 0--360 degrees.
y_data_min = min(delta_extracted(valid_delta,2));
y_data_max = max(delta_extracted(valid_delta,2));
y_span = y_data_max - y_data_min;

if y_span > 0
    y_margin = max(0.05*y_span, 1);
else
    y_margin = 5;
end

y_lower = max(angle_lower_limit, y_data_min-y_margin);
y_upper = min(angle_upper_limit, y_data_max+y_margin);

if y_upper <= y_lower
    y_lower = max(angle_lower_limit, y_data_min-5);
    y_upper = min(angle_upper_limit, y_data_max+5);
end

ylim([y_lower, y_upper]);

ax = gca;

% Only the complete 0--360-degree range is treated specially. For every
% other range, MATLAB selects the y ticks automatically according to the
% axes height and font size.
% Use 90-degree ticks when the complete angular range is displayed.
if y_lower == angle_lower_limit && y_upper == angle_upper_limit
    yticks(ax, angle_lower_limit:90:angle_upper_limit);
else
    ax.YTickMode = 'auto';
end

% Shade the fault interval over the complete y-axis range. The patch is
% explicitly placed below every plotted curve, while grid lines stay on top.
if exist('t_c','var') && isscalar(t_c) && isfinite(t_c) && t_c > 0
    hold_state = ishold(ax);
    hold(ax,'on');
    h_fault_area = fill(ax, [0 t_c t_c 0], ...
        [y_lower y_lower y_upper y_upper], ...
        [.9805 .7031 .6797], ...
        'LineStyle','none', ...
        'FaceAlpha',0.5, ...
        'HandleVisibility','off');
    uistack(h_fault_area,'bottom');

    if ~hold_state
        hold(ax,'off');
    end
end

ax.Layer = 'top';
set(ax, 'TickLabelInterpreter','latex', 'FontSize',18);
ylabel(ax,'$\delta\,(\,^{\circ}\,)$', ...
    'Interpreter','latex','FontSize',26);

% Move ylabel slightly closer to the y-axis. Increase 0.015 to move it
% further right (closer); decrease it to move the label further left.
ax.YLabel.Units = 'normalized';
ylabel_position = ax.YLabel.Position;
ylabel_position(1) = ylabel_position(1) + 0.005;
ax.YLabel.Position = ylabel_position;
box on;
drawnow;


