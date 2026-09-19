%% Plot experimental delta-omega phase trajectory on the existing f1 figure
% Required workspace variables:
%   delta      : N-by-2 array, [absolute timestamp, phase angle in rad]
%   w_psc      : K-by-2 array, [absolute timestamp, angular speed in rad/s]
%   fault_time : M-by-2 array, [absolute timestamp, fault flag]
%   t_start    : displayed prefault duration used by the time-domain plots
%   t_end      : total displayed duration used by the time-domain plots
%   Wbase      : base angular frequency in rad/s
%   f1         : handle of the existing phase-portrait figure (optional)
%
% The last rising edge of fault_time(:,2) is defined as t = 0. The phase
% trajectory is extracted from t = 0 to t = t_end-t_start, matching the
% endpoint of the existing time-domain plots. The phase angle remains in
% radians, while w_psc is converted to pu by division by Wbase.

%% ===================== Input checks =====================
if ~isnumeric(delta) || size(delta,2) < 2 || isempty(delta)
    error('delta must be a nonempty numeric N-by-2 array.');
end

if ~isnumeric(w_psc) || size(w_psc,2) < 2 || isempty(w_psc)
    error('w_psc must be a nonempty numeric K-by-2 array.');
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

if ~exist('Wbase','var') || ~isnumeric(Wbase) || ...
        ~isscalar(Wbase) || ~isfinite(Wbase) || Wbase <= 0
    error('Wbase must be an existing positive finite scalar.');
end

if any(~isfinite(delta(:,1))) || any(diff(delta(:,1)) < 0)
    error('delta(:,1) must contain finite, nondecreasing timestamps.');
end

if any(~isfinite(w_psc(:,1))) || any(diff(w_psc(:,1)) < 0)
    error('w_psc(:,1) must contain finite, nondecreasing timestamps.');
end

if any(~isfinite(fault_time(:,1))) || any(diff(fault_time(:,1)) < 0)
    error('fault_time(:,1) must contain finite, nondecreasing timestamps.');
end

%% ===================== Locate the last rising edge =====================
fault_signal = fault_time(:,2);
valid_fault_signal = isfinite(fault_signal);

if ~any(valid_fault_signal)
    error('fault_time(:,2) contains no finite samples.');
end

fault_low = min(fault_signal(valid_fault_signal));
fault_high_level = max(fault_signal(valid_fault_signal));

if fault_high_level <= fault_low
    error('No rising edge can be detected because fault_time(:,2) is constant.');
end

% Mid-level threshold works for both 0/1 and other two-level fault signals.
fault_threshold = (fault_low + fault_high_level)/2;
fault_logic = fault_signal > fault_threshold;
idx_fault_rise = find(~fault_logic(1:end-1) & fault_logic(2:end), ...
                      1, 'last') + 1;

if isempty(idx_fault_rise)
    error('No rising edge was found in fault_time(:,2).');
end

fault_start_time = fault_time(idx_fault_rise,1);

%% ===================== Extract the requested interval =====================
window_start_time = fault_start_time;
window_end_time = fault_start_time + (t_end-t_start);

if delta(1,1) > window_start_time || delta(end,1) < window_end_time
    error(['The delta recording does not fully cover the requested interval ' ...
           '[0, t_end-t_start].']);
end

if w_psc(1,1) > window_start_time || w_psc(end,1) < window_end_time
    error(['The w_psc recording does not fully cover the requested interval ' ...
           '[0, t_end-t_start].']);
end

% Select samples nearest to the requested interval boundaries.
[~, idx_delta_start] = min(abs(delta(:,1)-window_start_time));
[~, idx_delta_end] = min(abs(delta(:,1)-window_end_time));
[~, idx_w_start] = min(abs(w_psc(:,1)-window_start_time));
[~, idx_w_end] = min(abs(w_psc(:,1)-window_end_time));

if idx_delta_end < idx_delta_start || idx_w_end < idx_w_start
    error('The extracted delta/w_psc interval is invalid. Check timestamps.');
end

delta_segment = delta(idx_delta_start:idx_delta_end,1:2);
w_segment = w_psc(idx_w_start:idx_w_end,1:2);

% Use relative timestamps to avoid loss of precision for large absolute
% timestamps. The delta timestamps are used as the common plotting grid.
t_delta = delta_segment(:,1)-fault_start_time;
delta_phase = delta_segment(:,2);
t_w = w_segment(:,1)-fault_start_time;
w_value = w_segment(:,2);

phase_duration = t_end-t_start;
valid_delta = isfinite(t_delta) & isfinite(delta_phase) & ...
              t_delta >= 0 & t_delta <= phase_duration;
valid_w = isfinite(t_w) & isfinite(w_value) & ...
          t_w >= 0 & t_w <= phase_duration;
t_delta = t_delta(valid_delta);
delta_phase = delta_phase(valid_delta);
t_w = t_w(valid_w);
w_value = w_value(valid_w);

if isempty(t_delta) || isempty(t_w)
    error('The extracted delta or w_psc interval contains no finite samples.');
end

% interp1 requires unique sample locations. Retain the last value when
% duplicate timestamps occur in the experimental record.
[t_w_unique, idx_w_unique] = unique(t_w,'last');
w_unique = w_value(idx_w_unique);

if numel(t_w_unique) < 2
    error('At least two unique w_psc timestamps are required for alignment.');
end

w_aligned = interp1(t_w_unique, w_unique, t_delta, 'linear', NaN);
valid_common = isfinite(w_aligned);

if ~any(valid_common)
    error('delta and w_psc have no overlapping finite samples.');
end

t_phase = t_delta(valid_common);
delta_phase = delta_phase(valid_common);          % rad
omega_phase_pu = w_aligned(valid_common)/Wbase;   % pu

% Unified output: [relative time, delta (rad), w_psc (pu)].
deltaomega_extracted = [t_phase, delta_phase, omega_phase_pu];

%% ===================== Overlay on the existing f1 =====================
if ~exist('f1','var') || ~isgraphics(f1,'figure')
    f1 = figure(1);
else
    figure(f1);
end

ax = gca;
hold_state = ishold(ax);
hold(ax,'on');

h_deltaomega_exp = plot(ax, delta_phase, omega_phase_pu, ...
    'LineStyle','-', ...
    'LineWidth',1.8, ...
    'Color',[0 0.4470 0.7410], ...
    'DisplayName','Experiment');

uistack(h_deltaomega_exp,'top');

if ~hold_state
    hold(ax,'off');
end

drawnow;
