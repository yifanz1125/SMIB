%% Plot experimental delta on the existing f2 figure
% Required workspace variables:
%   delta      : N-by-2 array, [absolute timestamp, phase angle]
%   fault_time : M-by-2 array, [absolute timestamp, fault flag]
%   t_start    : displayed prefault duration
%   t_end      : total displayed duration
%   f2         : handle of the existing phase-angle figure (optional)
%
% Both timestamp columns must use the same time unit as t_start and t_end.
% No angle-unit conversion is applied. Therefore, delta(:,2) must use the
% same unit as the existing f2 vertical axis (normally degrees).

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

%% ===================== Extract and align delta =====================
window_start_time = fault_start_time - t_start;
window_end_time = fault_start_time + (t_end - t_start);

if delta(1,1) > window_start_time || delta(end,1) < window_end_time
    error(['The delta recording does not fully cover the requested interval ' ...
           '[-t_start, t_end-t_start].']);
end

% Use the samples nearest to the two requested boundaries. This avoids
% losing an endpoint because large absolute timestamps have finite precision.
[~, idx_delta_start] = min(abs(delta(:,1) - window_start_time));
[~, idx_delta_end] = min(abs(delta(:,1) - window_end_time));

if idx_delta_end < idx_delta_start
    error('The extracted delta interval is invalid. Check the timestamps.');
end

idx_delta_window = idx_delta_start:idx_delta_end;

% Unified experimental format: [time relative to fault start, delta].
delta_extracted = delta(idx_delta_window,1:2);
delta_extracted(:,1) = delta_extracted(:,1) - fault_start_time;

%% ===================== Overlay on the existing f2 =====================
if ~exist('f2','var') || ~isgraphics(f2,'figure')
    f2 = figure(2);
else
    figure(f2);
end
delta_extracted(:,2) = delta_extracted(:,2)*180/pi;
hold on;
h_delta_exp = plot(delta_extracted(:,1), delta_extracted(:,2), ...
    'LineStyle','-', ...
    'LineWidth',1.8, ...
    'Color',[0 0.4470 0.7410], ...
    'DisplayName','Experiment');
uistack(h_delta_exp,'top');
drawnow;
