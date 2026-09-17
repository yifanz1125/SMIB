clc

%% Fundamental parameters        
f_switching = 20e3;             % (kHz)
Fs = f_switching*50; %f_switching*2e2;
Ts = 1/Fs;
Tc = 1/f_switching;


%% Base values
Wbase = 2*pi*50;    % (rad/s)
Vbase = 89.2;   %L-L RMS
Sbase = 1.5e3;
Ibase = Sbase/Vbase;
Zbase = Vbase/Ibase;
Ybase = 1/Zbase;

%% AC filter parameters
Lf = 0.0592;
Cf = 0.0083; 
Lc = 1e-9;
%%
Vdc_ref = 2.5;


%% Rated line impedance1


%% Grid-forming inverter
E = Vgfm;
Pm = Pm;
m_gfm = 1/D;       
w_droop = D/J;  

% Current loop
w_i_GFM = 1200*2*pi;

% Voltage loop
w_v_GFM = 250 *2*pi;
Scale_ki_v = 20;


f_notch  = 50;   % Hz
BW_notch = 15;   % Hz



%% fault
% voltage sag1
switch fault_type
    case "voltage_sag"
    %voltage sag
    t_sim_start = 2;
    t0_sag = t_sim_start +t_start;
    dt_sag = t_c;%+0.005
    v_sag= Ug_fault;
    case "line_cut"
    %line cutting
end




