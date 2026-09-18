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
% Xg = 0.307; 
% Rg = 0.0338;
% Xv0 = 0.124;
% Rv0 = 0.0162;

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

%% EVA capacitor-voltage feedback: second-order LPF

EVA_LPF2_f     = 30;              % -3 dB bandwidth (Hz)
EVA_LPF2_w     = 2*pi*EVA_LPF2_f; % rad/s
EVA_LPF2_zeta  = 1/sqrt(2);       % Butterworth damping ratio
EVA_LPF2_Ts    = Tc;              % change to Tc later if required

% Prewarped Tustin transformation
EVA_LPF2_k  = tan(EVA_LPF2_w*EVA_LPF2_Ts/2);
EVA_LPF2_a0 = 1 ...
            + 2*EVA_LPF2_zeta*EVA_LPF2_k ...
            + EVA_LPF2_k^2;

% Coefficients in descending powers of z
EVA_LPF2_num = [EVA_LPF2_k^2, ...
                2*EVA_LPF2_k^2, ...
                EVA_LPF2_k^2] / EVA_LPF2_a0;

EVA_LPF2_den = [1, ...
                2*(EVA_LPF2_k^2-1)/EVA_LPF2_a0, ...
                (1-2*EVA_LPF2_zeta*EVA_LPF2_k ...
                +EVA_LPF2_k^2)/EVA_LPF2_a0];



%% fault
% voltage sag1
switch fault_type
    case "voltage_sag"
    %voltage sag
    t_sim_start = 2;
    t0_sag = t_sim_start +t_start;
    t0_jump = t_sim_start + t_end + 1;
    dt_sag = t_c;%+0.005
    v_sag= Ug_fault;
    case "line_cut"
    %line cutting
    case "phase_jump"
    t_sim_start = 2;
    t0_sag = t_sim_start + t_end + 1;
    t0_jump = t_sim_start +t_start;
    jump_value = -delta_jump_initial;  
end




