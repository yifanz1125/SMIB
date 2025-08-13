clc

%% Fundamental parameters            
Fs = 1e5;
Ts = 1/Fs;

%% Base values
Wbase = 2*pi*50;    % (rad/s)
Vbase = 1;
Sbase = 1;
Ibase = Sbase/Vbase;
Zbase = Vbase/Ibase;
Ybase = 1/Zbase;

%% AC filter parameters
Lf = 0.05;
Cf = 0.01; %ideal for no Cf 
Lc = 1e-9;


%% Rated line impedance1
switch fault_type
    case "line_cut"
        Lgg = Xgg;
        Rgg = Rgg;
        Lg1 = X1;
        Rg1 = R1;
    otherwise
        Rg11 = Rg;
        Lg11 = Xg;
end


%% Voltage-forming inverter
% Controller
w_droop = 100*2*pi;   %10 1 0
Vdc_ref;
Y_dc;
C_dc;
Kpp;
Kip;
Vvfm;
Pin;

% Current loop
w_i_GFM = 5000*2*pi;

% Voltage loop
w_v_GFM = 600 *2*pi;
Scale_ki_v = 20;

%% fault
% voltage sag1
switch fault_type
    case "voltage_sag"
    %voltage sag
    t_sim_start = 10;
    t0_sag = t_sim_start +t_start;
    dt_sag = t_c;
    v_sag= Ug_fault;
    case "line_cut"
    %line cutting
end




