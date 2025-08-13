clc

%% Fundamental parameters
f_switching = 10e3;             % (kHz)
Fs = f_switching*2e2;
Ts = 1/Fs;
Tc = 1/f_switching;

%% Base values
Wbase = 2*pi*50;    % (rad/s)
Vbase = 1;
Sbase = 1;
Ibase = Sbase/Vbase;
Zbase = Vbase/Ibase;
Ybase = 1/Zbase;



%% Rated line impedance1
Lgg = Xgg;
Rgg = Rgg;
Lg1 = X1;
Rg1 = R1;
Rg11 = Rg;
Lg11 = Xg;




%% Grid-following inverter1
% PLL1
w_pll1 = 10 *2*pi;   % (rad/s)
w_tau1 = 1000 *2*pi;  % (rad/s)
kp_pll1= kp;
ki_pll1= ki;

% Current loop
w_i_GFL1 = 1000 *2*pi;    % (rad/s)



% AC filter parameters
Lf = 0.05;
Cf1 = 0.01;%0.02
Cf2 = 0.01;%0.02
Lc = 1e-9;

Vdc = 2.5;
C_dc = 2.25;

Id0 = Id;
Iq01 = (Ug*cos(asin((Xg*Id+Rg*Iq)/Ug))+Id*Rg-Iq*Xg)*Cf1+Iq;
Iq02 = (Ug*cos(asin((Xg*Id+Rg*Iq)/Ug))+Id*Rg-Iq*Xg)*Cf2+Iq;

%% Fault
switch fault_type
    case "voltage_sag"
    %voltage sag
    t_sim_start = 0.5;
    t0_sag=t_sim_start +t_start;
    t0=10;
    dt_sag=t_c;
    v_sag=Ug_fault;
    tc = t_c;
    Rf = 1e-3;
    %phase jump
    t0_jump=10;
    dt_jump=2;
    phase_jump = -1.2;
    case "line_cut"
    t0_sag=20;
    t0_jump=20;
    dt_sag=t_c;
    v_sag=Ug_fault;
    %line cutting
    t_sim_start = 0.5;
    t0 = t_sim_start +t_start;
    tc = t_c;
end


