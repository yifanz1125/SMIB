%%
%clear all
%close all
clc

%% Fundamental parameters
Fs = 100;             % (kHz)
Fs = Fs*1e3;
Ts = 1/Fs;

%% Base values
Wbase = 2*pi*50;    % (rad/s)
Vbase = 1;
Sbase = 1;
Ibase = Sbase/Vbase;
Zbase = Vbase/Ibase;
Ybase = 1/Zbase;



%% Rated line impedance1
Lg = Xg*2;
Rg = Rg*2;
position = 0.5;
Rf=1e-3;

Wmax = w_limit;



%% Grid-following inverter1
% PLL1
w_pll1 = 10 *2*pi;   % (rad/s)
w_tau1 = 1e3 *2*pi;  % (rad/s)
kp_pll1= kp;
ki_pll1= ki;

% Current loop
w_i_GFL1 = 4000 *2*pi;    % (rad/s)



% AC filter parameters
Lf = 0.2;
Cf = 0.01;
Lc = 1e-9;


%% Fault

%voltage sag
t0_sag=1.1;
dt_sag=t_c;
v_sag=Ug_fault;

%phase jump
t0_jump=10;
dt_jump=2;
phase_jump = -1.2;

%line cutting
t0 = 10;
tc = 0.1;


