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
Lg1 = imag(Z1);
Rg1 = real(Z1);
Lg2 = imag(Z2);
Rg2 = real(Z2);
Lgl = imag(Zl);
Rgl = real(Zl);

%% Grid-forming inverter
% Droop
Pm2 = Pm2;
m_gfm = 1/D_sg;       
w_droop = D_sg/J_sg;  



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




