Wbase = 2*pi*50;    % (rad/s)
Vbase = 89.2;   %L-L RMS
Sbase = 1.5e3;
Ibase = Sbase/Vbase;
Zbase = Vbase/Ibase;
Ybase = 1/Zbase;

Lf = 1e-3/Zbase*Wbase
Cf = 5e-6*Zbase*Wbase

Vdc = 200/Vbase

Lg1 = 5.95e-3/Zbase*Wbase
Rg1 = 116e-3/Zbase;  

Lg2 = 3.1e-3/Zbase*Wbase
Rg2 = 90e-3/Zbase

Lg3 = 14.3e-3/Zbase*Wbase
Rg3 = 264e-3/Zbase

Lg1+Lg2

Rg1+Rg2