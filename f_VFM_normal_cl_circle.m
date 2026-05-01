function dfdt = f_VFM_normal_cl_circle(x)

    Ug =evalin('base','Ug');
    Pin =evalin('base','Pin');
        
    Kpp = evalin('base','Kpp');
    Kip = evalin('base','Kip');

    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Lg = Xg/Ws;
    Rg = evalin('base','Rg');

    Vvfm = evalin('base','Vvfm');
    C_dc = evalin('base','C_dc');

    Ilim = evalin('base','Ilim');


    delta = x(1);
    y = x(2);  %voltage
    
    
    % 公共项
    Den = Vvfm^2 + Ug^2 - 2*Vvfm*Ug*cos(delta);
    if Den <= 1e-12
        Den = 1e-12;
    end

    Rad = Den/Ilim^2 - Xg^2;
    if Rad < 0
        Rad = 0;
    end

    Sroot = sqrt(Rad);

    Imag = sqrt((Vvfm^2 + Ug^2 - 2*Vvfm*Ug*cos(delta)) / (Xg^2 + Rg^2));


    if Imag <= Ilim 
        P = Rg*(Vvfm^2-Vvfm*Ug*cos(delta))/(Rg^2+Xg^2)+Xg*Vvfm*Ug*sin(delta)/(Rg^2+Xg^2);
    else
        P = Sroot/Den * Ilim^2 * (Vvfm*Ug*cos(delta) - Ug^2) ...
          + Xg/Den * Ilim^2 * Vvfm*Ug*sin(delta) ...
          + Ilim^2*Rg;
    end

    dfdt(1) = Kpp*(Pin-P)+Kip*y;
    dfdt(2) = 2/C_dc*(Pin-P);
    
 
    dfdt = dfdt.';

    end