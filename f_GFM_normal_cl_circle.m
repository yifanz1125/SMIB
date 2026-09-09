function dfdt = f_GFM_normal_cl_circle(x)

    Ug =evalin('base','Ug');
    Pm =evalin('base','Pm');
        
    D = evalin('base','D');

    J = evalin('base','J');

    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Lg = Xg/Ws;
    Rg = evalin('base','Rg');

    Vgfm = evalin('base','Vgfm');

    Ilim = evalin('base','Ilim');


    delta = x(1);
    omega = x(2);
    
    
    % 公共项
    Den = Vgfm^2 + Ug^2 - 2*Vgfm*Ug*cos(delta);
    if Den <= 1e-12
        Den = 1e-12;
    end

    Rad = Den/Ilim^2 - Xg^2;
    if Rad < 0
        Rad = 0;
    end

    Sroot = sqrt(Rad);

    Imag = sqrt((Vgfm^2 + Ug^2 - 2*Vgfm*Ug*cos(delta)) / (Xg^2 + Rg^2));


    if Imag <= Ilim 
        P = Rg*(Vgfm^2-Vgfm*Ug*cos(delta))/(Rg^2+Xg^2)+Xg*Vgfm*Ug*sin(delta)/(Rg^2+Xg^2);
    else
        P = Sroot/Den * Ilim^2 * (Vgfm*Ug*cos(delta) - Ug^2) ...
          + Xg/Den * Ilim^2 * Vgfm*Ug*sin(delta) ...
          + Ilim^2*Rg;
    end

    dfdt(1) = omega*Ws;
    dfdt(2) = (Pm-P)/J - D/J*omega;
    
 
    dfdt = dfdt.';

    end