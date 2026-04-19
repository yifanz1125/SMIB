function dfdt = f_VFM_normal_cl(x)

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
    Phi = evalin('base','Phi');


    delta = x(1);
    y = x(2);  %voltage
    
    
    deltac = acos((Vvfm^2+Ug^2-Ilim^2*(Xg^2+Rg^2))/(2*Vvfm*Ug));
    delta_wrap = mod(delta + pi, 2*pi) - pi;
    if abs(delta_wrap) <= deltac
       P = Rg*(Vvfm^2-Vvfm*Ug*cos(delta))/(Rg^2+Xg^2)+Xg*Vvfm*Ug*sin(delta)/(Rg^2+Xg^2);
    else
       P = Ilim*Ug*cos(delta+Phi) + Ilim^2*Rg;
    end

    dfdt(1) = Kpp*(Pin-P)+Kip*y;
    dfdt(2) = 2/C_dc*(Pin-P);
    
 
    dfdt = dfdt.';

    end