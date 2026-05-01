function dfdt = f_VOC_normal_cl_circle(x)

    Ug =evalin('base','Ug');
    Pref =evalin('base','Pref');
    kv = evalin('base','kv'); 
    ki = evalin('base','ki'); 
    C = evalin('base','C');
    xi = evalin('base','xi');
    Ilim = evalin('base','Ilim');
        


    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Lg = Xg/Ws;
    Rg = evalin('base','Rg');

    VN = evalin('base','VN');
    Qref = evalin('base','Qref2');


    delta = x(1);
    Vgfm = x(2);
    if Vgfm<=1e-3
        Vgfm = 1e-3;
    end


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
    
    Q = Xg*Ilim^2/Den * (Vgfm^2 - Vgfm*Ug*cos(delta)) ...
      - Ilim^2*Sroot/Den * Vgfm*Ug*sin(delta);

    P = Sroot/Den * Ilim^2 * (Vgfm*Ug*cos(delta) - Ug^2) ...
      + Xg/Den * Ilim^2 * Vgfm*Ug*sin(delta) ...
      + Ilim^2*Rg;
    

    dfdt(1) = kv*ki/C/Vgfm^2*(Pref - P);  %delta

    if Vgfm<=1e-6 && (xi*(2*VN^2-2*Vgfm^2)*Vgfm/3/kv^2-ki*kv/C/Vgfm*(Q-Qref))<=0
        dfdt(2) = 1e-6;
    else
        dfdt(2) = xi*(2*VN^2-2*Vgfm^2)*Vgfm/3/kv^2-ki*kv/C/Vgfm*(Q-Qref);  %voltage
    end
    dfdt = dfdt.';

    end