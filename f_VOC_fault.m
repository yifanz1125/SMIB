function dfdt = f_VOC_fault(x)

    Ug =evalin('base','Ug_fault');
    Pref =evalin('base','Pref');
    kv = evalin('base','kv'); 
    ki = evalin('base','ki'); 
    C = evalin('base','C');
    xi = evalin('base','xi');
    VN = evalin('base','VN');
        


    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Lg = Xg/Ws;
    Rg = evalin('base','Rg');

    VN = evalin('base','VN');
    Qref = evalin('base','Qref2');


    delta = x(1);
    Vgfm = x(2);
    if Vgfm<=1e-6
        Vgfm = 1e-6;
    end
    
    P = Rg*(Vgfm^2-Vgfm*Ug*cos(delta))/(Rg^2+Xg^2)+Xg*Vgfm*Ug*sin(delta)/(Rg^2+Xg^2);
    Q = Xg*(Vgfm^2 - Vgfm*Ug*cos(delta))/(Rg^2+Xg^2) - Rg*Vgfm*Ug*sin(delta)/(Rg^2+Xg^2); 
    

    dfdt(1) = kv*ki/C/Vgfm^2*(Pref - P);  %delta

    if Vgfm<=1e-6 && (xi*(2*VN^2-2*Vgfm^2)*Vgfm/3/kv^2-ki*kv/C/Vgfm*(Q-Qref))<=0
        dfdt(2) = 1e-6;
    else
        dfdt(2) = xi*(2*VN^2-2*Vgfm^2)*Vgfm/3/kv^2-ki*kv/C/Vgfm*(Q-Qref);  %voltage
    end
    dfdt = dfdt.';

    end