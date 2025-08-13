function dfdt = f_GFM_normal_rp(x)

    kgfm = evalin('base','kgfm');
    w_droop = evalin('base','w_droop');
    Ug =evalin('base','Ug');
    Pm =evalin('base','Pm');
        
    kp = evalin('base','kp');
    ki = evalin('base','ki');

    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Lg = Xg/Ws;
    Rg = evalin('base','Rg');

    Vgfm = evalin('base','Vgfm');


    delta = x(1);
    
    P = Rg*(Vgfm^2-Vgfm*Ug*cos(delta))/(Rg^2+Xg^2)+Xg*Vgfm*Ug*sin(delta)/(Rg^2+Xg^2);
    eee = 1/w_droop;
    dPddelta = Rg*(Vgfm*Ug*sin(delta))/(Rg^2+Xg^2)+Xg*Vgfm*Ug*cos(delta)/(Rg^2+Xg^2);
    
    omega = (Pm-P)/(1/kgfm-dPddelta*eee);

    dfdt(1) = omega;

 
    dfdt = dfdt.';

    end