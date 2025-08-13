%% grid dynamics_ reduced

function dfdt = f_GFM_normal_gy(x)
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
    omega = x(2);
    
    P = Rg*(Vgfm^2-Vgfm*Ug*cos(delta))/(Rg^2+Xg^2)+Xg*Vgfm*Ug*sin(delta)/(Rg^2+Xg^2);
    
    P_e = 2*Xg*Rg*Vgfm*Lg/(Rg^2+Xg^2)^2*omega;

   


    dfdt(1) = omega;
    dfdt(2) = (kgfm*(Pm-P+P_e)-omega)*w_droop;
 
    dfdt = dfdt.';

    end