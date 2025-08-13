function dfdt = f_GFM_fault_gy(x)

    kgfm = evalin('base','kgfm');
    w_droop = evalin('base','w_droop');
    Ug =evalin('base','Ug_fault');
    Pm =evalin('base','Pm');
        
    kp = evalin('base','kp');
    ki = evalin('base','ki');

    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Lg = Xg/Ws;
    Rg = evalin('base','Rg');

    Vgfm = evalin('base','Vgfm');

    fault_type = evalin('base','fault_type');

    if fault_type == "line_cut"
        Xg = evalin('base','Xg_f');
        Rg = evalin('base','Rg_f');
        Lg = Xg/Ws;
    end
    


    delta = x(1);
    omega = x(2);
    
    P = Rg*(Vgfm^2-Vgfm*Ug*cos(delta))/(Rg^2+Xg^2)+Xg*Vgfm*Ug*sin(delta)/(Rg^2+Xg^2);
    
    P_e = 2*Xg*Rg*Vgfm*Lg/(Rg^2+Xg^2)^2*omega;

   


    dfdt(1) = omega;
    dfdt(2) = (kgfm*(Pm-P+P_e)-omega)*w_droop;
 
    dfdt = dfdt.';

    end