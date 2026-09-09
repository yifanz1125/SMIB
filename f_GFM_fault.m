function dfdt = f_GFM_fault(x)
    Ug =evalin('base','Ug_fault');
    Pm =evalin('base','Pm');
    

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
    D = evalin('base','D');

    J = evalin('base','J');


    delta = x(1);
    omega = x(2);
    
    P = Rg*(Vgfm^2-Vgfm*Ug*cos(delta))/(Rg^2+Xg^2)+Xg*Vgfm*Ug*sin(delta)/(Rg^2+Xg^2);
    
    dfdt(1) = omega*Ws;
    dfdt(2) = (Pm-P)/J - D/J*omega;%(kgfm*(Pm-P)-omega)*w_droop;
 
    dfdt = dfdt.';

    end