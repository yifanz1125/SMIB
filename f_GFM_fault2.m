%% grid dynamic
%%
function dfdt = f_GFM_fault2(x)
 
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
    id = x(3);
    iq = x(4);
    
    P = id*Vgfm;
    


    dfdt(1) = omega;
    dfdt(2) = (kgfm*(Pm-P)-omega)*w_droop;

    %grid dynamic
    dfdt(3) = (Vgfm-Ug*cos(delta)+Xg*iq+Lg*omega*iq-Rg*id)/Lg;
    dfdt(4) = (Ug*sin(delta)-Xg*id-Lg*omega*id-Rg*iq)/Lg;


 
    dfdt = dfdt.';

    end