function dfdt = f_GFL_prefault2(x)

    model = evalin('base','model');
    Id = evalin('base','Id');
    Iq = evalin('base','Iq');
    Ug =evalin('base','Ug');
        
    kp = evalin('base','kp');
    ki = evalin('base','ki');

    Ws = evalin('base','Ws');
    W_g = evalin('base','W_g');

    Xg = evalin('base','Xg');
    
    Lg = Xg/Ws;
    Rg = evalin('base','Rg');
    
    fault_type = evalin('base','fault_type');

    if fault_type == "line_cut"
        Xg0 = evalin('base','Xg0');
        Rg0 = evalin('base','Rg0');
        Xg = Xg0;
        Rg = Rg0;
        Lg = Xg/Ws;
    end



    delta = x(1);
    Int = x(2);


    if model == "original"

        Vq = (Xg*Id+Rg*Iq-Ug*sin(delta)+Id*Lg*Int)/(1-Id*Lg*kp);
    
        dfdt(1) = kp*Vq+Int - W_g; %ddelta
        dfdt(2) = ki*Vq; %dint

    elseif model == "no_grid_dy" % ignore w2 & w1   
        Vq = Xg*Id+Rg*Iq-Ug*sin(delta);
        
        dfdt(1) = kp*Vq+Int; %ddelta;
        dfdt(2) = ki*Vq;

    else 
        dfdt(1) = 0;
        dfdt(2) = 0;
    end
    dfdt = dfdt.';

    end