function dfdt = f_GFL2(x)

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

    delta = x(1);
    Int = x(2);




    if model == "original"
        Vq = (Xg*Id+Rg*Iq-Ug*sin(delta)+Id*Lg*Int)/(1-Id*Lg*kp);
    
        dfdt(1) = kp*Vq+Int - W_g; %ddelta
        dfdt(2) = ki*Vq; %dint

    elseif model == "no_grid_dy" 
        Vq = Xg*Id+Rg*Iq-Ug*sin(delta);
        dfdt(1) = kp*Vq+Int; %ddelta;
        dfdt(2) = ki*Vq;

    else 
        dfdt(1) = 0;
        dfdt(2) = 0;
    end
    dfdt = dfdt.';

    end