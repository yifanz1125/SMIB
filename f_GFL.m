function dfdt = f_GFL(x)

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
    omega = x(2);

    M = (1-kp*Lg*Id)/ki;
    D = kp/ki*Ug*cos(delta)-Id*Lg;


    if model == "original"
    
        dfdt(1) = omega; %ddelta
        dfdt(2) = (Xg*Id+Rg*Iq-Ug*sin(delta) - D*omega+Lg*Id*W_g)/M; %domega

    elseif model == "no_grid_dy" % ignore w2 & w1   
        dfdt(1) = omega;
        dfdt(2) = (Xg*Id+Rg*Iq-Ug*sin(delta) - kp/ki*Ug*cos(delta)*omega)*ki;

    else 
        dfdt(1) = 0;
        dfdt(2) = 0;
    end
    dfdt = dfdt.';

    end