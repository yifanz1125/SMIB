function dfdt = f_GFL_limit(x)

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

    Wmax = evalin('base','w_max');
    Wmin = evalin('base','w_min');

    delta = x(1);
    omega = x(2);

    M = (1-kp*Lg*Id)/ki;
    D = kp/ki*Ug*cos(delta)-Id*Lg;


    if omega>=Wmax
        dfdt(1) = Wmax; %ddelta
        dfdt(2) = (Xg*Id+Rg*Iq-Ug*sin(delta) - kp/ki*Ug*cos(delta)*Wmax)*ki; %domega
    elseif omega<=Wmin
        dfdt(1) = Wmin; %ddelta
        dfdt(2) = (Xg*Id+Rg*Iq-Ug*sin(delta) - kp/ki*Ug*cos(delta)*Wmin)*ki; %domega
    else
        dfdt(1) = omega; %ddelta
        dfdt(2) = (Xg*Id+Rg*Iq-Ug*sin(delta) - D*omega)/M; %domega
    end


    dfdt = dfdt.';

    end