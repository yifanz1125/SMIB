function dfdt = f_GFMPI_normal(x)

    Ug =evalin('base','Ug');
    Pm =evalin('base','Pm');
        
    Kpp = evalin('base','Kpp');
    Kip = evalin('base','Kip');

    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Lg = Xg/Ws;
    Rg = evalin('base','Rg');

    Vvfm = evalin('base','Vvfm');


    delta = x(1);
    y = x(2);  %int
    
    P = Rg*(Vvfm^2-Vvfm*Ug*cos(delta))/(Rg^2+Xg^2)+Xg*Vvfm*Ug*sin(delta)/(Rg^2+Xg^2);

    dfdt(1) = Kpp*(Pm-P)+y;
    dfdt(2) = Kip*(Pm-P);
 
    dfdt = dfdt.';

    end