function dfdt = f_VFM_normal(x)

    Ug =evalin('base','Ug');
    Pin =evalin('base','Pin');
        
    Kpp = evalin('base','Kpp');
    Kip = evalin('base','Kip');

    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Lg = Xg/Ws;
    Rg = evalin('base','Rg');

    Vvfm = evalin('base','Vvfm');
    C_dc = evalin('base','C_dc');


    delta = x(1);
    y = x(2);  %voltage
    
    P = Rg*(Vvfm^2-Vvfm*Ug*cos(delta))/(Rg^2+Xg^2)+Xg*Vvfm*Ug*sin(delta)/(Rg^2+Xg^2);
    

   


    dfdt(1) = Kpp*(Pin-P)+Kip*y;
    dfdt(2) = 2/C_dc*(Pin-P);
 
    dfdt = dfdt.';

    end