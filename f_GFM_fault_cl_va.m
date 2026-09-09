function dfdt = f_GFM_fault_cl_va(x)

    Ug =evalin('base','Ug_fault');
    Pm =evalin('base','Pm');
    

    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Lg = Xg/Ws;
    Rg = evalin('base','Rg');

    Vgfm = evalin('base','Vgfm');

    Ilim = evalin('base','Ilim');


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
    
   % 公共项
    Den = Vgfm^2 + Ug^2 - 2*Vgfm*Ug*cos(delta);
    if Den <= 1e-12
        Den = 1e-12;
    end

    Rad = Den/Ilim^2 - Rg^2;
    if Rad < 0
        Rad = 0;
    end

    Xvar = sqrt(Rad);

    Imag = sqrt((Vgfm^2 + Ug^2 - 2*Vgfm*Ug*cos(delta)) / (Xg^2 + Rg^2));


    if Imag <= Ilim 
        P = Rg*(Vgfm^2-Vgfm*Ug*cos(delta))/(Rg^2+Xg^2)+Xg*Vgfm*Ug*sin(delta)/(Rg^2+Xg^2);
    else
        P = Rg / Den * Ilim^2 * (Vgfm^2 - Vgfm*Ug*cos(delta))+ Xvar / Den * Ilim^2 * Vgfm*Ug*sin(delta);
    end

    dfdt(1) = omega*Ws;
    dfdt(2) = (Pm-P)/J - D/J*omega;%(kgfm*(Pm-P)-omega)*w_droop
 
    dfdt = dfdt.';

    end