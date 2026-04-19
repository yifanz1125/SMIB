function dfdt = f_VFM_fault_cl(x)

    Ug =evalin('base','Ug_fault');
    Pin =evalin('base','Pin');
        
    Kpp = evalin('base','Kpp');
    Kip = evalin('base','Kip');

    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Lg = Xg/Ws;
    Rg = evalin('base','Rg');

    Vvfm = evalin('base','Vvfm');
    C_dc = evalin('base','C_dc');

    y_lim = evalin('base','y_lim');

    fault_type = evalin('base','fault_type');


        Ilim = evalin('base','Ilim');
    Phi = evalin('base','Phi');

    if fault_type == "line_cut"
        Xg = evalin('base','Xg_f');
        Rg = evalin('base','Rg_f');
        Lg = Xg/Ws;
    end


    delta = x(1);
    y = x(2);  %voltage
    
    if y>= y_lim 
        y= y_lim;
    end

   
    deltac = acos((Vvfm^2+Ug^2-Ilim^2*(Xg^2+Rg^2))/(2*Vvfm*Ug));
    delta_wrap = mod(delta + pi, 2*pi) - pi;
    if abs(delta_wrap) <= deltac
       P = Rg*(Vvfm^2-Vvfm*Ug*cos(delta))/(Rg^2+Xg^2)+Xg*Vvfm*Ug*sin(delta)/(Rg^2+Xg^2);
    else
       P = Ilim*Ug*cos(delta+Phi) + Ilim^2*Rg;
    end
    

   


    dfdt(1) = Kpp*(Pin-P)+Kip*y;
    if y>= y_lim  && (Pin-P)>=0
        dfdt(2) = 0;  % No change in voltage if below limit
    else
        dfdt(2) = 2/C_dc*(Pin-P);
    end
 
    dfdt = dfdt.';

    end