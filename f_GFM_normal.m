function dfdt = f_GFM_normal(x)
    global limit_type

    Ug =evalin('base','Ug');
    Pm =evalin('base','Pm');
    
    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Lg = Xg/Ws;
    Rg = evalin('base','Rg');

    Vgfm = evalin('base','Vgfm');

    D = evalin('base','D');

    J = evalin('base','J');


    delta = x(1);
    omega = x(2);

    if ~isempty(limit_type) && limit_type == "EVA"
        % EVA: Rg+jXg is the physical grid impedance and Rv0+jXv0 is
        % present even before current limiting. P is measured at PCC/POC.
        Xv0 = evalin('base','Xv0');
        Rv0 = evalin('base','Rv0');
        Rsum = Rg + Rv0;
        Xsum = Xg + Xv0;
        Zsum2 = Rsum^2 + Xsum^2;

        P = (Xsum*Vgfm*Ug*sin(delta) ...
            + Rg*(Vgfm^2 - Vgfm*Ug*cos(delta)) ...
            + Rv0*(Vgfm*Ug*cos(delta) - Ug^2)) / Zsum2;
    else
        % Original definition retained for all existing modes. In the old
        % VA case, Xg already represents grid plus nominal virtual reactance.
        P = Rg*(Vgfm^2-Vgfm*Ug*cos(delta))/(Rg^2+Xg^2) ...
          + Xg*Vgfm*Ug*sin(delta)/(Rg^2+Xg^2);
    end
    


    dfdt(1) = omega*Ws;
    dfdt(2) = (Pm-P)/J - D/J*omega;%(kgfm*(Pm-P)-omega)*w_droop;
 
    dfdt = dfdt.';

    end