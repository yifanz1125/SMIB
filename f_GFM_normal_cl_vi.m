function dfdt = f_GFM_normal_cl_vi(x)
% GFM with virtual impedance current limiting
% Xv = 0, if I0 < Ilim
% Xv = Kvi*(|I|-Ilim), if I0 >= Ilim
% |I|^2 = (E^2+Ug^2-2E*Ug*cos(delta))/(Rg^2+(Xg+Xv)^2)

    Ug   = evalin('base','Ug');
    Pm   = evalin('base','Pm');
    Ws   = evalin('base','Ws');
    Xg   = evalin('base','Xg');
    Rg   = evalin('base','Rg');
    E    = evalin('base','Vgfm');
    D    = evalin('base','D');
    J    = evalin('base','J');
    Ilim = evalin('base','Ilim');

    if evalin('base','exist(''Kvi'',''var'')')
        Kvi = evalin('base','Kvi');
    else
        Kvi = 0.2;
    end

    delta = x(1);
    omega = x(2);

    [P, ~] = vi_power(delta,E,Ug,Rg,Xg,Ilim,Kvi);

    dfdt(1) = omega*Ws;
    dfdt(2) = (Pm - P)/J - D/J*omega;

    dfdt = dfdt.';
end

function [P, Xv] = vi_power(delta,E,Ug,Rg,Xg,Ilim,Kvi)

    Den = E^2 + Ug^2 - 2*E*Ug*cos(delta);
    Den = max(Den,1e-12);

    I0 = sqrt(Den/(Rg^2 + Xg^2));

    if I0 < Ilim
        Xv = 0;
    else
        Xv = solve_Xv(Den,Rg,Xg,Ilim,Kvi);
    end

    Xt = Xg + Xv;

    P = Rg/(Rg^2 + Xt^2) * (E^2 - E*Ug*cos(delta)) ...
      + Xt/(Rg^2 + Xt^2) * E*Ug*sin(delta);
end

function Xv = solve_Xv(Den,Rg,Xg,Ilim,Kvi)

    fun = @(z) z - Kvi*(sqrt(Den/(Rg^2 + (Xg+z).^2)) - Ilim);

    lo = 0;
    hi = max(Kvi*(sqrt(Den/(Rg^2 + Xg^2)) - Ilim),1e-6);

    while fun(hi) < 0
        hi = 2*hi;
        if hi > 1e4
            error('solve_Xv:NoBracket','Cannot bracket Xv solution.');
        end
    end

    Xv = fzero(fun,[lo hi]);
    Xv = max(real(Xv),0);
end
