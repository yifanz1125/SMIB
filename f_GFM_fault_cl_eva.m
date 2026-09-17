function dfdt = f_GFM_fault_cl_eva(x)
% Fault-on GFM model for EVA current limiting during a voltage sag.
% Rg+jXg remains the physical grid impedance, and active power is the
% PCC/POC power. Phase-jump cases have no fault-on interval and therefore
% use f_GFM_normal_cl_eva directly.

    Ug = evalin('base','Ug_fault');
    Pm = evalin('base','Pm');
    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Rg = evalin('base','Rg');
    Xv0 = evalin('base','Xv0');
    Rv0 = evalin('base','Rv0');

    Vgfm = evalin('base','Vgfm');
    D = evalin('base','D');
    J = evalin('base','J');
    Ilim = evalin('base','Ilim');

    delta = x(1);
    omega = x(2);

    P = eva_fault_pcc_power(delta, Ug, Vgfm, Rg, Xg, Rv0, Xv0, Ilim);

    dfdt(1) = omega*Ws;
    dfdt(2) = (Pm-P)/J - D/J*omega;
    dfdt = dfdt.';
end

function P = eva_fault_pcc_power(delta, Ug, Vgfm, Rg, Xg, Rv0, Xv0, Ilim)
    DeltaV2 = Vgfm^2 + Ug^2 - 2*Vgfm*Ug*cos(delta);
    DeltaV2 = max(DeltaV2, 0);

    Rsum0 = Rg + Rv0;
    Xsum0 = Xg + Xv0;
    Zsum20 = Rsum0^2 + Xsum0^2;
    I_unlimited = sqrt(DeltaV2/Zsum20);

    if I_unlimited <= Ilim
        Rv = Rv0;
        Xv = Xv0;
    else
        % Retain the actual angle of the user-defined Rv0+jXv0.
        lambda0 = hypot(Rv0, Xv0);
        if lambda0 <= eps
            error(['EVA requires a nonzero nominal virtual impedance ' ...
                   'Rv0+jXv0 to define its limiting direction.']);
        end
        cos_phi_v = Rv0/lambda0;
        sin_phi_v = Xv0/lambda0;
        projection = Rg*cos_phi_v + Xg*sin_phi_v;
        radicand = projection^2 + DeltaV2/Ilim^2 - (Rg^2 + Xg^2);
        lambda = -projection + sqrt(max(radicand, 0));
        lambda = max(lambda, lambda0);
        Rv = lambda*cos_phi_v;
        Xv = lambda*sin_phi_v;
    end

    Rsum = Rg + Rv;
    Xsum = Xg + Xv;
    Zsum2 = Rsum^2 + Xsum^2;

    % PCC/POC power: P = Ug*real(I) + Rg*abs(I)^2.
    P = (Xsum*Vgfm*Ug*sin(delta) ...
        + Rg*(Vgfm^2 - Vgfm*Ug*cos(delta)) ...
        + Rv*(Vgfm*Ug*cos(delta) - Ug^2)) / Zsum2;
end
