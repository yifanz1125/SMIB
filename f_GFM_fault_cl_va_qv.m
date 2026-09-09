function dfdt = f_GFM_fault_cl_va_qv(x)
% Fault-on GFM + virtual admittance current limiting + Q-V droop

    Ug   = evalin('base','Ug_fault');
    Pm   = evalin('base','Pm');
    Ws   = evalin('base','Ws');
    Xg   = evalin('base','Xg');
    Rg   = evalin('base','Rg');
    UN   = evalin('base','Vgfm');
    D    = evalin('base','D');
    J    = evalin('base','J');
    Ilim = evalin('base','Ilim');

    fault_type = evalin('base','fault_type');
    if fault_type == "line_cut"
        Xg = evalin('base','Xg_f');
        Rg = evalin('base','Rg_f');
    end

    if evalin('base','exist(''kq'',''var'')')
        kq = evalin('base','kq');
    else
        kq = 1;
    end

    if evalin('base','exist(''Qref1'',''var'')')
        Qref = evalin('base','Qref1');
    else
        Ug_nom = evalin('base','Ug');
        Xg_nom = evalin('base','Xg');
        Rg_nom = evalin('base','Rg');
        Qref = get_Qref_va_qv(UN,Ug_nom,Rg_nom,Xg_nom,Pm);
    end

    delta = x(1);
    omega = x(2);
    [P, ~, ~, ~] = va_qv_power(delta,Ug,Rg,Xg,UN,Ilim,kq,Qref);

    dfdt(1) = omega*Ws;
    dfdt(2) = (Pm - P)/J - D/J*omega;
    dfdt = dfdt.';
end

function Qref = get_Qref_va_qv(UN,Ug,Rg,Xg,Pm)
    fsep = @(delta) Rg*(UN^2 - UN*Ug*cos(delta))/(Rg^2 + Xg^2) ...
                  + Xg*UN*Ug*sin(delta)/(Rg^2 + Xg^2) - Pm;
    opts = optimoptions('fsolve','Display','off');
    deltas = fsolve(fsep,0,opts);
    Qref = Xg*(UN^2 - UN*Ug*cos(deltas))/(Rg^2+Xg^2) ...
         - Rg*UN*Ug*sin(deltas)/(Rg^2+Xg^2);
end

function [P,Qe,Uf,Xv] = va_qv_power(delta,Ug,Rg,Xg,UN,Ilim,kq,Qref)
    Uf = solve_Uf_qv(delta,Ug,Rg,Xg,UN,kq,Qref);
    Den0 = Uf^2 + Ug^2 - 2*Uf*Ug*cos(delta);
    Den0 = max(Den0,1e-12);
    I0 = sqrt(Den0/(Rg^2 + Xg^2));

    if I0 <= Ilim
        Xv = 0; Xt = Xg;
        P = Rg/(Rg^2 + Xt^2)*(Uf^2 - Uf*Ug*cos(delta)) ...
          + Xt/(Rg^2 + Xt^2)*Uf*Ug*sin(delta);
        Qe = Xt/(Rg^2 + Xt^2)*(Uf^2 - Uf*Ug*cos(delta)) ...
           - Rg/(Rg^2 + Xt^2)*Uf*Ug*sin(delta);
    else
        Xv = sqrt(Den0/Ilim^2 - Rg^2) - Xg;
        Xv = max(real(Xv),0); Xt = Xg + Xv;
        P = Rg/Den0*Ilim^2*(Uf^2 - Uf*Ug*cos(delta)) ...
          + Xt/Den0*Ilim^2*Uf*Ug*sin(delta);
        Qe = Xt/Den0*Ilim^2*(Uf^2 - Uf*Ug*cos(delta)) ...
           - Rg/Den0*Ilim^2*Uf*Ug*sin(delta);
    end
end

function Uf = solve_Uf_qv(delta,Ug,Rg,Xg,UN,kq,Qref)
    Z2 = Rg^2 + Xg^2;
    a = 1 - Xg*kq*Ug*cos(delta)/Z2 + Rg*kq*Ug*sin(delta)/Z2;
    b = kq*Xg/Z2;
    rhs = kq*Qref + UN;
    if abs(b) < 1e-12
        Uf = rhs/a;
    else
        disc = a^2 + 4*b*rhs;
        disc = max(real(disc),0);
        Uf = (-a + sqrt(disc))/(2*b);
    end
    if ~isfinite(Uf) || Uf <= 0
        fun = @(U) U - (kq*(Qref - (Xg*(U.^2-U*Ug*cos(delta))/Z2 ...
                    - Rg*U*Ug*sin(delta)/Z2)) + UN);
        lo = 1e-6; hi = max(2*UN,2);
        while fun(lo)*fun(hi) > 0
            hi = 2*hi;
            if hi > 1e4
                error('solve_Uf_qv:NoBracket','Cannot bracket Uf solution.');
            end
        end
        Uf = fzero(fun,[lo hi]);
    end
end
