function [P,Qe,Uf,Re,Imag,limited,saturated] = circle_qv_power(delta,Ug,Rg,Xg,UN,Ilim,kq,Qref)
% Circular current limiter + Q-V droop with Uf >= 0 saturation.
%
% If the limited-region Q-V algebraic equation has no positive real root,
% Uf is clamped to Uf_min = 0. Electrical quantities are then recomputed
% using the actual current at Uf = 0.
%
% Outputs:
%   P, Qe       active/reactive power
%   Uf          internal voltage magnitude
%   Re          added virtual resistance
%   Imag        actual current magnitude
%   limited     current limiter active flag
%   saturated   Uf lower saturation active flag

    Uf_min = 0;
    saturated = false;

    %% ===== 1. Unlimited Q-V candidate =====
    Uf0 = solve_Uf_qv(delta,Ug,Rg,Xg,UN,kq,Qref);
    Uf0 = max(Uf0,Uf_min);

    Den0 = max(Uf0^2 + Ug^2 - 2*Uf0*Ug*cos(delta),1e-12);
    I0 = sqrt(Den0/(Rg^2 + Xg^2));

    if I0 <= Ilim
        Uf = Uf0;
        limited = false;
        Re = 0;
        Imag = I0;

        P = Rg/(Rg^2 + Xg^2)*(Uf^2 - Uf*Ug*cos(delta)) ...
          + Xg/(Rg^2 + Xg^2)*Uf*Ug*sin(delta);

        Qe = Xg/(Rg^2 + Xg^2)*(Uf^2 - Uf*Ug*cos(delta)) ...
           - Rg/(Rg^2 + Xg^2)*Uf*Ug*sin(delta);
        return;
    end

    %% ===== 2. Search limited-region positive Uf root =====
    limited = true;

    Umin_search = max(Uf_min,1e-8);
    Umax = max([4*UN,4*Uf0,4]);

    bracket = find_valid_bracket(...
        Umin_search,Umax,2000,delta,Ug,Rg,Xg,UN,Ilim,kq,Qref);

    if isempty(bracket)
        Umax = 20*max([UN,Uf0,1]);
        bracket = find_valid_bracket(...
            Umin_search,Umax,5000,delta,Ug,Rg,Xg,UN,Ilim,kq,Qref);
    end

    if isempty(bracket)
        %% ===== 3. Uf lower saturation =====
        Uf = Uf_min;
        saturated = true;

        % Re-evaluate the actual mode at the saturated voltage.
        Den = max(Uf^2 + Ug^2 - 2*Uf*Ug*cos(delta),1e-12);
        I_unlimited = sqrt(Den/(Rg^2 + Xg^2));

        if I_unlimited <= Ilim
            limited = false;
            Re = 0;
            Imag = I_unlimited;

            P = Rg/(Rg^2 + Xg^2)*(Uf^2 - Uf*Ug*cos(delta)) ...
              + Xg/(Rg^2 + Xg^2)*Uf*Ug*sin(delta);

            Qe = Xg/(Rg^2 + Xg^2)*(Uf^2 - Uf*Ug*cos(delta)) ...
               - Rg/(Rg^2 + Xg^2)*Uf*Ug*sin(delta);
        else
            limited = true;
            Rt = sqrt(max(Den/Ilim^2 - Xg^2,0));
            Re = max(Rt-Rg,0);
            Imag = Ilim;

            P = Rt/Den*Ilim^2*(Uf*Ug*cos(delta)-Ug^2) ...
              + Xg/Den*Ilim^2*Uf*Ug*sin(delta) ...
              + Ilim^2*Rg;

            Qe = Xg/Den*Ilim^2*(Uf^2-Uf*Ug*cos(delta)) ...
               - Rt/Den*Ilim^2*Uf*Ug*sin(delta);
        end
        return;
    end

    %% ===== 4. Positive limited-region root =====
    if bracket(1) == bracket(2)
        Uf = bracket(1);
    else
        Uf = fzero(@(U) limited_residual(...
            U,delta,Ug,Rg,Xg,UN,Ilim,kq,Qref),bracket);
    end

    Uf = max(Uf,Uf_min);

    Den = max(Uf^2 + Ug^2 - 2*Uf*Ug*cos(delta),1e-12);
    Rt = sqrt(max(Den/Ilim^2 - Xg^2,0));

    Re = max(Rt-Rg,0);
    Imag = Ilim;

    P = Rt/Den*Ilim^2*(Uf*Ug*cos(delta)-Ug^2) ...
      + Xg/Den*Ilim^2*Uf*Ug*sin(delta) ...
      + Ilim^2*Rg;

    Qe = Xg/Den*Ilim^2*(Uf^2-Uf*Ug*cos(delta)) ...
       - Rt/Den*Ilim^2*Uf*Ug*sin(delta);
end

function bracket = find_valid_bracket(...
    Umin,Umax,N,delta,Ug,Rg,Xg,UN,Ilim,kq,Qref)

    Ugrid = linspace(Umin,Umax,N);
    Fgrid = NaN(size(Ugrid));

    for kk = 1:N
        Fgrid(kk) = limited_residual(...
            Ugrid(kk),delta,Ug,Rg,Xg,UN,Ilim,kq,Qref);
    end

    bracket = [];

    for kk = 1:N-1
        if ~isfinite(Fgrid(kk)) || ~isfinite(Fgrid(kk+1))
            continue;
        end

        if Fgrid(kk) == 0
            bracket = [Ugrid(kk),Ugrid(kk)];
            return;
        end

        if Fgrid(kk)*Fgrid(kk+1) < 0
            bracket = [Ugrid(kk),Ugrid(kk+1)];
            return;
        end
    end
end

function F = limited_residual(Uf,delta,Ug,Rg,Xg,UN,Ilim,kq,Qref)

    if ~isfinite(Uf) || Uf < 0
        F = NaN;
        return;
    end

    Den = max(Uf^2 + Ug^2 - 2*Uf*Ug*cos(delta),1e-12);
    Rt2 = Den/Ilim^2 - Xg^2;

    if Rt2 < 0
        F = NaN;
        return;
    end

    Rt = sqrt(Rt2);

    Qe = Xg/Den*Ilim^2*(Uf^2-Uf*Ug*cos(delta)) ...
       - Rt/Den*Ilim^2*Uf*Ug*sin(delta);

    F = Uf - (kq*(Qref-Qe)+UN);
end
