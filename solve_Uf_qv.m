function Uf = solve_Uf_qv(delta,Ug,Rg,Xg,UN,kq,Qref)
% solve_Uf_qv
% Solve the algebraic Q-V droop voltage Uf:
%
%   Uf = kq*(Qref - Qe) + UN
%
% with
%
%   Qe = Xg*(Uf^2 - Uf*Ug*cos(delta))/(Rg^2+Xg^2) ...
%      - Rg*Uf*Ug*sin(delta)/(Rg^2+Xg^2)
%
% The equation is quadratic in Uf. The positive-voltage branch is used.

    Z2 = Rg^2 + Xg^2;

    % Equation:
    % (kq*Xg/Z2)*Uf^2 ...
    % + (1 - kq*Xg*Ug*cos(delta)/Z2 + kq*Rg*Ug*sin(delta)/Z2)*Uf ...
    % - (kq*Qref + UN) = 0

    a2 = kq*Xg/Z2;
    a1 = 1 - kq*Xg*Ug*cos(delta)/Z2 + kq*Rg*Ug*sin(delta)/Z2;
    a0 = -(kq*Qref + UN);

    if abs(a2) < 1e-12
        Uf = -a0/a1;
    else
        disc = a1.^2 - 4*a2*a0;
        disc = max(real(disc),0);

        Uf1 = (-a1 + sqrt(disc))/(2*a2);
        Uf2 = (-a1 - sqrt(disc))/(2*a2);

        if Uf1 > 0
            Uf = Uf1;
        else
            Uf = Uf2;
        end
    end

    if ~isfinite(Uf) || Uf <= 0
        error('solve_Uf_qv:InvalidSolution','Uf solution is invalid.');
    end
end
