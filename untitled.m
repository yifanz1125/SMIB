delta = -2*pi:0.01:2*pi;
Kv=0.1;
Xv = 1/Kv;
Xt = Xg+Xv;
Pv = Xt*Vvfm*Ug*sin(delta)./(Xt^2);
figure;
hold on;
plot(delta, Pv,  'b--','linewidth',2);
plot(delta, Xg*Vvfm*Ug*sin(delta)./(Xg^2),  'b-','linewidth',2);


Pgfl = Pv + Ug^2/2*sin(2*delta)*(1/Xg-1/(Xg+Xv));
plot(delta, Pgfl,  'r--','linewidth',2);

Pgf2 = (Pv + Ug^2/2*sin(2*delta)*(1/Xg-1/(Xg+Xv)))*Xg./Ug./(cos(delta)+1e-3);
plot(delta, Pgf2,  'r-','linewidth',2);

axis([-2*pi 2*pi -3 3])