%%
wcc = 500;               % L_Sigma
Xf = 0.05;
Lf = Xf/Ws;
Xsum = Xf+Xg;
Lsum = Xsum/Ws;
kpcc = wcc*Lf;              % k_p_cc
kicc = wcc^2*Lf/4;             % k_i_cc



%%
num1 = [Lsum, (kpcc + Rg), ki, 0];  
den = [
    Lsum^2, ...
    2*Lsum*(kpcc + Rg), ...
    (kpcc + Rg)^2 + 2*Lsum*kicc + (Xsum)^2, ...
    2*(kpcc + Rg)*kicc, ...
    kicc^2
];

G1 = tf(num1, den);

%%
figure;
pzmap(G1);
grid on;

[z, p, k] = zpkdata(G1, 'v')