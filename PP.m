y1 = [0.01;0.05;0.1;0.2;0.3;0.4;0.45;0.5;0.6;0.7;0.8;0.9;0.95;1];   
x1 = [7740;1070;409;140;69;40;30;24;14.3;8.4;4.43;1.8;0.7;0];

y2 = 0:0.01:1;
delta = asin(y2);
x2 = cos(delta).*100*pi./y2;

y3 = [0.1;0.2;0.3;0.4;0.5;0.6;0.8;0.9;0.99;1];   
x3=[1405;789;558;429;341;274;164;109.1;32.8;0];

yq1 = linspace(min(y1), max(y1), 100);  
xq_spline1 = interp1(y1, x1, yq1, 'pchip');


yq3 = linspace(min(y3), max(y3), 100);  
xq_spline3 = interp1(y3, x3, yq3, 'pchip');

% xq2 = linspace(min(x2), max(x2), 100);  
% yq_spline2 = interp1(x2, y12, xq2, 'spline');
% 
% xq3 = linspace(min(x3), max(x3), 100);  
% yq_spline3 = interp1(x3, y3, xq3, 'spline');

figure;
plot(xq_spline1, yq1, 'k-','LineWidth',1.5);hold on;
% plot(x1, y1, 'ko','LineWidth',1.5);hold on;
plot(x2, y2, 'k-', 'LineWidth',1.5);
% plot(xq3, yq_spline3, 'k:', 'DisplayName', 'Hopf Bifurcation','LineWidth',1.5);
plot(xq_spline3, yq3, 'k-','LineWidth',1.5);hold on;


xticks(0:100:600);
yticks(0:0.1:1);
xticklabels({'0', '', '200', '','400', '','600'});
yticklabels({'0', '', '0.2', '','0.4', '', '0.6', '','0.8', '', '1.0'});
grid on
axis([0 600 0 1.1])


