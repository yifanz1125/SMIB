y1 = [0.35;0.37;0.38;0.385;0.389;0.39;0.4;0.5;0.6;0.7;0.8];   
x1 = [0;205;287.2;317;330;337;390;640;816;957;1130];

y2 = [0.385;0.389;0.39;0.4;0.5;0.6;0.7;0.8;0.9;1];
x2 = [317;193;188;153;55;25.8;12.9;6;2;0];
% 
% y3 = [0.1;0.2;0.3;0.4;0.5;0.6;0.8;0.9;0.99;1];   
% x3=[1405;789;558;429;341;274;164;109.1;32.8;0];

yq1 = linspace(min(y1), max(y1), 100);  
xq_spline1 = interp1(y1, x1, yq1, 'pchip');


% yq3 = linspace(min(y3), max(y3), 100);  
% xq_spline3 = interp1(y3, x3, yq3, 'pchip');

yq2 = linspace(min(y2), max(y2), 100);  
xq_spline2 = interp1(y2, x2, yq2, 'pchip');
% 
% xq3 = linspace(min(x3), max(x3), 100);  
% yq_spline3 = interp1(x3, y3, xq3, 'spline');

figure;
plot(xq_spline1, yq1, 'k-','LineWidth',1.5);hold on;
plot(xq_spline2, yq2, 'k-','LineWidth',1.5);hold on;
% % plot(x1, y1, 'ko','LineWidth',1.5);hold on;
% plot(x2, y2, 'k-', 'LineWidth',1.5);
% % plot(xq3, yq_spline3, 'k:', 'DisplayName', 'Hopf Bifurcation','LineWidth',1.5);
% plot(xq_spline3, yq3, 'k-','LineWidth',1.5);hold on;


xticks(0:100:700);
yticks(0:0.1:1);
xticklabels({'0', '', '200', '','400', '','600'});
yticklabels({'0', '', '0.2', '','0.4', '', '0.6', '','0.8', '', '1.0'});
grid on
axis([0 700 0 1.1])


