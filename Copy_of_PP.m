y1 = [0.5;0.45;0.4;0.39;0.3856;0.3848;0.3837;0.38];   
x1 = [55;84;154;189;217;226;250;287];
yq1 = linspace(min(y1), max(y1), 100);  
xq_spline1 = interp1(y1, x1, yq1, 'pchip');


y2 = [0.5;0.45;0.4;0.39;0.3856;0.38];   
x2 = [640;541;389;336;300;287];
yq2 = linspace(min(y1), max(y1), 100);  
xq_spline2 = interp1(y2, x2, yq2, 'pchip');

% y2 = 0:0.01:1;
% delta = asin(y2);
% x2 = cos(delta).*100*pi./y2;
% 
% y3 = [0.1;0.2;0.3;0.4;0.5;0.6;0.8;0.9;0.99;1];   
% x3=[1405;789;558;429;341;274;164;109.1;32.8;0];
% 

% 
% 
% yq3 = linspace(min(y3), max(y3), 100);  
% xq_spline3 = interp1(y3, x3, yq3, 'pchip');

% xq2 = linspace(min(x2), max(x2), 100);  
% yq_spline2 = interp1(x2, y12, xq2, 'spline');
% 
% xq3 = linspace(min(x3), max(x3), 100);  
% yq_spline3 = interp1(x3, y3, xq3, 'spline');

figure;
plot(xq_spline1, yq1, 'k-','LineWidth',1.5);hold on;
plot(x1,y1,'o');

plot(xq_spline2, yq2, 'k-','LineWidth',1.5);hold on;
plot(x2,y2,'o');

xticks(0:100:600);
yticks(0:0.1:1);
xticklabels({'0', '', '200', '','400', '','600'});
yticklabels({'0', '', '0.2', '','0.4', '', '0.6', '','0.8', '', '1.0'});
grid on
axis([0 600 0 1.1])


