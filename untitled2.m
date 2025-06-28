% f2=figure(2);
% figure(f2);
hold on;
ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
trange=[0.5,0.5+tc2,0.5+tc2,0.5];
hold on;

fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);