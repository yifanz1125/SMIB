%% parameter

%grid
Rg = 0.01;
Xg = 1/2.1;
Ug = 1;
Ws = 2*pi*50; 
Lg= Xg/Ws;
W_g = 0;


%GFL
Id = 1;
Iq = 0; %negative q
kp = 10*2*pi;
ki = kp*100;
w_max = Ws*100;
w_min = -w_max;

%system
global system;
global fault_type; %line_cut voltage_sag frequency
fault_type = "voltage_sag"; %"voltage_sag";%"line_cut";%"line_cut";
system = "GFL";
model = "original";% "original"

switch fault_type
    case "voltage_sag"
        %fault sag
        Ug_fault = 0.2;
        t_c = 0.01;
        X1 = 0.1;
        R1 = 0.01;
        Xgg = (Xg - X1)*2;
        Rgg = (Rg - R1)*2;
    case "line_cut"
    %fault line cut 
        t_c = 0.06;
        X1 = 0.1;
        R1 = 1e-3;
        Xgg = Xg-X1;
        Rgg = Rg - R1;
        Xg0 = Xgg/2+X1;
        Rg0 = Rgg/2+R1;
        position = 1; %fault to inf bus
        Rf = 1e-5/(690^2/15e6);
        Im_temp = Rf*(Xgg*position*1j+Rgg*position)/(Rf + Xgg*position*1j+Rgg*position)+Xgg*(1-position)*1j+Rgg*(1-position);
        Imgf = Im_temp *(Xgg*1j+Rgg)/(Xgg*1j+Rgg+Im_temp); %(Rf//Xg*location+Xg*(1-location))//Xg
        Xg_f = imag(Imgf)+X1;
        Lg_f = Xg_f/Ws;
        Rg_f = real(Imgf)+R1;
        Iop_tmp = 1/( ((2-position)*(Xgg*1j+Rgg)*(position)*(Xgg*1j+Rgg))/(2*Xgg*1j+2*Rgg) + Rf); %Xg*(2-location)//Xg*(location)+Rf
        Ug_fault = Iop_tmp*Rf+Iop_tmp/2*position*(1-position)*(Xgg*1j+Rgg);
        Ug_fault_angle = angle(Ug_fault);
        Ug_fault= abs(Ug_fault);
end



%% 
try
    system;
catch
    system = "GFL";
end
%% SEP
for ki = ki
x=(0:0.1:1)*2*pi;
n = length(x);
x_set = zeros(2,n);
x_set(1,:) = x;

torralence = 1e-2; 
mm = 1;
ep_set = [];
options = optimoptions('fsolve','FunctionTolerance',1e-10,'MaxIterations',100000,'OptimalityTolerance',1e-10);
for n = 1:length(x_set(1,:))
    xep = x_set(:,n);
    [xep,ferr,~,~,A] = fsolve(@f,xep,options);
    
    if maxabs(ferr) < torralence
        if isnewxep(ep_set,xep,torralence)
           
            [V,Lambda]=eig(A);
            Lambda = diag(Lambda);
            sig = sign(sign(real(Lambda))+0.1); % zero counted as positive
            sig = (sig + 1)/2;                  % [0,1]
            flag = sum(sig);                    % number of non-negative eigenvalues

            v = V(:,~sig);                      % the stable sub-space
            vv = V(:,~(~sig)); 
            
            ep_set(mm).xep = xep; %#ok<*SAGROW> 
            ep_set(mm).A = A;
            ep_set(mm).Lambda = Lambda;
            ep_set(mm).V = V;   
            ep_set(mm).v = v;     % stable eigenvectors of unstable ep 
            ep_set(mm).vv = vv;
            ep_set(mm).flag = flag;
           
            mm = mm+1;
            if flag == 0
               postfault_SEP=xep;
            end
        end
    end
end

%prefault
clear xep flag v V Lambda A sig m
ep_set_pre = [];
m = 1;
for n = 1:length(x_set(1,:))
    xep = x_set(:,n);
    [xep,ferr,~,~,A] = fsolve(@(x)f_prefault(0,x),xep,options);
    if maxabs(ferr) < torralence
        if isnewxep(ep_set_pre,xep,torralence)
           
            [V,Lambda]=eig(A);
            Lambda = diag(Lambda);
            sig = sign(sign(real(Lambda))+0.1); % zero counted as positive
            sig = (sig + 1)/2;                  % [0,1]
            flag = sum(sig);                    % number of non-negative eigenvalues

            v = V(:,~sig);                      % the stable sub-space
            
            ep_set_pre(m).xep = real(xep); %#ok<*SAGROW> 
            ep_set_pre(m).A = A;
            ep_set_pre(m).Lambda = Lambda;
            ep_set_pre(m).V = V;   
            ep_set_pre(m).v = v;     % stable eigenvectors of unstable ep 
            ep_set_pre(m).flag = flag;
           
            m = m+1;
            if flag == 0
            jacob=A;
            prefault_SEP=xep;
            end
        end
    end
end

for mm = 1:length(ep_set)
    disp_v('Index',mm);
    disp_v('Equilibrium',ep_set(mm).xep);
    disp_v('Eigenvalue', ep_set(mm).Lambda);
    disp_v('Eigenvector',ep_set(mm).V);
end

clear ep_set_ext;
for n = 1:length(ep_set)
    mm = (n-1)*10;
    ep_set_ext(mm+1)=ep_set(n); %#ok<*AGROW> 
    ep_set_ext(mm+2)=ep_set(n);
    ep_set_ext(mm+3)=ep_set(n); 
    ep_set_ext(mm+4)=ep_set(n); 
    ep_set_ext(mm+5)=ep_set(n);
    ep_set_ext(mm+6)=ep_set(n);
    ep_set_ext(mm+7)=ep_set(n);
    ep_set_ext(mm+8)=ep_set(n);
    ep_set_ext(mm+9)=ep_set(n);
    ep_set_ext(mm+10)=ep_set(n);
    ep_set_ext(mm+2).xep(1) = ep_set(n).xep(1) - 2*pi;
    ep_set_ext(mm+3).xep(1) = ep_set(n).xep(1) + 2*pi;
    ep_set_ext(mm+4).xep(1) = ep_set(n).xep(1) + 4*pi;
    ep_set_ext(mm+5).xep(1) = ep_set(n).xep(1) + 6*pi;
    ep_set_ext(mm+6).xep(1) = ep_set(n).xep(1) + 10*pi;
    ep_set_ext(mm+7).xep(1) = ep_set(n).xep(1) + 18*pi;
    ep_set_ext(mm+8).xep(1) = ep_set(n).xep(1) + 26*pi;
    ep_set_ext(mm+9).xep(1) = ep_set(n).xep(1) + 34*pi;
    ep_set_ext(mm+10).xep(1) = ep_set(n).xep(1) - 4*pi;
end

f1 = figure(1);
hold on;
grid on;
ymin=-400;
ymax=200;
color_code = {'black','magenta','red','black'};
axis([-2*pi,1*pi,ymin,ymax]);
xticks(-2*pi:pi/2:2*pi);
xticklabels({'$-2\pi$', '', '$-\pi$', '','$0$', '','$\pi$', '','$2\pi$'});
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 14);
rangex=[-acos(Id*Lg*ki/Ug/kp),acos(Id*Lg*ki/Ug/kp),acos(Id*Lg*ki/Ug/kp),-acos(Id*Lg*ki/Ug/kp)];  rangey=[ymin,ymin,ymax,ymax];
%fill(rangex,rangey,[.9 .9 .9], 'linestyle', 'none', 'FaceAlpha',0.6);
for mm = 1 : length(ep_set_ext)
    xep = ep_set_ext(mm).xep;
    flag= ep_set_ext(mm).flag;
    scatter(xep(1),xep(2),color_code{flag+1},'LineWidth', 1.5);
    if flag == 1
        v = ep_set_ext(mm).v;
        vv = ep_set_ext(mm).vv;
        perturb = 1e-3;
        switch system
            case "GFL"
            [~ , x_p] = ode78(@f_backward,[0,1],xep+v*perturb,odeset('RelTol',1e-5));
            [~ , x_n] = ode78(@f_backward,[0,1],xep-v*perturb,odeset('RelTol',1e-5));
            %[~ , x_pp] = ode45(@f_forward,[0,0.1],xep+vv*perturb,odeset('RelTol',1e-5));
            %[~ , x_nn] = ode45(@f_forward,[0,0.1],xep-vv*perturb,odeset('RelTol',1e-5));
        end
        x_all = [flip(x_n,1);x_p];
        %x_allall = [flip(x_nn,1);x_pp];
        plot(x_all(:,1),x_all(:,2),'k-','linewidth',1.5);%scatter(x_all(:,1),x_all(:,2),'.');
        %plot(x_allall(:,1),x_allall(:,2),'r-','linewidth',1.5);%scatter(x_all(:,1),x_all(:,2),'.');
    end
end


%nullcline

% xdelta=-2*pi:0.01:2*pi;
% yomega= (-Ug*sin(xdelta)+Id*Xg+Iq*Rg)./(kp/ki*Ug*cos(xdelta)-Id*Lg);
% plot(xdelta,yomega,'r--','linewidth',1.5);%scatter(x_all(:,1),x_all(:,2),'.');


% %trajectory drawing
% [tt , x_all] = ode78(@f_forward,[0,70],[2*pi,-200],odeset('RelTol',1e-10));
% plot(x_all(:,1),x_all(:,2),'g-','linewidth',1.5);
% %trajectory drawing
% [tt , x_all] = ode78(@f_forward,[0,5],[20*pi,-300],odeset('RelTol',1e-6));
% plot(x_all(:,1),x_all(:,2),'g-','linewidth',1.5);
% 
% [tt , x_all] = ode78(@f_forward,[0,10],[2*pi,-100],odeset('RelTol',1e-5));
% plot(x_all(:,1),x_all(:,2),'b-','linewidth',1.5);

% [tt , x_all] = ode78(@f_backward,[0,1],[-10*pi,300],odeset('RelTol',1e-5));
% plot(x_all(:,1),x_all(:,2),'g-','linewidth',1.5);
switch fault_type
    case "voltage_sag"
    model = "original";

    t_start = 0.1;
    t_end = 0.4;
    delta_pre = [prefault_SEP(1); prefault_SEP(1)];
    omega_pre = [prefault_SEP(2); prefault_SEP(2)];
    t_prefault = [0;0.1];

    [t_fault , x_all] = ode78(@f_fault,[t_start,t_start+t_c],[prefault_SEP(1);0],odeset('RelTol',1e-6));
    Vq_fault = (Xg*Id+Rg*Iq-Ug_fault*sin(x_all(:,1))+Id*Lg*x_all(:,2))./(1-Id*Lg*kp);
    omega_fault_ori=kp.*Vq_fault+x_all(:,2);
    Vq_fault(find(omega_fault_ori>w_max))=(Xg*Id+Rg*Iq-Ug_fault*sin(x_all(find(omega_fault_ori>w_max),1))+Id*Lg*w_max);
    Vq_fault(find(omega_fault_ori<w_min))=(Xg*Id+Rg*Iq-Ug_fault*sin(x_all(find(omega_fault_ori<w_min),1))+Id*Lg*w_min);
    omega_fault=kp.*Vq_fault+x_all(:,2);
    omega_fault(find(omega_fault_ori>w_max))=w_max;
    omega_fault(find(omega_fault_ori<w_min))=w_min;
    delta_fault =  x_all(:,1);

    [t_postfault , x_all2] = ode78(@f_post,[t_fault(end),t_end],x_all(end,:),odeset('RelTol',1e-10));
    Vq_post = (Xg*Id+Rg*Iq-Ug*sin(x_all2(:,1))+Id*Lg*x_all2(:,2))./(1-Id*Lg*kp);
    omega_post_ori=kp.*Vq_post+x_all2(:,2);
    Vq_post(find(omega_post_ori>w_max))=(Xg*Id+Rg*Iq-Ug*sin(x_all2(find(omega_post_ori>w_max),1))+Id*Lg*w_max);
    Vq_post(find(omega_post_ori<w_min))=(Xg*Id+Rg*Iq-Ug*sin(x_all2(find(omega_post_ori<w_min),1))+Id*Lg*w_min);
    omega_post=kp.*Vq_post+x_all2(:,2);
    omega_post(find(omega_post>w_max))=w_max;
    omega_post(find(omega_post<w_min))=w_min;
    delta_post = x_all2(:,1);

    plot(delta_fault,omega_fault,'k-','linewidth',1.5);
    %plot(x_all(:,1),x_all(:,2),'m-','linewidth',1.5);
%     plot(x_all2(1,1),omega_post(1),'k.','MarkerSize',15);
%     plot(x_all(1,1),omega_fault(1),'k.','MarkerSize',15);
%     plot(x_all(end,1),omega_fault(end),'k.','MarkerSize',15);
%     plot(x_all2(:,1),omega_post,'y-','linewidth',1.5);
%     plot(x_all2(:,1),x_all2(:,2),'m-','linewidth',1.5);

    plot(delta_post(1),omega_post(1),'k.','MarkerSize',15);
    plot(delta_fault(1),omega_fault(1),'k.','MarkerSize',15);
    plot(delta_post,omega_post,'k-','linewidth',1.5)

%     [tt , xx] = ode78(@f_post,[0,1],[0.3877,-66.86],odeset('RelTol',1e-10));
%     plot(xx(:,1),xx(:,2),'g-','linewidth',1.5);
    model = "no_grid_dy";

    [t_fault2 , x_all] = ode78(@f_fault,[t_start,t_start+t_c],[prefault_SEP(1);0],odeset('RelTol',1e-6));
    Vq_fault2 = (Xg*Id+Rg*Iq-Ug_fault*sin(x_all(:,1))+Id*Lg*x_all(:,2))./(1-Id*Lg*kp);
    omega_fault_ori2=kp.*Vq_fault2+x_all(:,2);
    Vq_fault2(find(omega_fault_ori2>w_max))=(Xg*Id+Rg*Iq-Ug_fault*sin(x_all(find(omega_fault_ori2>w_max),1))+Id*Lg*w_max);
    Vq_fault2(find(omega_fault_ori2<w_min))=(Xg*Id+Rg*Iq-Ug_fault*sin(x_all(find(omega_fault_ori2<w_min),1))+Id*Lg*w_min);
    omega_fault2=kp.*Vq_fault2+x_all(:,2);
    omega_fault2(find(omega_fault_ori2>w_max))=w_max;
    omega_fault2(find(omega_fault_ori2<w_min))=w_min;
    delta_fault2 =  x_all(:,1);

    [t_postfault2 , x_all2] = ode78(@f_post,[t_fault(end),t_end],x_all(end,:),odeset('RelTol',1e-10));
    Vq_post2 = (Xg*Id+Rg*Iq-Ug*sin(x_all2(:,1))+Id*Lg*x_all2(:,2))./(1-Id*Lg*kp);
    omega_post_ori2=kp.*Vq_post2+x_all2(:,2);
    Vq_post2(find(omega_post_ori2>w_max))=(Xg*Id+Rg*Iq-Ug*sin(x_all2(find(omega_post_ori2>w_max),1))+Id*Lg*w_max);
    Vq_post2(find(omega_post_ori2<w_min))=(Xg*Id+Rg*Iq-Ug*sin(x_all2(find(omega_post_ori2<w_min),1))+Id*Lg*w_min);
    omega_post2=kp.*Vq_post2+x_all2(:,2);
    omega_post2(find(omega_post2>w_max))=w_max;
    omega_post2(find(omega_post2<w_min))=w_min;
    delta_post2 = x_all2(:,1);

    %plot(delta_fault2,omega_fault2,'y-','linewidth',1.5);
    %plot(delta_post2,omega_post2,'y-','linewidth',1.5);


    model = "original";



    case "line_cut"
    model = "original";

    
    t_start = 0.1;
    t_end = 0.4;
    delta_pre = [prefault_SEP(1); prefault_SEP(1)];
    omega_pre = [prefault_SEP(2); prefault_SEP(2)];
    t_prefault = [0;0.1];

    [t_fault , x_all] = ode78(@f_fault,[t_start,t_start+t_c],[prefault_SEP(1)-Ug_fault_angle;0],odeset('RelTol',1e-6));
    Vq_fault = (Xg_f*Id+Rg_f*Iq-Ug_fault*sin(x_all(:,1))+Id*Lg_f*x_all(:,2))./(1-Id*Lg_f*kp);
    omega_fault_ori=kp.*Vq_fault+x_all(:,2);
    Vq_fault(find(omega_fault_ori>w_max))=(Xg_f*Id+Rg_f*Iq-Ug_fault*sin(x_all(find(omega_fault_ori>w_max),1))+Id*Lg_f*w_max);
    Vq_fault(find(omega_fault_ori<w_min))=(Xg_f*Id+Rg_f*Iq-Ug_fault*sin(x_all(find(omega_fault_ori<w_min),1))+Id*Lg_f*w_min);
    omega_fault=kp.*Vq_fault+x_all(:,2);
    omega_fault(find(omega_fault_ori>w_max))=w_max;
    omega_fault(find(omega_fault_ori<w_min))=w_min;
    delta_fault =  x_all(:,1)+Ug_fault_angle;

    [t_postfault , x_all2] = ode78(@f_post,[t_fault(end),t_end],[delta_fault(end),x_all(end,2)],odeset('RelTol',1e-10));
    Vq_post = (Xg*Id+Rg*Iq-Ug*sin(x_all2(:,1))+Id*Lg*x_all2(:,2))./(1-Id*Lg*kp);
    omega_post_ori=kp.*Vq_post+x_all2(:,2);
    Vq_post(find(omega_post_ori>w_max))=(Xg*Id+Rg*Iq-Ug*sin(x_all2(find(omega_post_ori>w_max),1))+Id*Lg*w_max);
    Vq_post(find(omega_post_ori<w_min))=(Xg*Id+Rg*Iq-Ug*sin(x_all2(find(omega_post_ori<w_min),1))+Id*Lg*w_min);
    omega_post=kp.*Vq_post+x_all2(:,2);
    omega_post(find(omega_post>w_max))=w_max;
    omega_post(find(omega_post<w_min))=w_min;
    delta_post = x_all2(:,1);

    plot(prefault_SEP(1),0,'k.','MarkerSize',15);
    plot(delta_fault,omega_fault_ori,'k-','linewidth',1.5);
    plot(delta_post(1),omega_post(1),'k.','MarkerSize',15);
    plot(delta_fault(1),omega_fault(1),'k.','MarkerSize',15);
    plot(delta_post,omega_post_ori,'k-','linewidth',1.5);

%     plot(delta_fault,x_all(:,2),'g-','linewidth',1.5);
%     plot(delta_post,x_all2(:,2),'y-','linewidth',1.5);
%     [tt , xx] = ode78(@f_post,[0,1],[1.9649,38.1636],odeset('RelTol',1e-10));

    model = "no_grid_dy";

    [t_fault2 , x_all] = ode78(@f_fault,[t_start,t_start+t_c],[prefault_SEP(1)-Ug_fault_angle;0],odeset('RelTol',1e-6));
    Vq_fault2 = (Xg_f*Id+Rg_f*Iq-Ug_fault*sin(x_all(:,1)));
    omega_fault_ori2=kp.*Vq_fault2+x_all(:,2);
    omega_fault2=kp.*Vq_fault+x_all(:,2);
    omega_fault2(find(omega_fault_ori2>w_max))=w_max;
    omega_fault2(find(omega_fault_ori2<w_min))=w_min;
    delta_fault2 =  x_all(:,1)+Ug_fault_angle;

    [t_postfault2 , x_all2] = ode78(@f_post,[t_fault(end),t_end],[delta_fault2(end),x_all(end,2)],odeset('RelTol',1e-10));
    Vq_post2 = (Xg*Id+Rg*Iq-Ug*sin(x_all2(:,1)));
    omega_post_ori2=kp.*Vq_post2+x_all2(:,2);
    omega_post2=kp.*Vq_post2+x_all2(:,2);
    omega_post2(find(omega_post2>w_max))=w_max;
    omega_post2(find(omega_post2<w_min))=w_min;
    delta_post2 = x_all2(:,1);

    plot(delta_fault2,omega_fault_ori2,'y-','linewidth',1.5);
    plot(delta_post2,omega_post_ori2,'y-','linewidth',1.5);


    model = "original";



    case "frequency"
    [t1 , x_all] = ode78(@f_post,[0,0.4],prefault_SEP(1)-[0;0],odeset('RelTol',1e-5));
    Vq_post = (Xg*Id+Rg*Iq-Ug*sin(x_all(:,1))+Id*Lg*x_all(:,2))./(1-Id*Lg*kp);
    omegapll_post=kp.*Vq_post+x_all(:,2);
    omegapll_post(find(omegapll_post>w_max))=w_max;
    omegapll_post(find(omegapll_post<w_min))=w_min;
    omega_post = omegapll_post-W_g;
    plot(x_all(1,1),omega_post(1),'k.','MarkerSize',15);
    plot(x_all(:,1),omega_post,'r-','linewidth',1.5);
end

%% Vq=0
% syms deltax omegax;
% VVq = (Xg*Id+Rg*Iq-Ug*sin(deltax)+Id*Lg*omegax);
% VVq=vpa(VVq);
% VVq=matlabFunction(VVq);
% x1=-2*pi:0.01*pi:4*pi;
% x2=ymin:1:ymax;
% [y1,y2]=meshgrid(x1,x2);
% zz = zeros(length(x2),length(x1));
% for a = 1: length(x1)
%     for b = 1: length(x2)
%         V = VVq(y1(b,a), y2(b,a));
%         zz(b,a) = V;
%     end
% end
% contour(y1,y2,zz,[0 0],'g','linewidth',1.5,"ShowText",true);




%%
% model = "no_grid_dy";
% for mm = 1 : length(ep_set_ext)
%     xep = ep_set_ext(mm).xep;
%     flag= ep_set_ext(mm).flag;
%     if flag == 1
%         v = ep_set_ext(mm).v;
%         perturb = 1e-3;
%         switch system
%             case "GFL"
%             [~ , x_p] = ode45(@f_backward,[0,0.5],xep+v*perturb,odeset('RelTol',1e-5));
%             [~ , x_n] = ode45(@f_backward,[0,0.5],xep-v*perturb,odeset('RelTol',1e-5));
%         end
%         x_all = [flip(x_n,1);x_p];
%         plot(x_all(:,1),x_all(:,2),'b-','linewidth',1);%scatter(x_all(:,1),x_all(:,2),'.');
%     end
% end
%trajectory drawing
% [tt , x_all] = ode78(@f_forward,[0,0.5],[0; -500],odeset('RelTol',1e-6));
% plot(x_all(:,1),x_all(:,2),'k:','linewidth',1.5);

% x1=-2*pi:0.05*pi:2*pi;
% x2=-300:25:300;
% [y1,y2]=meshgrid(x1,x2);
% yy1 = zeros(length(x2),length(x1));
% yy2 = zeros(length(x2),length(x1));
% for a = 1: length(x1)
%     for b = 1: length(x2)
%     switch system
%         case "GFL"
%         yy = f_GFL([y1(b,a) y2(b,a)]);
%     end
%         yy1(b,a) = yy(1);
%         yy2(b,a) = yy(2);
%     end
% end
% quiver(y1,y2,yy1,yy2,'color',[0.8 0.8 0.8],'AutoScale','on','AutoScaleFactor',0.5);
% xlabel('\delta (rad)');
% ylabel('\omega (rad/s)');
% %trajectory drawing
% [t1 , x_all] = ode78(@f_fault,[0,t_c],ep_set(1).xep-[0;0],odeset('RelTol',1e-5));
% [t2 , x_all2] = ode78(@f_post,[t1(end),1],x_all(end,:),odeset('RelTol',1e-5));
% Vq = (Xg*Id+Rg*Iq-Ug*sin(x_all(:,1)));
% omega=kp.*Vq+x_all(:,2);
% 
% plot(x_all(:,1),omega,'b:','linewidth',1.5)
% Vq = (Xg*Id+Rg*Iq-Ug*sin(x_all2(:,1)));
% omega=kp.*Vq+x_all2(:,2);
% 
% plot(x_all2(1,1),omega(1),'b.','MarkerSize',15);
% plot(x_all2(:,1),omega,'b-.','linewidth',1.5);
model = "original";
%%
% energy function 1
syms deltax omegax;
M = (1-kp*Lg*Id)/ki;
V1 = -Xg*Id*deltax-Rg*Iq*deltax-Ug*cos(deltax) + 1/2*M*omegax^2;
V1=vpa(V1);
VV1=matlabFunction(V1);
V1d = jacobian(V1);
VV1d = matlabFunction(V1d);

x1=-2*pi:0.02*pi:2*pi;
x2=ymin:5:ymax;
[y1,y2]=meshgrid(x1,x2);
zz = zeros(length(x2),length(x1));
dzz = zeros(length(x2),length(x1));
for a = 1: length(x1)
    for b = 1: length(x2)
        V = VV1(y1(b,a), y2(b,a));
        dV = VV1d(y1(b,a), y2(b,a))*f_GFL([y1(b,a) y2(b,a)]);
        zz(b,a) = V;
        dzz(b,a)=dV;
    end
end
Vcr1 = VV1(ep_set(2).xep(1),ep_set(2).xep(2));
contour(y1,y2,zz,[Vcr1 Vcr1],'b-','linewidth',1,"ShowText",false);
Vcr2 = VV1(acos(Id*Lg*ki/Ug/kp),0);
contour(y1,y2,zz,[Vcr2 Vcr2],'b-','linewidth',1.5,"ShowText",false);

% energy function 2
syms deltax omegax;
M = (1-kp*Lg*Id)/ki;
delta_s = ep_set(1).xep(1);
V2 = 0.5*(M*omegax-kp/ki*(Xg*Id-Ug*sin(deltax))-Lg*Id*(deltax-delta_s))^2-M*(Xg*Id*deltax+Ug*cos(deltax));
V2=vpa(V2);
VV2=matlabFunction(V2);
V2d = jacobian(V2);
VV2d = matlabFunction(V2d);
fun = @(x)(Xg*Id-Ug*sin(x)+ki/kp*Lg*Id*(x-ep_set(1).xep(1)));
x0 = ep_set(2).xep(1);
[x_critical fval exitflag output] = fzero(fun,x0);
plot(x_critical,0,'m*','MarkerSize',10);
x1=-2*pi:0.02*pi:2*pi;
x2=ymin:5:ymax;
[y1,y2]=meshgrid(x1,x2);
zz = zeros(length(x2),length(x1));
dzz = zeros(length(x2),length(x1));
for a = 1: length(x1)
    for b = 1: length(x2)
        V = VV2(y1(b,a), y2(b,a));
        dV = VV2d(y1(b,a), y2(b,a))*f_GFL([y1(b,a) y2(b,a)]);
        zz(b,a) = V;
        dzz(b,a)=dV;
    end
end
V2cr = VV2(x_critical,0);
contour(y1,y2,zz,[V2cr V2cr],'m-','linewidth',1.5,"ShowText",false);
% % model = "no_grid_dy";
% % energy function 3
syms deltax omegax;
vvq = (Xg*Id-Ug*sin(deltax)+Id*Lg*omegax);
V3 = 1/ki/2*(omegax-kp*vvq)^2 - (Ug*cos(deltax)-Ug*cos(delta_s)+Xg*Id*(deltax-delta_s)+1/2*Id*Lg*omegax*(deltax-delta_s));%
V3=vpa(V3);
VV3=matlabFunction(V3);
V3d = jacobian(V3);
VV3d = matlabFunction(V3d);
x1=-2*pi:0.01*pi:2*pi;
x2=ymin:2:ymax;
[y1,y2]=meshgrid(x1,x2);
zz = zeros(length(x2),length(x1));
dzz = zeros(length(x2),length(x1));
for a = 1: length(x1)
    for b = 1: length(x2)
        V = VV3(y1(b,a), y2(b,a));
        dV = VV3d(y1(b,a), y2(b,a))*f_GFL([y1(b,a) y2(b,a)]);
        zz(b,a) = V;
        dzz(b,a)=dV;
    end
end
V3cr = VV3(ep_set(2).xep(1),0);
contour(y1,y2,zz,[V3cr V3cr],'r-','linewidth',1.5,"ShowText",false);
contour(y1,y2,dzz,[-10 0 10],'r:','linewidth',0.5,"ShowText",true);

% energy function 4
syms deltax omegax;
M = (1-kp*Lg*Id)/ki;
D = kp/ki*Ug*cos(deltax)-Id*Lg;
delta_s = ep_set(1).xep(1); D_s = kp/ki*Ug*cos(delta_s)-Id*Lg;
beta = D_s/(D_s^2/4+M*Ug*cos(delta_s))/1.2;
V4 = 0.5*M*omegax^2-Xg*Id*deltax-Ug*cos(deltax) - beta*M*omegax*(Xg*Id-Ug*sin(deltax));
V4=vpa(V4);
VV4=matlabFunction(V4);
V4d = jacobian(V4);
VV4d = matlabFunction(V4d);
x1=-2*pi:0.01*pi:2*pi;
x2=ymin:2:ymax;
[y1,y2]=meshgrid(x1,x2);
zz = zeros(length(x2),length(x1));
dzz = zeros(length(x2),length(x1));
for a = 1: length(x1)
    for b = 1: length(x2)
        V = VV4(y1(b,a), y2(b,a));
        dV = VV4d(y1(b,a), y2(b,a))*f_GFL([y1(b,a) y2(b,a)]);
        zz(b,a) = V;
        dzz(b,a)=dV;
    end
end
V4cr = VV4(1.64,-90);
contour(y1,y2,zz,[V4cr V4cr],'c-','linewidth',1.5,"ShowText",false);
contour(y1,y2,dzz,[0 0],'c:','linewidth',1.5,"ShowText",false);

beta = D_s/(D_s^2/4+M*Ug*cos(delta_s))/2;
V4 = 0.5*M*omegax^2-Xg*Id*deltax-Ug*cos(deltax) - beta*M*omegax*(Xg*Id-Ug*sin(deltax));
V4=vpa(V4);
VV4=matlabFunction(V4);
V4d = jacobian(V4);
VV4d = matlabFunction(V4d);
x1=-2*pi:0.01*pi:2*pi;
x2=ymin:2:ymax;
[y1,y2]=meshgrid(x1,x2);
zz = zeros(length(x2),length(x1));
dzz = zeros(length(x2),length(x1));
for a = 1: length(x1)
    for b = 1: length(x2)
        V = VV4(y1(b,a), y2(b,a));
        dV = VV4d(y1(b,a), y2(b,a))*f_GFL([y1(b,a) y2(b,a)]);
        zz(b,a) = V;
        dzz(b,a)=dV;
    end
end
V4cr = VV4(1.68,-60);
contour(y1,y2,zz,[V4cr V4cr],'g-','linewidth',1.5,"ShowText",false);
contour(y1,y2,dzz,[0 0],'g:','linewidth',1.5,"ShowText",false);

beta = D_s/(D_s^2/4+M*Ug*cos(delta_s))/3;
V4 = 0.5*M*omegax^2-Xg*Id*deltax-Ug*cos(deltax) - beta*M*omegax*(Xg*Id-Ug*sin(deltax));
V4=vpa(V4);
VV4=matlabFunction(V4);
V4d = jacobian(V4);
VV4d = matlabFunction(V4d);
x1=-2*pi:0.01*pi:2*pi;
x2=ymin:2:ymax;
[y1,y2]=meshgrid(x1,x2);
zz = zeros(length(x2),length(x1));
dzz = zeros(length(x2),length(x1));
for a = 1: length(x1)
    for b = 1: length(x2)
        V = VV4(y1(b,a), y2(b,a));
        dV = VV4d(y1(b,a), y2(b,a))*f_GFL([y1(b,a) y2(b,a)]);
        zz(b,a) = V;
        dzz(b,a)=dV;
    end
end
V4cr = VV4(1.66,-48);
contour(y1,y2,zz,[V4cr V4cr],'y-','linewidth',1.5,"ShowText",false);
contour(y1,y2,dzz,[0 0],'y:','linewidth',1.5,"ShowText",false);
 end

%% function
function yes = isnewxep(ep_set,xep,torr)
    if isempty(ep_set)
        yes = 1;
        return;
    end
    minerr = inf;
    for m = 1 : length(ep_set)
        err = abs(xep - ep_set(m).xep);
        err = min(err, abs(2*pi-err));
        err = max(err);
        if minerr > err
            minerr = err;
        end
    end
    if(minerr>torr)
        yes = 1;
    else
        yes = 0;
    end
end



function dfdt = f(x)
global system;
switch system 
    case "GFL"
        dfdt = f_GFL(x);
end
end
function dfdt = f_normal(x)
global system;
switch system 
    case "GFL"
        dfdt = f_GFL_normal(x);
end
end

function out = maxabs(in)

    out = abs(in);
    
    while length(out) > 1
        out = max(out);
    end

end

function disp_v(msg,v)
    disp([msg '=']);
    disp(v);
end

function dfdt = f_backward(t,x)
    dfdt = -f(x);
end

function dfdt = f_forward(t,x)
    dfdt = f(x);
end

function dfdt = f_forward_limit(t,x)
    dfdt = f_GFL_limit(x);
end

function dfdt = f_backward_limit(t,x)
    dfdt = -f_GFL_limit(x);
end

function dfdt = f_fault(t,x)
global system;
    switch system 
        case "GFL"
          dfdt = f_GFL_fault(x);
    end
end
function dfdt = f_post(t,x)
global system;
    switch system 
        case "GFL"
          dfdt = f_GFL_normal(x);
    end
end
function dfdt = f_prefault(t,x)
global system;
    switch system 
        case "GFL"
          dfdt = f_GFL_prefault(x);
    end
end