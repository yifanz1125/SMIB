%% parameter


t_end = 1;

%grid

Xg = 0.5;
Rg = 0.02;%0.08
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

Lf = 0.05/Ws;
Xf = 0.05;
Xsum = 0.05+Xg;
Lsum = Xsum/Ws;
beta = 800*2*pi;
kpcc = beta*Lf;
kiccc = beta^2*Lf/4;


%GFM
m_gfm = 0.05; 
w_droop = 0.8*2*pi;
1/(w_droop*m_gfm)
D = 20; %1/m_gfm;
J = 4;  %
Vgfm = 1;
Pm = 1;


%GFM with Q-V droop
m_gfm = 0.05;
Qref = -0.2;
Vgfm = 1;
Pm = 1;
k_q = 0.5;
tau_q = 1/(50*2*pi);


%VFM
Vdc_ref = 2.5;
Y_dc = 12.5;%12.5;  %12.5
C_dc = Y_dc/Ws;
Kpp = 2*2*pi;
Kip = 15;%15
Vvfm = 1;
Pin = 1;
Phi = -pi/4;
Ilim = 2;

y_lim = 3;

%system
global system;
global fault_type; %line_cut voltage_sag frequency
fault_type = "voltage_sag"; %"voltage_sag";%"line_cut";%"line_cut";
system = "GFM";  %GFMQ
model = "original";% "original"

if system == "VFM2"
    Kip =  Kip*2*Vdc_ref;
end

switch fault_type
    case "voltage_sag"
        %fault sag
        Ug_fault = 0.1;%0.1;
        X1 = 0.1;
        R1 = 0.01;
        Xgg = (Xg - X1)*2;
        Rgg = (Rg - R1)*2;
        t_c = 0.05;%0.072;%0.086;0.0795
    case "line_cut"
    %fault line cut 
        t_c = 0.06;
        X1 = 0.1;
        R1 = 0.01;
        Xgg = Xg-X1;
        Rgg = Rg - R1;
        Xg0 = Xgg/2+X1;
        Rg0 = Rgg/2+R1;
        position = 1; %fault to inf bus
        Rf = 1e-5;%1e-5/(690^2/1e6);
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
    system = "VFM";
end
%% SEP
if system == "GFMQ"
    x=(0:0.1:1)*2*pi;
    y=(0:0.1:1)*2;
    n = length(x);
    x_set = zeros(2,n^2);
    for m = 0:(n^2 - 1)
        index = floor(m);
        index = mod(index,n)+1;
        x_set(1,m+1) = x(index);
        index = floor(m/n);
        index = mod(index,n)+1;
        x_set(2,m+1) = y(index);
    end
else
x=(0:0.1:1)*2*pi;
n = length(x);
x_set = zeros(2,n);
x_set(1,:) = x;
end


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
if system == "GFMQ"
    ymin = 0;
    ymax = 1.5;
elseif system == "GFM"
ymin=-0.2;
ymax=0.15;
else
ymin=-300;
ymax=300;
end
color_code = {'black','magenta','red','black'};
set(gcf,'position',[1 1 800 500]);
axis([-1*pi,3*pi/2,ymin,ymax]);
xticks(-2*pi:pi/2:2*pi);
yticks(-0.2:0.1:0.2);
xticklabels({'$-2\pi$', '', '$-\pi$', '','$0$', '','$\pi$', '','$2\pi$'});
yticklabels({'-0.2', '-0.1', '0', '0.1','$0.2$'});
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 20);
if system == "GFL"   
rangex=[-acos(Id*Lg*ki/Ug/kp),acos(Id*Lg*ki/Ug/kp),acos(Id*Lg*ki/Ug/kp),-acos(Id*Lg*ki/Ug/kp)];  rangey=[ymin,ymin,ymax,ymax];
fill(rangex,rangey,[.9 .9 .9], 'linestyle', 'none', 'FaceAlpha',0.6);
end
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
            [~ , x_p] = ode78(@f_backward,[0,2],xep+v*perturb,odeset('RelTol',1e-5));
            [~ , x_n] = ode78(@f_backward,[0,2],xep-v*perturb,odeset('RelTol',1e-5));
            case "GFM"
            [~ , x_p] = ode78(@f_backward,[0,0.8],xep+v*perturb,odeset('RelTol',1e-5));
            [~ , x_n] = ode78(@f_backward,[0,0.8],xep-v*perturb,odeset('RelTol',1e-5));
            %[~ , x_pp] = ode45(@f_forward,[0,0.1],xep+vv*perturb,odeset('RelTol',1e-5));
            %[~ , x_nn] = ode45(@f_forward,[0,0.1],xep-vv*perturb,odeset('RelTol',1e-5));
            case "GFMQ"
            [~ , x_p] = ode78(@f_backward,[0,0.4],xep+v*perturb,odeset('RelTol',1e-5));
            [~ , x_n] = ode78(@f_backward,[0,0.4],xep-v*perturb,odeset('RelTol',1e-5));
            case "VFM"
            figure(f1)
            ymin=-4;
            ymax=6;
            axis([-1*pi,3/2*pi,ymin,ymax]);
            
            yticks(ymin:2:ymax);
            [~ , x_p] = ode78(@f_backward,[0,2],xep+v*perturb,odeset('RelTol',1e-5));
            [~ , x_n] = ode78(@f_backward,[0,2],xep-v*perturb,odeset('RelTol',1e-5));  
            case "VFM2"
            figure(f1)
            ymin=-1;
            ymax=1;
            axis([-1*pi,3/2*pi,ymin,ymax]);
            yticks(ymin:0.5:ymax);
            [~ , x_p] = ode78(@f_backward,[0,2],xep+v*perturb,odeset('RelTol',1e-5));
            [~ , x_n] = ode78(@f_backward,[0,2],xep-v*perturb,odeset('RelTol',1e-5));  
        end
        x_all = [flip(x_n,1);x_p];
        plot(x_all(:,1),x_all(:,2),'k-','linewidth',1.5);%scatter(x_all(:,1),x_all(:,2),'.');
    end
end

%%

if system == "GFL"
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

% % energy function 2
% syms deltax omegax;
% M = (1-kp*Lg*Id)/ki;
% delta_s = ep_set(1).xep(1);
% V2 = 0.5*(M*omegax-kp/ki*(Xg*Id-Ug*sin(deltax))-Lg*Id*(deltax-delta_s))^2-M*(Xg*Id*deltax+Ug*cos(deltax));
% V2=vpa(V2);
% VV2=matlabFunction(V2);
% V2d = jacobian(V2);
% VV2d = matlabFunction(V2d);
% fun = @(x)(Xg*Id-Ug*sin(x)+ki/kp*Lg*Id*(x-ep_set(1).xep(1)));
% x0 = ep_set(2).xep(1);
% [x_critical fval exitflag output] = fzero(fun,x0);
% plot(x_critical,0,'m*','MarkerSize',10);
% x1=-2*pi:0.02*pi:2*pi;
% x2=ymin:5:ymax;
% [y1,y2]=meshgrid(x1,x2);
% zz = zeros(length(x2),length(x1));
% dzz = zeros(length(x2),length(x1));
% for a = 1: length(x1)
%     for b = 1: length(x2)
%         V = VV2(y1(b,a), y2(b,a));
%         dV = VV2d(y1(b,a), y2(b,a))*f_GFL([y1(b,a) y2(b,a)]);
%         zz(b,a) = V;
%         dzz(b,a)=dV;
%     end
% end
% V2cr = VV2(x_critical,0);
% contour(y1,y2,zz,[V2cr V2cr],'m-','linewidth',1.5,"ShowText",false);
% % model = "no_grid_dy";
% energy function 3
syms deltax omegax;
delta_s = ep_set(1).xep(1);
vvq = (Xg*Id-Ug*sin(deltax)+Id*Lg*omegax);
deltauep = ep_set(2).xep(1);
V3 = 1/ki/2*(omegax-kp*vvq)^2 - (Ug*cos(deltax)-Ug*cos(delta_s)+Xg*Id*(deltax-delta_s)-1/2*Id*Lg*omegax*(deltauep-deltax));%
%V3 = 1/ki/2*(omegax-kp*vvq)^2 - (Ug*cos(deltax)-Ug*cos(delta_s)+Xg*Id*(deltax-delta_s))+1/2*Id*Lg*abs(omegax)*(deltauep-deltax)-vvq^2*kp^2/2/ki;%
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
%contour(y1,y2,dzz,[-10 0 10],'r:','linewidth',0.5,"ShowText",true);

% % energy function 4
% syms deltax omegax;
% M = (1-kp*Lg*Id)/ki;
% D = kp/ki*Ug*cos(deltax)-Id*Lg;
% delta_s = ep_set(1).xep(1); D_s = kp/ki*Ug*cos(delta_s)-Id*Lg;
% beta = D_s/(D_s^2/4+M*Ug*cos(delta_s))/1.2;
% V4 = 0.5*M*omegax^2-Xg*Id*deltax-Ug*cos(deltax) - beta*M*omegax*(Xg*Id-Ug*sin(deltax));
% V4=vpa(V4);
% VV4=matlabFunction(V4);
% V4d = jacobian(V4);
% VV4d = matlabFunction(V4d);
% x1=-2*pi:0.01*pi:2*pi;
% x2=ymin:2:ymax;
% [y1,y2]=meshgrid(x1,x2);
% zz = zeros(length(x2),length(x1));
% dzz = zeros(length(x2),length(x1));
% for a = 1: length(x1)
%     for b = 1: length(x2)
%         V = VV4(y1(b,a), y2(b,a));
%         dV = VV4d(y1(b,a), y2(b,a))*f_GFL([y1(b,a) y2(b,a)]);
%         zz(b,a) = V;
%         dzz(b,a)=dV;
%     end
% end
% V4cr = VV4(1.64,-90);
% contour(y1,y2,zz,[V4cr V4cr],'c-','linewidth',1.5,"ShowText",false);
% contour(y1,y2,dzz,[0 0],'c:','linewidth',1.5,"ShowText",false);
% 
% beta = D_s/(D_s^2/4+M*Ug*cos(delta_s))/2;
% V4 = 0.5*M*omegax^2-Xg*Id*deltax-Ug*cos(deltax) - beta*M*omegax*(Xg*Id-Ug*sin(deltax));
% V4=vpa(V4);
% VV4=matlabFunction(V4);
% V4d = jacobian(V4);
% VV4d = matlabFunction(V4d);
% x1=-2*pi:0.01*pi:2*pi;
% x2=ymin:2:ymax;
% [y1,y2]=meshgrid(x1,x2);
% zz = zeros(length(x2),length(x1));
% dzz = zeros(length(x2),length(x1));
% for a = 1: length(x1)
%     for b = 1: length(x2)
%         V = VV4(y1(b,a), y2(b,a));
%         dV = VV4d(y1(b,a), y2(b,a))*f_GFL([y1(b,a) y2(b,a)]);
%         zz(b,a) = V;
%         dzz(b,a)=dV;
%     end
% end
% V4cr = VV4(1.68,-60);
% contour(y1,y2,zz,[V4cr V4cr],'g-','linewidth',1.5,"ShowText",false);
% contour(y1,y2,dzz,[0 0],'g:','linewidth',1.5,"ShowText",false);
% 
% beta = D_s/(D_s^2/4+M*Ug*cos(delta_s))/3;
% V4 = 0.5*M*omegax^2-Xg*Id*deltax-Ug*cos(deltax) - beta*M*omegax*(Xg*Id-Ug*sin(deltax));
% V4=vpa(V4);
% VV4=matlabFunction(V4);
% V4d = jacobian(V4);
% VV4d = matlabFunction(V4d);
% x1=-2*pi:0.01*pi:2*pi;
% x2=ymin:2:ymax;
% [y1,y2]=meshgrid(x1,x2);
% zz = zeros(length(x2),length(x1));
% dzz = zeros(length(x2),length(x1));
% for a = 1: length(x1)
%     for b = 1: length(x2)
%         V = VV4(y1(b,a), y2(b,a));
%         dV = VV4d(y1(b,a), y2(b,a))*f_GFL([y1(b,a) y2(b,a)]);
%         zz(b,a) = V;
%         dzz(b,a)=dV;
%     end
% end
% V4cr = VV4(1.66,-48);
% contour(y1,y2,zz,[V4cr V4cr],'y-','linewidth',1.5,"ShowText",false);
% contour(y1,y2,dzz,[0 0],'y:','linewidth',1.5,"ShowText",false);

elseif system == "VFM"

    % energy function 1
    % syms deltax yx;
    % delta_s = prefault_SEP(1);
    % M = C_dc/2*Kip;
    % ppp = Rg*(Vvfm^2-Vvfm*Ug*cos(deltax))/(Rg^2+Xg^2)+Xg*Vvfm*Ug*sin(deltax)/(Rg^2+Xg^2);
    % V1 = - Pin*(deltax-delta_s) + Rg*(Vvfm^2*(deltax-delta_s)-Vvfm*Ug*(sin(deltax)-sin(delta_s)))/(Rg^2+Xg^2)-Xg*Vvfm*Ug*(cos(deltax)-cos(delta_s))/(Rg^2+Xg^2) + 1/2*M*(Kip*yx+Kpp*(Pin - ppp))^2;
    % V1=vpa(V1);
    % VV1=matlabFunction(V1);
    % V1d = jacobian(V1);
    % VV1d = matlabFunction(V1d);
    % 
    % x1=-2*pi:0.01*pi:2*pi;
    % x2=-4:0.02:6;%-8:0.1:8;
    % [y1,y2]=meshgrid(x1,x2);
    % zz = zeros(length(x2),length(x1));
    % dzz = zeros(length(x2),length(x1));
    % for a = 1: length(x1)
    %     for b = 1: length(x2)
    %         V = VV1(y1(b,a), y2(b,a));
    %         dV = VV1d(y1(b,a), y2(b,a))*f_GFL([y1(b,a) y2(b,a)]);
    %         zz(b,a) = V;
    %         dzz(b,a)=dV;
    %     end
    % end
    % Vcr1 = VV1(ep_set(2).xep(1),0);
    % contour(y1,y2,zz,[Vcr1 Vcr1],'b-','linewidth',1,"ShowText",false);
%     Vcr2 = VV1(acos(Id*Lg*ki/Ug/kp),0);
%     contour(y1,y2,zz,[Vcr2 Vcr2],'b-','linewidth',1.5,"ShowText",false);


    % energy function 3
    syms deltax yx;
    delta_s = prefault_SEP(1);
    deltauep = ep_set(2).xep(1);
    %V3 = 1/ki/2*(omegax-kp*vvq)^2 - (Ug*cos(deltax)-Ug*cos(delta_s)+Xg*Id*(deltax-delta_s)+1/2*Id*Lg*omegax*(deltax-delta_s));%
    V3 = C_dc/4*Kip*yx^2 - Pin*(deltax-delta_s) + Rg*(Vvfm^2*(deltax-delta_s)-Vvfm*Ug*(sin(deltax)-sin(delta_s)))/(Rg^2+Xg^2)-Xg*Vvfm*Ug*(cos(deltax)-cos(delta_s))/(Rg^2+Xg^2) ;%- C_dc/2*Kpp*yx*(Rg*(Vvfm^2-Vvfm*Ug*cos(deltax))/(Rg^2+Xg^2)+Xg*Vvfm*Ug*sin(deltax)/(Rg^2+Xg^2) - Pin)/2;%
    V3_2 = V3-C_dc/2*Kpp*yx*(Rg*(Vvfm^2-Vvfm*Ug*cos(deltax))/(Rg^2+Xg^2)+Xg*Vvfm*Ug*sin(deltax)/(Rg^2+Xg^2) - Pin)/2;
    V3=vpa(V3);
    VV3=matlabFunction(V3);
    VV3_2=matlabFunction(V3_2);
    V3d = jacobian(V3);
    VV3d = matlabFunction(V3d);
    x1=-2*pi:0.01*pi:2*pi;
    x2=-4:0.02:6;%-8:0.1:8;
    [y1,y2]=meshgrid(x1,x2);
    zz = zeros(length(x2),length(x1));
    zz_2 = zeros(length(x2),length(x1));
    dzz = zeros(length(x2),length(x1));
    for a = 1: length(x1)
        for b = 1: length(x2)
            V = VV3(y1(b,a), y2(b,a));
            dV = VV3d(y1(b,a), y2(b,a))*f_VFM_normal([y1(b,a) y2(b,a)]);
            zz(b,a) = V;
            zz_2(b,a) = VV3_2(y1(b,a), y2(b,a));
            dzz(b,a)=dV;
        end
    end
    V3cr = VV3(ep_set(2).xep(1),0);
    contour(y1,y2,zz,[V3cr V3cr],'b-','linewidth',1.5,"ShowText",false);
    % contour(y1,y2,zz_2,[V3cr V3cr],'m-','linewidth',1.5,"ShowText",false);
    % contour(y1,y2,dzz,[-10 -5 0 5 10],'r:','linewidth',0.5,"ShowText",true);



    % % energy function 4
    % syms deltax yx;
    % delta_s = prefault_SEP(1);
    % deltauep = ep_set(2).xep(1);
    % %V3 = 1/ki/2*(omegax-kp*vvq)^2 - (Ug*cos(deltax)-Ug*cos(delta_s)+Xg*Id*(deltax-delta_s)+1/2*Id*Lg*omegax*(deltax-delta_s));%
    % V3 = C_dc/4*Kip*yx^2 - Pin*(deltax-delta_s) + Rg*(Vvfm^2*(deltax-delta_s)-Vvfm*Ug*(sin(deltax)-sin(delta_s)))/(Rg^2+Xg^2)-Xg*Vvfm*Ug*(cos(deltax)-cos(delta_s))/(Rg^2+Xg^2) ;%- C_dc/2*Kpp*yx*(Rg*(Vvfm^2-Vvfm*Ug*cos(deltax))/(Rg^2+Xg^2)+Xg*Vvfm*Ug*sin(deltax)/(Rg^2+Xg^2) - Pin)/2;%
    % V3_2 = C_dc/4*Kip*yx^2 - Pin*(deltax-delta_s) + Ilim*Ug*sin(deltax+Phi) - Ilim*Ug*sin(delta_s+Phi) + Ilim^2*Rg*(deltax-delta_s);
    % V3=vpa(V3);
    % VV3=matlabFunction(V3);
    % VV3_2=matlabFunction(V3_2);
    % V3d = jacobian(V3);
    % VV3d = matlabFunction(V3d);
    % x1=-2*pi:0.01*pi:2*pi;
    % x2=-4:0.02:6;%-8:0.1:8;
    % [y1,y2]=meshgrid(x1,x2);
    % zz = zeros(length(x2),length(x1));
    % zz_2 = zeros(length(x2),length(x1));
    % dzz = zeros(length(x2),length(x1));
    % for a = 1: length(x1)
    %     for b = 1: length(x2)
    %         V = VV3(y1(b,a), y2(b,a));
    %         dV = VV3d(y1(b,a), y2(b,a))*f_VFM_normal([y1(b,a) y2(b,a)]);
    %         zz(b,a) = V;
    %         zz_2(b,a) = VV3_2(y1(b,a), y2(b,a));
    %         dzz(b,a)=dV;
    %     end
    % end
    % V3cr = VV3_2(ep_set(2).xep(1),0);
    % contour(y1,y2,zz_2,[V3cr V3cr],'r-','linewidth',1.5,"ShowText",false);
    % deltacc = acos((Vvfm^2+Ug^2-Ilim^2*(Xg^2+Rg^2))/(2*Vvfm*Ug));
    % V32cr = VV3_2(deltacc,0);
    % V322cr = VV3(deltacc,0);
    % contour(y1,y2,zz,[V3cr-V32cr+V322cr V3cr-V32cr+V322cr],'y-','linewidth',1.5,"ShowText",false);
    % 
    % contour(y1,y2,zz_2,[V3cr V3cr],'y-','linewidth',1.5,"ShowText",false);
    %contour(y1,y2,dzz,[-10 -5 0 5 10],'y:','linewidth',0.5,"ShowText",true);
 elseif system == "VFM2"
    % energy function 3
    syms deltax yx;
    delta_s = prefault_SEP(1);
    deltauep = ep_set(2).xep(1);
    %V3 = 1/ki/2*(omegax-kp*vvq)^2 - (Ug*cos(deltax)-Ug*cos(delta_s)+Xg*Id*(deltax-delta_s)+1/2*Id*Lg*omegax*(deltax-delta_s));%
    V3 = C_dc/2*Kip*Vdc_ref*yx^2 + C_dc/3*Kip*yx^3- Pin*(deltax-delta_s) + Rg*(Vvfm^2*(deltax-delta_s)-Vvfm*Ug*(sin(deltax)-sin(delta_s)))/(Rg^2+Xg^2)-Xg*Vvfm*Ug*(cos(deltax)-cos(delta_s))/(Rg^2+Xg^2) ;%- C_dc/2*Kpp*yx*(Rg*(Vvfm^2-Vvfm*Ug*cos(deltax))/(Rg^2+Xg^2)+Xg*Vvfm*Ug*sin(deltax)/(Rg^2+Xg^2) - Pin)/2;%
    V3_2 = V3-C_dc/2*Kpp*yx*(Rg*(Vvfm^2-Vvfm*Ug*cos(deltax))/(Rg^2+Xg^2)+Xg*Vvfm*Ug*sin(deltax)/(Rg^2+Xg^2) - Pin)/2;
    V3=vpa(V3);
    VV3=matlabFunction(V3);
    VV3_2=matlabFunction(V3_2);
    V3d = jacobian(V3);
    VV3d = matlabFunction(V3d);
    x1=-2*pi:0.01*pi:2*pi;
    x2=-4:0.02:6;%-8:0.1:8;
    [y1,y2]=meshgrid(x1,x2);
    zz = zeros(length(x2),length(x1));
    zz_2 = zeros(length(x2),length(x1));
    dzz = zeros(length(x2),length(x1));
    for a = 1: length(x1)
        for b = 1: length(x2)
            V = VV3(y1(b,a), y2(b,a));
            dV = VV3d(y1(b,a), y2(b,a))*f_VFM_normal([y1(b,a) y2(b,a)]);
            zz(b,a) = V;
            zz_2(b,a) = VV3_2(y1(b,a), y2(b,a));
            dzz(b,a)=dV;
        end
    end
    V3cr = VV3(ep_set(2).xep(1),0);
    contour(y1,y2,zz,[V3cr V3cr],'b-','linewidth',1.5,"ShowText",false);

elseif system == "GFM"
% syms deltax omegax lamda;
% 
% delta_s = prefault_SEP(1);
% deltauep = ep_set(2).xep(1);
% 
% J_ori = J/Ws;
% 
% V3 = J_ori/2*(omegax*Ws)^2 ...
%     - Pm*(deltax-delta_s) ...
%     + Rg*(Vgfm^2*(deltax-delta_s) ...
%     - Vgfm*Ug*(sin(deltax)-sin(delta_s)))/(Rg^2+Xg^2) ...
%     - Xg*Vgfm*Ug*(cos(deltax)-cos(delta_s))/(Rg^2+Xg^2) ...
%     + lamda*D*(omegax)*(deltax-delta_s) ...
%     + lamda/2*D^2/J/Ws*(deltax-delta_s)^2;
% 
% V3_2 = J_ori/2*(omegax*Ws)^2 ...
%     - Pm*(deltax-delta_s) ...
%     + Rg*(Vgfm^2*(deltax-delta_s) ...
%     - Vgfm*Ug*(sin(deltax)-sin(delta_s)))/(Rg^2+Xg^2) ...
%     - Xg*Vgfm*Ug*(cos(deltax)-cos(delta_s))/(Rg^2+Xg^2);
% 
% V3 = vpa(V3);
% 
% VV3 = matlabFunction(V3, 'Vars', [deltax, omegax, lamda]);
% VV3_2 = matlabFunction(V3_2, 'Vars', [deltax, omegax]);
% 
% x1 = -2*pi:0.01*pi:2*pi;
% x2 = -0.2:0.0005:0.2;
% 
% [y1,y2] = meshgrid(x1,x2);
% 
% lambda_list = 0:0.2:1;
% colors = lines(length(lambda_list));
% 
% delta_uep = ep_set(2).xep(1);
% 
% hold on
% 
% for k = 1:length(lambda_list)
% 
%     lam = lambda_list(k);
% 
%     zz = zeros(length(x2), length(x1));
% 
%     for a = 1:length(x1)
%         for b = 1:length(x2)
%             zz(b,a) = VV3(y1(b,a), y2(b,a), lam);
%         end
%     end
% 
%     V_line = zeros(size(x2));
% 
%     for b = 1:length(x2)
%         V_line(b) = VV3(delta_uep, x2(b), lam);
%     end
% 
%     V3cr = VV3(delta_uep, 0, lam);%min(V_line);%
% 
%     contour(y1, y2, zz, [V3cr V3cr], ...
%         'Color', colors(k,:), ...
%         'LineWidth', 1.5, ...
%         'ShowText', false);
% 
% end
    % energy function 4
    syms deltax omegax;

delta_s = prefault_SEP(1);
deltauep = ep_set(2).xep(1);

J_ori = J/Ws;

V3 = J_ori/2*(omegax*Ws)^2 ...
    - Pm*(deltax-delta_s) ...
    + Rg*(Vgfm^2*(deltax-delta_s) ...
    - Vgfm*Ug*(sin(deltax)-sin(delta_s)))/(Rg^2+Xg^2) ...
    - Xg*Vgfm*Ug*(cos(deltax)-cos(delta_s))/(Rg^2+Xg^2);% -(deltauep-deltax)/2*D*omegax;


V3 = vpa(V3);

VV3 = matlabFunction(V3, 'Vars', [deltax, omegax]);

x1 = -2*pi:0.01*pi:2*pi;
x2 = -0.2:0.0005:0.2;

[y1,y2] = meshgrid(x1,x2);


delta_uep = ep_set(2).xep(1);

hold on


zz = zeros(length(x2), length(x1));

for a = 1:length(x1)
    for b = 1:length(x2)
        zz(b,a) = VV3(y1(b,a), y2(b,a));
    end
end

    V3cr = VV3(delta_uep, 0);%min(V_line);%

    contour(y1, y2, zz, [V3cr V3cr], ...
        'Color', 'blue', ...
        'LineWidth', 1.5, ...
        'ShowText', false);




    % energy function 4
    syms deltax omegax;

delta_s = prefault_SEP(1);
deltauep = ep_set(2).xep(1);
pdp = Xg*Vgfm*Ug/(Xg^2+Rg^2)*cos(delta_s)+Rg*Vgfm*Ug/(Xg^2+Rg^2)*sin(delta_s);

PP = Rg*(Vgfm^2-Vgfm*Ug*cos(deltax))/(Rg^2+Xg^2)+Xg*Vgfm*Ug*sin(deltax)/(Rg^2+Xg^2);
J_ori = J/Ws;
beta = D*Ws/(D^2/4+J*Ws*pdp)*0.91;
V4 = J_ori/2*(omegax*Ws)^2 ...
    - Pm*(deltax-delta_s) ...
    + Rg*(Vgfm^2*(deltax-delta_s) ...
    - Vgfm*Ug*(sin(deltax)-sin(delta_s)))/(Rg^2+Xg^2) ...
    - Xg*Vgfm*Ug*(cos(deltax)-cos(delta_s))/(Rg^2+Xg^2) - beta*J*omegax*(Pm-PP);


V4 = vpa(V4);
V4d = jacobian(V4);
VV4 = matlabFunction(V4, 'Vars', [deltax, omegax]);
VV4d = matlabFunction(V4d);


x1 = -2*pi:0.01*pi:2*pi;
x2 = -0.2:0.0005:0.2;

[y1,y2] = meshgrid(x1,x2);


delta_uep = ep_set(2).xep(1);

hold on


zz = zeros(length(x2), length(x1));
dzz = zeros(length(x2),length(x1));

for a = 1:length(x1)
    for b = 1:length(x2)
        zz(b,a) = VV4(y1(b,a), y2(b,a));
        dzz(b,a)= VV4d(y1(b,a), y2(b,a))*f_GFM_normal([y1(b,a) y2(b,a)]);
    end
end

    V4cr = VV4(delta_uep, 0);%min(V_line);%

    contour(y1, y2, zz, [V4cr/2 V4cr], ...
        'Color', 'green', ...
        'LineWidth', 1.5, ...
        'ShowText', false);

   contour(y1,y2,dzz,[-1 -0.1 0 1],'g-','linewidth',1.5,"ShowText",true);


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
    case "GFM"
        dfdt = f_GFM_normal(x);
    case "VFM"
        dfdt = f_VFM_normal(x);%%
    case "VFM2"
        dfdt = f_VFM2_normal(x);%%
    case "GFMQ"
        dfdt = f_GFMQ_normal(x);
end
end
function dfdt = f_normal(x)
global system;
switch system 
    case "GFL"
        dfdt = f_GFL_normal(x);
    case "GFM"
        dfdt = f_GFM_normal(x);
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
        case "GFM"
          dfdt = f_GFM_fault(x);
        case "VFM"
          dfdt = f_VFM_fault_cl(x);%
    end
end
function dfdt = f_fault_g(t,x)
global system;
    switch system 
        case "GFM"
          dfdt = f_GFM_fault2(x);
        case "GFL"
          dfdt = f_GFL_fault_ful(x);
    end
end
function dfdt = f_fault_rp(t,x)
global system;
    switch system 
        case "GFM"
          dfdt = f_GFM_fault_rp(x);
    end
end
function dfdt = f_fault_gy(t,x)
global system;
    switch system 
        case "GFM"
          dfdt = f_GFM_fault_gy(x);
    end
end
function dfdt = f_fault_1(t,x)
global system;
    switch system 
        case "GFM"
          dfdt = f_GFM_fault_1(x);
    end
end
function dfdt = f_post(t,x)
global system;
    switch system 
        case "GFL"
          dfdt = f_GFL_normal(x);
        case "GFM"
          dfdt = f_GFM_normal(x);
        case "VFM"
          dfdt = f_VFM_normal(x); %
    end
end
function dfdt = f_post_g(t,x)
global system;
    switch system 
        case "GFM"
          dfdt = f_GFM_normal2(x);
        case "GFL"
          dfdt = f_GFL_normal_ful(x);
    end
end
function dfdt = f_post_rp(t,x)
global system;
    switch system 
        case "GFM"
          dfdt = f_GFM_normal_rp(x);
    end
end
function dfdt = f_post_1(t,x)
global system;
    switch system 
        case "GFM"
          dfdt = f_GFM_normal_1(x);
    end
end
function dfdt = f_post_gy(t,x)
global system;
    switch system 
        case "GFM"
          dfdt = f_GFM_normal_gy(x);
    end
end
function dfdt = f_prefault(t,x)
global system;
    switch system 
        case "GFL"
          dfdt = f_GFL_prefault(x);
        case "GFM"
          dfdt = f_GFM_prefault(x);
        case "VFM"
          dfdt = f_VFM_prefault(x);
        case "VFM2"
          dfdt = f_VFM_prefault(x);
        case "GFMQ"
          dfdt = f_GFMQ_normal(x);
    end
end