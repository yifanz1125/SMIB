%% parameter

%grid
Rg = 0;
Xg = 0.4;
Ug = 0.9;
Ws = 2*pi*50; 
Lg= Xg/Ws;
Wbase = 50*2*pi;

%fault
Ug_fault = 0.1;
t_c = 0.0271;

%GFL
%Id = 1;
Iq = 0; %negative q
w_pll = 10*2*pi;     % (rad/s)
kp = w_pll;
ki = kp*400;%w_pll^2/2;
w_limit = 500;

Vdc_ref = 2.5;
C_dc = 2.25e2;
w_vdc   = 1*2*pi;
kp_v_dc	= Vdc_ref*C_dc*w_vdc/Wbase;
ki_v_dc	= kp_v_dc*w_vdc^2/2;%kp_v_dc*w_vdc^2*3;
Pref = 1;


%system
global system;  
system = "DVCslow";
model = "original";% "original"


%% 
try
    system;
catch
    system = "DVCslow";
end
%% SEP
switch system
    case "DVCslow"
        x=0:0.02:5;
        n = length(x);
        x_set = zeros(2,n);
        x_set(1,:) = x;
    case"PLLslow"
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
            
            ep_set(mm).xep = xep; %#ok<*SAGROW> 
            ep_set(mm).A = A;
            ep_set(mm).Lambda = Lambda;
            ep_set(mm).V = V;   
            ep_set(mm).v = v;     % stable eigenvectors of unstable ep 
            ep_set(mm).flag = flag;
           
            mm = mm+1;
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
switch system
    case "DVCslow"
            ep_set_ext = ep_set; 
    case "PLLslow"
        for n = 1:length(ep_set)
            mm = (n-1)*3;
            ep_set_ext(mm+1)=ep_set(n); %#ok<*AGROW> 
            ep_set_ext(mm+2)=ep_set(n);
            ep_set_ext(mm+3)=ep_set(n); 
            ep_set_ext(mm+2).xep(1) = ep_set(n).xep(1) - 2*pi;
            ep_set_ext(mm+3).xep(1) = ep_set(n).xep(1) + 2*pi;
        end
end

figure;
hold on;
grid on;
ymin=-0.05;
ymax=0.05;
color_code = {'blue','magenta','red','black'};
axis([0,2.5,ymin,ymax]);

for mm = 1 : length(ep_set_ext)
    xep = ep_set_ext(mm).xep;
    flag= ep_set_ext(mm).flag;
    scatter(xep(1),xep(2),color_code{flag+1});
    if flag == 1
        v = ep_set_ext(mm).v;
        perturb = 1e-3;
        switch system
            case "DVCslow"
            [~ , x_p] = ode78(@f_backward,[0,1],xep+v*perturb,odeset('RelTol',1e-5));
            [~ , x_n] = ode78(@f_backward,[0,0.8],xep-v*perturb,odeset('RelTol',1e-5));
        end
        x_all = [flip(x_n,1);x_p];
        plot(x_all(:,1),x_all(:,2),'k-','linewidth',1.5);
    end
end

 %trajectory drawing
[t1 , x_all] = ode78(@f_forward,[0,1],[1.2;0],odeset('RelTol',1e-5));
% Vq_fault = (Xg*Id+Rg*Iq-Ug_fault*sin(x_all(:,1))+Id*Lg*x_all(:,2))./(1-Id*Lg*kp);
% omega_fault=kp.*Vq_fault+x_all(:,2);
% omega_fault(find(omega_fault>w_limit))=w_limit;
% [t2 , x_all2] = ode78(@f_post,[t1(end),0.4],x_all(end,:),odeset('RelTol',1e-5));
% Vq_post = (Xg*Id+Rg*Iq-Ug*sin(x_all2(:,1))+Id*Lg*x_all2(:,2))./(1-Id*Lg*kp);
% omega_post=kp.*Vq_post+x_all2(:,2);
% omega_post(find(omega_post>w_limit))=w_limit;
plot(x_all(:,1),x_all(:,2),'k-.','linewidth',1.5);
% plot(x_all2(1,1),omega_post(1),'k.','MarkerSize',15);
% plot(x_all(1,1),omega_fault(1),'k.','MarkerSize',15);
% plot(x_all(end,1),omega_fault(end),'k.','MarkerSize',15);
% plot(x_all2(:,1),omega_post,'k-.','linewidth',1.5);

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
%         plot(x_all(:,1),x_all(:,2),'k:','linewidth',1);%scatter(x_all(:,1),x_all(:,2),'.');
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
% model = "original";
%%
% energy function 1
syms yx Idx;
V1 = (C_dc/Ws)/4*ki_v_dc*yx^2-Pref*Idx-(sqrt(Ug^2-Idx^2*Xg^2))^3/3/Xg^2;
V1=vpa(V1);
VV1=matlabFunction(V1);
V1d = jacobian(V1);
VV1d = matlabFunction(V1d);

x1=0:0.01:2.5;
x2=ymin:0.001:ymax;
[y1,y2]=meshgrid(x1,x2);
zz = zeros(length(x2),length(x1));
dzz = zeros(length(x2),length(x1));
for a = 1: length(x1)
    for b = 1: length(x2)
        V = VV1(y1(b,a), y2(b,a));
        dV = VV1d(y1(b,a), y2(b,a))*f_DVCslow([y1(b,a) y2(b,a)]);
        zz(b,a) = real(V);
        dzz(b,a)=dV;
    end
end
Vcr1 = VV1(ep_set(2).xep(1),ep_set(2).xep(2));
contour(y1,y2,zz,[Vcr1 Vcr1],'b-','linewidth',1,"ShowText",false);
%contour(y1,y2,dzz,[-0.5 -0.1 0 0.1],'r:','linewidth',0.5,"ShowText",true);



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
    case "DVCslow"
        dfdt = f_DVCslow(x);
end
end
function dfdt = f_normal(x)
global system;
switch system 
    case "DVCslow"
        dfdt = f_DVCslow(x);
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

function dfdt = f_fault(t,x)
global system;
    switch system 
        case "DVCslow"
          dfdt = f_DVCslow(x);
    end
end
function dfdt = f_post(t,x)
global system;
    switch system 
        case "DVCslow"
          dfdt = f_DVCslow(x);
    end
end