%% parameter


t_end = 0.5;

%grid

Xg = 0.5;
Rg = 0.05;%0.08
Ug = 1;
Ws = 2*pi*50; 
Lg= Xg/Ws;
W_g = 0;


%GFM
kgfm = 20*2*pi; 
m_gfm = kgfm/Ws; 
w_droop = 0.5*2*pi;
D = 10; %1/m_gfm;
J = 2; 
Vgfm = 1;
Pm = 1;


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



%system
global system;
global fault_type; %line_cut voltage_sag frequency
global limit_type
fault_type = "voltage_sag"; %"voltage_sag";%"line_cut";%"line_cut";
limit_type = "VA";   %"cir"   "VA"
system = "VFM";  %GFMQ   VOC
model = "original";% "original"

switch fault_type
    case "voltage_sag"
        %fault sag
        Ug_fault = 0.1;%0.1;
        X1 = 0.1;
        R1 = 0.01;
        Xgg = (Xg - X1)*2;
        Rgg = (Rg - R1)*2;
        t_c = 0.035;%0.072;%0.086;0.0795
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
    system = "GFM";
end
%% SEP

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
            elseif flag ==1
               postfault_UEP=xep;
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

if system == "VFM"
    ymin = -4;
    ymax = 6;
elseif system == "GFM"
ymin=-0.2;
ymax=0.2;
end
color_code = {'black','magenta','red','black'};

axis([-1*pi,3/2*pi,ymin,ymax]);
xticks(-2*pi:pi/2:2*pi);
xticklabels({'$-2\pi$', '', '$-\pi$', '','$0$', '','$\pi$', '','$2\pi$'});
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 14);
for mm = 1 : length(ep_set_ext)
    xep = ep_set_ext(mm).xep;
    flag= ep_set_ext(mm).flag;
    scatter(xep(1),xep(2),color_code{flag+1},'LineWidth', 1.5);
    if flag == 1
        v = ep_set_ext(mm).v;
        vv = ep_set_ext(mm).vv;
        perturb = 1e-3;
        switch system
            case "GFM"
            [~ , x_p] = ode78(@f_backward,[0,1],xep+v*perturb,odeset('RelTol',1e-5));
            [~ , x_n] = ode78(@f_backward,[0,1],xep-v*perturb,odeset('RelTol',1e-5)); 
            case "VFM"
            [~ , x_p] = ode78(@f_backward,[0,2],xep+v*perturb,odeset('RelTol',1e-5));
            [~ , x_n] = ode78(@f_backward,[0,2],xep-v*perturb,odeset('RelTol',1e-5)); 
        end
        x_all = [flip(x_n,1);x_p];
        plot(x_all(:,1),x_all(:,2),'k:','linewidth',1.5);
    end
end

%%
switch fault_type
    case "voltage_sag"
    
    if system == "GFM"
    t_start = 0.1;
    delta_pre = [prefault_SEP(1); prefault_SEP(1)];
    omega_pre = [prefault_SEP(2); prefault_SEP(2)];

    t_prefault = [0;0.1];

    [t_fault , x_all] = ode78(@f_fault,[t_start,t_start+t_c],[prefault_SEP(1);prefault_SEP(2)],odeset('RelTol',1e-6));
    omega_fault= x_all(:,2);
    delta_fault =  x_all(:,1);
    [t_postfault , x_all2] = ode78(@f_post,[t_fault(end),t_end],x_all(end,1:2),odeset('RelTol',1e-10));
    omega_post= x_all2(:,2);
    delta_post = x_all2(:,1);



    figure(f1)

    plot(delta_fault,omega_fault,'r-','linewidth',1.5);
    plot(delta_post(1),omega_post(1),'k.','MarkerSize',15);
    plot(delta_fault(1),omega_fault(1),'k.','MarkerSize',15);
    plot(delta_post,omega_post,'b-','linewidth',1.5)

    
    elseif system == "VFM"
        t_start = 0.1;
        delta_pre = [prefault_SEP(1); prefault_SEP(1)];
        omega_pre = [prefault_SEP(2); prefault_SEP(2)];
    
        t_prefault = [0;0.1];
    
        [t_fault , x_all] = ode78(@f_fault,[t_start,t_start+t_c],[prefault_SEP(1);prefault_SEP(2)],odeset('RelTol',1e-6));
        y_fault= x_all(:,2);
        delta_fault =  x_all(:,1);
        options = odeset('RelTol',1e-10);
        [t_postfault , x_all2] = ode78(@f_post,[t_fault(end),t_end],x_all(end,1:2),options);
        y_post= x_all2(:,2);
        delta_post = x_all2(:,1);
    
    
    
        figure(f1)
    
        plot(delta_fault,y_fault,'r:','linewidth',1.5);
        plot(delta_post(1),y_post(1),'k.','MarkerSize',15);
        plot(delta_fault(1),y_fault(1),'k.','MarkerSize',15);
        plot(delta_post,y_post,'b:','linewidth',1.5)



    end
 end
%% current limit
if limit_type == "cir"

%%  ======= circle limit ==========
% ===== 1. 找 equilibrium =====
torralence = 1e-2;
mm = 1;
ep_set_cl = [];

options = optimoptions('fsolve',...
    'FunctionTolerance',1e-10,...
    'MaxIterations',100000,...
    'OptimalityTolerance',1e-10,...
    'Display','off');

for n = 1:length(x_set(1,:))
    xep0 = x_set(:,n);

    [xep,ferr,exitflag] = fsolve(@f_VFM_normal_cl_circle, xep0, options);

    if exitflag > 0 && maxabs(ferr) < torralence
        if isnewxep(ep_set_cl,xep,torralence)

            % 数值 Jacobian
            A = numerical_jacobian(@f_VFM_normal_cl_circle, xep);

            [V,Lambda] = eig(A);
            Lambda = diag(Lambda);

            sig = sign(sign(real(Lambda))+0.1);
            sig = (sig + 1)/2;
            flag = sum(sig);

            v = V(:,~sig);

            ep_set_cl(mm).xep = xep;
            ep_set_cl(mm).A = A;
            ep_set_cl(mm).Lambda = Lambda;
            ep_set_cl(mm).V = V;
            ep_set_cl(mm).v = v;
            ep_set_cl(mm).flag = flag;

            mm = mm + 1;
        end
    end
end

% ===== 2. 扩展周期（和原代码一致） =====
clear ep_set_ext_cl
for n = 1:length(ep_set_cl)
    ep_set_ext_cl(n) = ep_set_cl(n);
end

% ===== 3. 画稳定流形 =====
figure(f1)
hold on

for mm = 1:length(ep_set_ext_cl)

    xep  = ep_set_ext_cl(mm).xep;
    flag = ep_set_ext_cl(mm).flag;

    if flag == 1   % UEP

        % 取稳定特征向量
        stable_idx = find(real(ep_set_ext_cl(mm).Lambda) < 0);
        v = ep_set_ext_cl(mm).V(:, stable_idx);

        if size(v,2) > 1
            v = v(1,:); % 简化（理论上2维只会1个）
        end

        v = real(v);
        v = v / norm(v);

        perturb = 1e-4;

        opt_back = odeset('RelTol',1e-6);

        [~, x_p] = ode78(@(t,x)-f_VFM_normal_cl_circle(x), [0,2], xep + perturb*v);
        [~, x_n] = ode78(@(t,x)-f_VFM_normal_cl_circle(x), [0,2], xep - perturb*v);

        x_all = [flip(x_n,1); x_p];

        plot(x_all(:,1), x_all(:,2),'k-','LineWidth',1.5);
    end
end

% ===== 4. 画限流分界 δc =====
Ug   = evalin('base','Ug');
Vvfm = evalin('base','Vvfm');
Ilim = evalin('base','Ilim');
Xg   = evalin('base','Xg');
Rg   = evalin('base','Rg');

deltac = acos((Vvfm^2 + Ug^2 - Ilim^2*(Xg^2 + Rg^2))/(2*Vvfm*Ug));

yl = ylim;

plot([deltac deltac], [yl(1) yl(2)], 'g-','LineWidth',2);
plot([-deltac -deltac], yl, 'g-','LineWidth',2);

delta_uep_cl = [];

for k = 1:length(ep_set_cl)
    if ep_set_cl(k).flag == 1   % UEP
        delta_uep_cl = ep_set_cl(k).xep(1);
        break;   % 如果只有一个UEP，直接取第一个
    end
end

if ~isempty(delta_uep_cl)
    yl = ylim;
    plot([delta_uep_cl delta_uep_cl], [yl(1) yl(2)], 'c-','LineWidth',2);
end


% =====  traj ======%
t_start = 0.1;
delta_pre_cl = [prefault_SEP(1); prefault_SEP(1)];
y_pre_cl     = [prefault_SEP(2); prefault_SEP(2)];

t_prefault_cl = [0; 0.1];

[t_fault_cl, x_fault_cl] = ode78(@(t,x) f_VFM_fault_cl_circle(x), ...
    [t_start, t_start + t_c], ...
    [prefault_SEP(1); prefault_SEP(2)], ...
    odeset('RelTol',1e-6));

delta_fault_cl = x_fault_cl(:,1);
y_fault_cl     = x_fault_cl(:,2);

options_cl = odeset('RelTol',1e-10);
[t_postfault_cl, x_post_cl] = ode78(@(t,x) f_VFM_normal_cl_circle(x), ...
    [t_fault_cl(end), t_end], ...
    x_fault_cl(end,1:2), ...
    options_cl);

delta_post_cl = x_post_cl(:,1);
y_post_cl     = x_post_cl(:,2);

figure(f1)
hold on

plot(delta_fault_cl, y_fault_cl, 'r-', 'LineWidth', 1.8);

plot(delta_post_cl(1), y_post_cl(1), 'k.', 'MarkerSize', 6, 'LineWidth', 1.2);
plot(delta_fault_cl(1), y_fault_cl(1), 'k.', 'MarkerSize', 6, 'LineWidth', 1.2);

% 故障后轨迹
plot(delta_post_cl, y_post_cl, 'b-', 'LineWidth', 1.8);

%%  =============virtuial admitance ===================
elseif limit_type == "VA"
%% ===== 1. 找 equilibrium =====
torralence = 1e-2;
mm = 1;
ep_set_va = [];

options = optimoptions('fsolve',...
    'FunctionTolerance',1e-10,...
    'MaxIterations',100000,...
    'OptimalityTolerance',1e-10,...
    'Display','off');

for n = 1:length(x_set(1,:))
    xep0 = x_set(:,n);

    % ★★★ 用 VA 模型 ★★★
    [xep,ferr,exitflag] = fsolve(@f_VFM_normal_cl_va, xep0, options);

    if exitflag > 0 && maxabs(ferr) < torralence
        if isnewxep(ep_set_va,xep,torralence)

            A = numerical_jacobian(@f_VFM_normal_cl_va, xep);

            [V,Lambda] = eig(A);
            Lambda = diag(Lambda);

            sig = sign(sign(real(Lambda))+0.1);
            sig = (sig + 1)/2;
            flag = sum(sig);

            v = V(:,~sig);

            ep_set_va(mm).xep = xep;
            ep_set_va(mm).A = A;
            ep_set_va(mm).Lambda = Lambda;
            ep_set_va(mm).V = V;
            ep_set_va(mm).v = v;
            ep_set_va(mm).flag = flag;

            mm = mm + 1;
        end
    end
end

% ===== 2. 扩展周期（简单版） =====
clear ep_set_ext_va
for n = 1:length(ep_set_va)
    ep_set_ext_va(n) = ep_set_va(n);
end

% ===== 3. 画稳定流形 =====
figure(f1)
hold on

for mm = 1:length(ep_set_ext_va)

    xep  = ep_set_ext_va(mm).xep;
    flag = ep_set_ext_va(mm).flag;

    if flag == 1   % UEP

        stable_idx = find(real(ep_set_ext_va(mm).Lambda) < 0);
        v = ep_set_ext_va(mm).V(:, stable_idx);

        if size(v,2) > 1
            v = v(1,:);
        end

        v = v / norm(v);

        perturb = 1e-4;

        [~, x_p] = ode78(@(t,x)-f_VFM_normal_cl_va(x), [0,2], xep + perturb*v);
        [~, x_n] = ode78(@(t,x)-f_VFM_normal_cl_va(x), [0,2], xep - perturb*v);

        x_all = [flip(x_n,1); x_p];

        plot(x_all(:,1), x_all(:,2),'k-','LineWidth',1.5);
    end
end

% ===== 4. 找 VA 模型 UEP 并画竖线 =====
delta_uep_va = [];

for k = 1:length(ep_set_va)
    if ep_set_va(k).flag == 1
        delta_uep_va = ep_set_va(k).xep(1);
        break;
    end
end

if ~isempty(delta_uep_va)
    yl = ylim;
    plot([delta_uep_va delta_uep_va], [yl(1) yl(2)], 'm-','LineWidth',2);
end
deltac = acos((Vvfm^2 + Ug^2 - Ilim^2*(Xg^2 + Rg^2))/(2*Vvfm*Ug));
plot([deltac deltac], [yl(1) yl(2)], 'g-','LineWidth',2);
plot([-deltac -deltac], yl, 'g-','LineWidth',2);
% ===== 5. fault / postfault trajectory =====
t_start = 0.1;

[t_fault_va, x_fault_va] = ode78(@(t,x) f_VFM_fault_cl_va(x), ...
    [t_start, t_start + t_c], ...
    [prefault_SEP(1); prefault_SEP(2)], ...
    odeset('RelTol',1e-6));

delta_fault_va = x_fault_va(:,1);
y_fault_va     = x_fault_va(:,2);

options_va = odeset('RelTol',1e-10);
[t_post_va, x_post_va] = ode78(@(t,x) f_VFM_normal_cl_va(x), ...
    [t_fault_va(end), t_end], ...
    x_fault_va(end,1:2), ...
    options_va);

delta_post_va = x_post_va(:,1);
y_post_va     = x_post_va(:,2);

figure(f1)
hold on

plot(delta_fault_va, y_fault_va, 'r-', 'LineWidth', 1.8);

plot(delta_post_va(1), y_post_va(1), 'k.', 'MarkerSize', 6);
plot(delta_fault_va(1), y_fault_va(1), 'k.', 'MarkerSize', 6);

plot(delta_post_va, y_post_va, 'b-', 'LineWidth', 1.8);

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
    case "GFM"
        dfdt = f_GFM_normal(x);
    case "VFM"
        dfdt = f_VFM_normal(x);
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
        case "GFM"
          dfdt = f_GFM_fault(x);
        case "VFM"
          dfdt = f_VFM_fault(x);
    end
end
function dfdt = f_post(t,x)
global system;
    switch system 
        case "GFM"
          dfdt = f_GFM_normal(x);
        case "VFM"
          dfdt = f_VFM_normal(x);
    end
end

function dfdt = f_prefault(t,x)
global system;
    switch system 
        case "GFM"
          dfdt = f_GFM_normal(x);
        case "VFM"
          dfdt = f_VFM_normal(x);
    end
end


function A = numerical_jacobian(fun, x)
    n = length(x);
    A = zeros(n,n);
    h = 1e-6;

    fx = fun(x);

    for k = 1:n
        xp = x;
        xm = x;

        dx = h * max(1, abs(x(k)));
        xp(k) = xp(k) + dx;
        xm(k) = xm(k) - dx;

        fp = fun(xp);
        fm = fun(xm);

        A(:,k) = (fp - fm) / (2*dx);
    end
end