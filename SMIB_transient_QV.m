%% parameter


t_end = 0.5;

%grid

Xg = 0.5;
Rg = 0.015;%0.08
Ug = 1;
Ws = 2*pi*50; 
Lg= Xg/Ws;
W_g = 0;


%GFM with Q-V droop
m_gfm = 0.05;
Vgfm = 1;
Pm = 1;
k_q = 1;
tau_q = 1/(0.5*2*pi);
fsep = @(delta) Rg*(Vgfm^2 - Vgfm*Ug*cos(delta))/(Rg^2 + Xg^2) ...
           + Xg*Vgfm*Ug*sin(delta)/(Rg^2 + Xg^2) ...
           - Pm;
deltas = fsolve(fsep, 0);
Qref = Xg*(Vgfm^2 - Vgfm*Ug*cos(deltas))/(Rg^2+Xg^2) - Rg*Vgfm*Ug*sin(deltas)/(Rg^2+Xg^2);


%GFM with VOC
kv = 1; 
ki = 1; 
C = 0.25; %0.2679
xi = 20;
VN = 1; 
Pref = 1;
fsep = @(delta) Rg*(VN^2 - VN*Ug*cos(delta))/(Rg^2 + Xg^2) ...
           + Xg*VN*Ug*sin(delta)/(Rg^2 + Xg^2) ...
           - Pref;
deltas = fsolve(fsep, 0);
Qref2 = Xg*(VN^2 - VN*Ug*cos(deltas))/(Rg^2+Xg^2) - Rg*VN*Ug*sin(deltas)/(Rg^2+Xg^2);


%system
global system;
global fault_type; %line_cut voltage_sag frequency
fault_type = "voltage_sag"; %"voltage_sag";%"line_cut";%"line_cut";
system = "VOC";  %GFMQ   VOC  "dVOC"
model = "original";% "original"

switch fault_type
    case "voltage_sag"
        %fault sag
        Ug_fault = 0.1;%0.1;
        X1 = 0.1;
        R1 = 0.01;
        Xgg = (Xg - X1)*2;
        Rgg = (Rg - R1)*2;
        t_c = 0.06;%0.0545
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
    system = "GFMQ";
end
%% SEP

x=(0:0.1:1)*2*pi;
y=(0:0.05:1)*2;
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

ymin = 0;
ymax = 1.5;
color_code = {'black','magenta','red','black'};

axis([-1/4*pi,1*pi,ymin,ymax]);
xticks(-2*pi:pi/2:2*pi);
xticklabels({'$-2\pi$', '', '$-\pi$', '','$0$', '','$\pi$', '','$2\pi$'});
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 18);
if system == "GFL"   
rangex=[-acos(Id*Lg*ki/Ug/kp),acos(Id*Lg*ki/Ug/kp),acos(Id*Lg*ki/Ug/kp),-acos(Id*Lg*ki/Ug/kp)];  rangey=[ymin,ymin,ymax,ymax];
fill(rangex,rangey,[.9 .9 .9], 'linestyle', 'none', 'FaceAlpha',0.6);
end
for mm = 1 : length(ep_set_ext)
    xep = ep_set_ext(mm).xep;
    flag= ep_set_ext(mm).flag;
    plot(xep(1), xep(2), 'o', ...
    'Color', color_code{flag+1}, ...
    'LineWidth', 2.5, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'none');
    if flag == 1
        v = ep_set_ext(mm).v;
        vv = ep_set_ext(mm).vv;
        perturb = 1e-3;
        switch system
            case "GFMQ"
            [~ , x_p] = ode78(@f_backward,[0,1],xep+v*perturb,odeset('RelTol',1e-5));
            [~ , x_n] = ode78(@f_backward,[0,1],xep-v*perturb,odeset('RelTol',1e-5)); 
            case "VOC"
            [~ , x_p] = ode78(@f_backward,[0,0.5],xep+v*perturb,odeset('RelTol',1e-5));
            [~ , x_n] = ode78(@f_backward,[0,0.5],xep-v*perturb,odeset('RelTol',1e-5)); 
        end
        x_all = [flip(x_n,1);x_p];
        %x_allall = [flip(x_nn,1);x_pp];
        plot(x_all(:,1),x_all(:,2),'k-','linewidth',1.5);%scatter(x_all(:,1),x_all(:,2),'.');
        %plot(x_allall(:,1),x_allall(:,2),'r-','linewidth',1.5);%scatter(x_all(:,1),x_all(:,2),'.');
        % p_traj= Rg*(Vvfm^2-Vvfm*Ug*cos(x_all(:,1)))/(Rg^2+Xg^2)+Xg*Vvfm*Ug*sin(x_all(:,1))/(Rg^2+Xg^2);
        % plot(p_traj-Pin,x_all(:,2),'y-','linewidth',1.5);
    end
end

%%
switch fault_type
    case "voltage_sag"
    
    if system == "GFMQ"
    t_start = 0.1;
    delta_pre = [prefault_SEP(1); prefault_SEP(1)];
    voltage_pre = [prefault_SEP(2); prefault_SEP(2)];

    t_prefault = [0;0.1];

    [t_fault , x_all] = ode78(@f_fault,[t_start,t_start+t_c],[prefault_SEP(1);prefault_SEP(2)],odeset('RelTol',1e-6));
    voltage_fault= x_all(:,2);
    delta_fault =  x_all(:,1);
    [t_postfault , x_all2] = ode78(@f_post,[t_fault(end),t_end],x_all(end,1:2),odeset('RelTol',1e-10));
    voltage_post= x_all2(:,2);
    delta_post = x_all2(:,1);



    figure(f1)

    plot(delta_fault,voltage_fault,'r-','linewidth',1.5);
    plot(delta_post(1),voltage_post(1),'k.','MarkerSize',15);
    plot(delta_fault(1),voltage_fault(1),'k.','MarkerSize',15);
    plot(delta_post,voltage_post,'b-','linewidth',1.5)

    %[t_postfault , test] = ode78(@f_post,[t_fault(end),t_end],x_all(end,1:2),odeset('RelTol',1e-10));
    elseif system == "VOC"
        % t_start = 0.1;
        % delta_pre = [prefault_SEP(1); prefault_SEP(1)];
        % voltage_pre = [prefault_SEP(2); prefault_SEP(2)];
        % 
        % t_prefault = [0;0.1];
        % 
        % [t_fault , x_all] = ode78(@f_fault,[t_start,t_start+t_c],[prefault_SEP(1);prefault_SEP(2)],odeset('RelTol',1e-6));
        % voltage_fault= x_all(:,2);
        % delta_fault =  x_all(:,1);
        % options = odeset('RelTol',1e-10,'Events',@stop_event);
        % [t_postfault , x_all2] = ode78(@f_post,...
        %     [t_fault(end),t_end],...
        %     x_all(end,1:2),...
        %     options);
        % voltage_post= x_all2(:,2);
        % delta_post = x_all2(:,1);
        % 
        % 
        % 
        % figure(f1)
        % 
        % plot(delta_fault,voltage_fault,'r:','linewidth',1.5);
        % plot(delta_post(1),voltage_post(1),'k.','MarkerSize',15);
        % plot(delta_fault(1),voltage_fault(1),'k.','MarkerSize',15);
        % plot(delta_post,voltage_post,'b:','linewidth',1.5)



    end



 end
%% energy function 
if system == "VOC" 
% 
delta_sep = postfault_SEP(1);
voltage_sep = postfault_SEP(2);
syms deltax voltagex;
V1 = -Pref*deltax-Xg/(Xg^2+Rg^2)*voltagex*Ug*cos(deltax)-Rg/(Xg^2+Rg^2)*voltagex*Ug*sin(deltax);
V2 = 1/2*Xg/(Xg^2+Rg^2)*voltagex^2-Qref2*log(voltagex)-C*xi/(3*kv^3*ki)*(VN^2*voltagex^2-1/2*voltagex^4);
V3 = 1/2*Rg/(Xg^2+Rg^2)*(deltax-delta_sep)*(voltage_sep^2+voltage_sep*voltagex+voltagex^2);
Vt= V1 +V2 +V3;

Vt=vpa(Vt);
VV=matlabFunction(Vt);
Vd = jacobian(Vt);
VVd = matlabFunction(Vd);

x1=-2*pi:0.02*pi:2*pi;
x2=ymin+1e-6:1e-2:ymax;
[y1,y2]=meshgrid(x1,x2);
zz = zeros(length(x2),length(x1));
dzz = zeros(length(x2),length(x1));
for a = 1: length(x1)
    for b = 1: length(x2)
        V = VV(y1(b,a), y2(b,a));
        dV = VVd(y1(b,a), y2(b,a))*f_VOC_normal([y1(b,a) y2(b,a)]);
        zz(b,a) = V;
        dzz(b,a)=dV;
    end
end
Vcr1 = VV(postfault_UEP(1),postfault_UEP(2));
contour(y1,y2,zz,[Vcr1 Vcr1],'m-','linewidth',1,"ShowText",false);
%contour(y1,y2,dzz,[-20 -10 -1 0 1],'y-','linewidth',0.5,"ShowText",true);




%% current limit

% Ilim = 2;
% 
% figure(f1);   % 调出已有图框
% hold on;
% 
% % 定义隐函数
% f_cl = @(deltax, voltagex) ...
%     sqrt((voltagex.^2 + Ug^2 - 2*voltagex.*Ug.*cos(deltax)) / (Xg^2 + Rg^2)) - Ilim;
% 
% % 绘制等值线（=0 对应电流等于 Ilim）
% fimplicit(f_cl, [-2*pi, 2*pi, ymin, ymax],'LineWidth', 1.8, 'Color', 'g');
% 
% % （可选）加图例
% % legend('SEP/轨迹','current limit');
% 
% %% 符号变量
% syms deltax voltagex;
% 
% % 公共分母
% Den = voltagex^2 + Ug^2 - 2*voltagex*Ug*cos(deltax);
% 
% % 限流下的功率表达式
% Pi = sqrt(Den/Ilim^2 - Xg^2) ...
%    ./ Den .* Ilim^2 .* (voltagex*Ug*cos(deltax)-Ug^2) ...
%    + Xg ./ Den .* Ilim^2 .* voltagex*Ug.*sin(deltax) ...
%    + Ilim^2*Rg;
% 
% % 转成数值函数
% Pi_fun = matlabFunction(Pi, 'Vars', [deltax, voltagex]);
% 
% % 网格范围
% x1 = -2*pi:0.01*pi:2*pi;      % deltax
% x2 = ymin+1e-4:1e-4:ymax;     % voltagex，避免分母/开方奇异
% [DX, VX] = meshgrid(x1, x2);
% 
% % 数值计算
% Den_num = VX.^2 + Ug^2 - 2*VX.*Ug.*cos(DX);
% rad_num = Den_num/Ilim^2 - Xg^2;
% 
% Pi_val = nan(size(DX));
% 
% % % 只在实数有定义区域计算
% % mask = (Den_num > 1e-10) & (rad_num >= 0);
% % Pi_val(mask) = sqrt(rad_num(mask)) ...
% %     ./ Den_num(mask) .* Ilim^2 .* (VX(mask).*Ug.*cos(DX(mask)) - Ug^2) ...
% %     + Xg ./ Den_num(mask) .* Ilim^2 .* VX(mask).*Ug.*sin(DX(mask)) ...
% %     + Ilim^2*Rg;
% 
% % 数值计算
% Den_num = VX.^2 + Ug^2 - 2*VX.*Ug.*cos(DX);
% rad_num = Den_num/Ilim^2 - Xg^2;
% Imag_num = sqrt(Den_num./(Xg^2 + Rg^2));
% 
% Pi_val = nan(size(DX));
% 
% % 只画 Imag >= Ilim 的区域
% mask = (Den_num > 1e-10) & (rad_num >= 0) & (Imag_num >= Ilim);
% 
% Pi_val(mask) = sqrt(rad_num(mask)) ...
%     ./ Den_num(mask) .* Ilim^2 .* (VX(mask).*Ug.*cos(DX(mask)) - Ug^2) ...
%     + Xg ./ Den_num(mask) .* Ilim^2 .* VX(mask).*Ug.*sin(DX(mask)) ...
%     + Ilim^2*Rg;
% 
% 
% % 新开图窗：三维图
% f2 = figure(2);
% surf(DX, VX, Pi_val, 'EdgeColor', 'none');
% grid on;
% xlabel('\delta_x');
% ylabel('voltage_x');
% zlabel('P_i');
% title('Power under current limit');
% colorbar;
% view(45,30);
% 
% % 可选：如果想让图更平滑一点
% shading interp;
% 
% %% Q under current limit (3D plot)
% 
% syms deltax voltagex;
% 
% % 公共项
% Den = voltagex^2 + Ug^2 - 2*voltagex*Ug*cos(deltax);
% Rad = Den/Ilim^2 - Xg^2;
% 
% % 图中公式：Q(E, delta)
% Qilim = Xg*Ilim^2 ./ Den .* (voltagex^2 - voltagex*Ug*cos(deltax)) ...
%       - Ilim^2 .* sqrt(Rad) ./ Den .* voltagex*Ug.*sin(deltax);
% 
% % 转数值函数（如果后面还想调用）
% Qilim_fun = matlabFunction(Qilim, 'Vars', [deltax, voltagex]);
% 
% % 网格
% x1 = -2*pi:0.02*pi:2*pi;      % deltax
% x2 = ymin+1e-4:1e-2:ymax;     % voltagex
% [DXq, VXq] = meshgrid(x1, x2);
% 
% % % 数值计算
% % Den_q = VXq.^2 + Ug^2 - 2*VXq.*Ug.*cos(DXq);
% % Rad_q = Den_q/Ilim^2 - Xg^2;
% % 
% % Q_val = nan(size(DXq));
% % 
% % % 只在有定义的区域计算
% % mask_q = (Den_q > 1e-8) & (Rad_q >= 0);
% % 
% % Q_val(mask_q) = ...
% %     Xg*Ilim^2 ./ Den_q(mask_q) .* ...
% %     (VXq(mask_q).^2 - VXq(mask_q).*Ug.*cos(DXq(mask_q))) ...
% %     - Ilim^2 .* sqrt(Rad_q(mask_q)) ./ Den_q(mask_q) .* ...
% %     VXq(mask_q).*Ug.*sin(DXq(mask_q));
% 
% % 数值计算
% Den_q = VXq.^2 + Ug^2 - 2*VXq.*Ug.*cos(DXq);
% Rad_q = Den_q/Ilim^2 - Xg^2;
% Imag_q = sqrt(Den_q./(Xg^2 + Rg^2));
% 
% Q_val = nan(size(DXq));
% 
% % 只画 Imag >= Ilim 的区域
% mask_q = (Den_q > 1e-8) & (Rad_q >= 0) & (Imag_q >= Ilim);
% 
% Q_val(mask_q) = ...
%     Xg*Ilim^2 ./ Den_q(mask_q) .* ...
%     (VXq(mask_q).^2 - VXq(mask_q).*Ug.*cos(DXq(mask_q))) ...
%     - Ilim^2 .* sqrt(Rad_q(mask_q)) ./ Den_q(mask_q) .* ...
%     VXq(mask_q).*Ug.*sin(DXq(mask_q));
% 
% 
% % 新开图窗画三维图
% f3 = figure(3);
% surf(DXq, VXq, Q_val, 'EdgeColor', 'none');
% grid on;
% view(45,30);
% colorbar;
% shading interp;
% 
% xlabel('\delta_x');
% ylabel('voltage_x');
% zlabel('Q');
% title('Reactive power under current limit');
% 
% %% SEP/UEP search for VOC current-limit circle model
% 
% clear xep flag v V Lambda A sig m
% ep_set = [];
% m = 1;
% 
% torralence = 1e-2;
% options = optimoptions('fsolve', ...
%     'FunctionTolerance',1e-10, ...
%     'MaxIterations',100000, ...
%     'OptimalityTolerance',1e-10, ...
%     'Display','off');
% 
% for n = 1:length(x_set(1,:))
%     xep0 = x_set(:,n);
% 
%     % 求零点
%     [xep, ferr, exitflag] = fsolve(@f_VOC_post_cl_circle, xep0, options);
% 
%     if exitflag > 0 && maxabs(ferr) < torralence
%         if isnewxep(ep_set, xep, torralence)
% 
%             % 数值雅可比
%             A = numerical_jacobian(@f_VOC_normal_cl_circle, xep);
% 
%             % 特征值 / 特征向量
%             [V, Lambda] = eig(A);
%             Lambda = diag(Lambda);
% 
%             sig = sign(sign(real(Lambda)) + 0.1);   % zero counted as positive
%             sig = (sig + 1)/2;                      % [0,1]
%             flag = sum(sig);                        % number of non-negative eigenvalues
% 
%             v  = V(:, ~sig);                        % stable sub-space
%             vv = V(:, logical(sig));                % unstable/center sub-space
% 
%             ep_set(m).xep    = real(xep);
%             ep_set(m).A      = A;
%             ep_set(m).Lambda = Lambda;
%             ep_set(m).V      = V;
%             ep_set(m).v      = v;
%             ep_set(m).vv     = vv;
%             ep_set(m).flag   = flag;
% 
%             m = m + 1;
%         end
%     end
% end
% 
% %
% for mm = 1:length(ep_set)
%     disp_v('Index', mm);
%     disp_v('Equilibrium', ep_set(mm).xep);
%     disp_v('Eigenvalue', ep_set(mm).Lambda);
%     disp_v('Eigenvector', ep_set(mm).V);
%     disp('------------------------------');
% end
% 
% %
% clear ep_set_ext;
% for n = 1:length(ep_set)
%     ep_set_ext(n)=ep_set(n);
% end
% 
% % stable manifold of UEPs for VOC current-limit circle model
% figure(f1)
% opt_back = odeset('RelTol',1e-6,'Events',@stop_event_cl);
% 
% for mm = 1:length(ep_set_ext)
%     xep  = ep_set_ext(mm).xep;
%     flag = ep_set_ext(mm).flag;
% 
%     scatter(xep(1),xep(2),color_code{flag+1},'LineWidth',1.5);
% 
%     if flag == 1
%         v = real(ep_set_ext(mm).v);
%         v = v / norm(v);
% 
%         perturb = 1e-4;
% 
%         [~, x_p] = ode78(@f_backward_cl_circle,[0,1],xep + perturb*v,opt_back);
%         [~, x_n] = ode78(@f_backward_cl_circle,[0,1],xep - perturb*v,opt_back);
% 
%         x_all = [flip(x_n,1); x_p];
%         plot(x_all(:,1),x_all(:,2),'m-','LineWidth',1.5);
%     end
% end
% 
%         figure(f1)
%         options = odeset('RelTol',1e-10,'Events',@stop_event_cl);
%         [t_fault , x_all] = ode78(@f_fault_cl_circle,[t_start,t_start+t_c],[prefault_SEP(1);prefault_SEP(2)],options);
%         voltage_fault= x_all(:,2);
%         delta_fault =  x_all(:,1);
% 
%         [t_postfault , x_all2] = ode78(@f_post_cl_circle,...
%             [t_fault(end),t_end],...
%             x_all(end,1:2),...
%             options);
%         voltage_post= x_all2(:,2);
%         delta_post = x_all2(:,1);
% 
% 
%         plot(delta_fault,voltage_fault,'r-','linewidth',1.5);
%         plot(delta_post(1),voltage_post(1),'k.','MarkerSize',15);
%         plot(delta_fault(1),voltage_fault(1),'k.','MarkerSize',15);
%         plot(delta_post,voltage_post,'b-','linewidth',1.5)















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
    case "GFMQ"
        dfdt = f_GFMQ_normal(x);

    case "VOC"
        dfdt = f_VOC_normal(x);
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
        case "GFMQ"
          dfdt = f_GFMQ_fault(x);
        case "VOC"
          dfdt = f_VOC_fault(x);
    end
end
function dfdt = f_post(t,x)
global system;
    switch system 
        case "GFMQ"
          dfdt = f_GFMQ_normal(x);
        case "VOC"
          dfdt = f_VOC_normal(x);
    end
end

function dfdt = f_fault_cl_circle(t,x)
global system;
    switch system 
        case "VOC"
          dfdt = f_VOC_fault_cl_circle(x);
    end
end
function dfdt = f_post_cl_circle(t,x)
global system;
    switch system 
        case "VOC"
          dfdt = f_VOC_post_cl_circle(x);
    end
end
function dfdt = f_prefault(t,x)
global system;
    switch system 
        case "GFMQ"
          dfdt = f_GFMQ_normal(x);
        case "VOC"
          dfdt = f_VOC_normal(x);
    end
end

function [value, isterminal, direction] = stop_event(t,x)
    value = x(2) - 1e-2;   % 当 x(2) = 1e-2 时触发
    isterminal = 1;        % 终止积分
    direction = -1;        % 只在 x(2) 从大到小穿过 1e-2 时触发
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

function dfdt = f_backward_cl_circle(t,x)
    dfdt = -f_VOC_post_cl_circle(x);
end

function [value, isterminal, direction] = stop_event_cl(t,x)
    value = x(2) - 1e-2;
    isterminal = 1;
    direction = -1;
end