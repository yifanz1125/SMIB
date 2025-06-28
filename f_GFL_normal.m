function dfdt = f_GFL_normal(x)

    model = evalin('base','model');
    Id = evalin('base','Id');
    Iq = evalin('base','Iq');
    Ug =evalin('base','Ug');
        
    kp = evalin('base','kp');
    ki = evalin('base','ki');

    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Lg = Xg/Ws;
    Rg = evalin('base','Rg');

    Wmax = evalin('base','w_max');
    Wmin = evalin('base','w_min');

    W_g = evalin('base','W_g');

    delta = x(1);
    Int = x(2);

    



    if model == "original"
        Vq = (Xg*Id+Rg*Iq-Ug*sin(delta)+Id*Lg*Int)/(1-Id*Lg*kp);
%         omega_pll=kp*Vq+Int;
%         if omega_pll>=Wmax
%             Vq = (Xg*Id+Rg*Iq-Ug*sin(delta)+Id*Lg*Wmax);
%         elseif omega_pll<=Wmin
%             Vq = (Xg*Id+Rg*Iq-Ug*sin(delta)+Id*Lg*Wmin);
%         else
%             Vq = (Xg*Id+Rg*Iq-Ug*sin(delta)+Id*Lg*Int)/(1-Id*Lg*kp);
%         end

%         if Int>=Wmax && Vq>=0
%             dfdt(2)=0;
%             dfdt(1)=Wmax;
%         elseif Int<=Wmin && Vq<=0
%             dfdt(2)=0;
%             dfdt(1)=Wmin;
%         else
%             omega=kp*Vq+Int;
%             if omega>=Wmax
%                 dfdt(1) = Wmax - W_g; %ddelta
%                 dfdt(2) = ki*Vq; %dint
%             elseif omega<=Wmin
%                 dfdt(1) = Wmin - W_g; %ddelta
%                 dfdt(2) = ki*Vq; %dint
%             else
                dfdt(1) = kp*Vq+Int - W_g; %ddelta
                dfdt(2) = ki*Vq; %dint
%             end
%         end

    elseif model == "no_grid_dy" % ignore w2 & w1   
        Vq = Xg*Id+Rg*Iq-Ug*sin(delta);
        %omega = kp*Vq+Int;
        if Int>=Wmax && Vq>=0
            dfdt(2)=0;
            dfdt(1)=Wmax;
        elseif Int<=Wmin && Vq<=0
            dfdt(2)=0;
            dfdt(1)=Wmin;
        else
            omega=kp*Vq+Int;
            if omega>=Wmax
                dfdt(1) = Wmax; %ddelta
            elseif omega<=Wmin
                dfdt(1) = Wmin; %ddelta
            else
                dfdt(1) = kp*Vq+Int; %ddelta
            end
            dfdt(2) = ki*Vq; %dint
        end

    else 
        dfdt(1) = 0;
        dfdt(2) = 0;
    end
    dfdt = dfdt.';

    end