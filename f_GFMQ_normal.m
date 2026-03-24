function dfdt = f_GFMQ_normal(x)

    Ug =evalin('base','Ug');
    Pm =evalin('base','Pm');

    m_gfm = evalin('base','m_gfm');
        
    k_q = evalin('base','k_q');

    tau_q = evalin('base','tau_q');

    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Lg = Xg/Ws;
    Rg = evalin('base','Rg');

    Vref = evalin('base','Vgfm');
    Qref = evalin('base','Qref');


    delta = x(1);
    Vgfm = x(2);
    
    P = Rg*(Vgfm^2-Vgfm*Ug*cos(delta))/(Rg^2+Xg^2)+Xg*Vgfm*Ug*sin(delta)/(Rg^2+Xg^2);
    Q = Xg*(Vgfm^2 - Vgfm*Ug*cos(delta))/(Rg^2+Xg^2) - Rg*Vgfm*Ug*sin(delta)/(Rg^2+Xg^2); 
    

  
    dfdt(1) = m_gfm(Pm - P)*Ws;  %delta
    dfdt(2) = k_q/tau_q*(Qref-Q)+1/tau_q*(Vref-Vgfm);  %voltage
 
    dfdt = dfdt.';

    end