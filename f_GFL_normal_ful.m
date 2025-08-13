function dfdt = f_GFL_normal_ful(x)

    Id = evalin('base','Id');
    Iq = evalin('base','Iq');
    Ug =evalin('base','Ug');
        
    kp = evalin('base','kp');
    ki = evalin('base','ki');

    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Lg = Xg/Ws;
    Rg = evalin('base','Rg');
    Lf = evalin('base','Lf');
    Xsum = evalin('base','Xsum');
    Lsum = evalin('base','Lsum');
    beta = evalin('base','beta');

    delta = x(1);
    Int = x(2);
    id = x(3);
    iq = x(4);
    yd = x(5);
    yq = x(6);

    eq = beta*Lf*(Iq-iq)+beta*Lf*yq;

    Vq=Lg/Lsum*eq+Lf/Lsum*Rg*iq-Lf/Lsum*Ug*sin(delta);

     omega = Vq*kp+Int;

     dfdt(1) = kp*Vq+Int; 
     dfdt(2) = ki*Vq; 

     dfdt(3) = beta*Lf/Lsum*(Id-id)+beta*Lf/Lsum*yd+Ws*iq+iq*omega-Rg/Lsum*id-1/Lsum*Ug*cos(delta); 
     dfdt(4) = beta*Lf/Lsum*(Iq-iq)+beta*Lf/Lsum*yq-Ws*id-id*omega-Rg/Lsum*iq+1/Lsum*Ug*sin(delta); 

     dfdt(5) = (Id-id)/4*beta; 
     dfdt(6) = (Iq-iq)/4*beta;  


    
    dfdt = dfdt.';

    end