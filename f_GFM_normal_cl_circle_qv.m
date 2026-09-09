function dfdt = f_GFM_normal_cl_circle_qv(x)
    Ug   = evalin('base','Ug');
    Pm   = evalin('base','Pm');
    Ws   = evalin('base','Ws');
    Xg   = evalin('base','Xg');
    Rg   = evalin('base','Rg');
    UN   = evalin('base','Vgfm');
    D    = evalin('base','D');
    J    = evalin('base','J');
    Ilim = evalin('base','Ilim');
    kq   = evalin('base','kq');
    Qref = evalin('base','Qref1');

    delta = x(1);
    omega = x(2);

    [P,~,~,~,~,~] = circle_qv_power(delta,Ug,Rg,Xg,UN,Ilim,kq,Qref);

    dfdt = [omega*Ws;
            (Pm-P)/J-D/J*omega];
end
