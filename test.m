        [t_postfault , x_all2] = ode78(@f_post,[0,0.0199],[2 0.25],odeset('RelTol',1e-10));
        voltage_post= x_all2(:,2);
        delta_post = x_all2(:,1);

        figure(f1)
    
        plot(delta_post(1),voltage_post(1),'k.','MarkerSize',15);
        plot(delta_post,voltage_post,'b-','linewidth',1.5)

function dfdt = f_post(t,x)
          dfdt = f_VOC_normal(x);
end