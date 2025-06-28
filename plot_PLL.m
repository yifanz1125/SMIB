%%
DeltaGFL=log_Delta.signals(1).values;
OmegaGFL=log_Omega.signals(1).values;


T_deta=Ts*10;
t_start=1.1;
t_end1 = t_start+t_c;
t_end2 = 1.5;


DeltaGFL_draw1=DeltaGFL(t_start/T_deta+1:t_end1/T_deta+1);
DeltaGFL_draw2=DeltaGFL(t_end1/T_deta+2:t_end2/T_deta+1);
OmegaGFL_draw1=OmegaGFL(t_start/T_deta+1:t_end1/T_deta+1);
OmegaGFL_draw2=OmegaGFL(t_end1/T_deta+2:t_end2/T_deta+1);


for n = 1:length(DeltaGFL_draw1)-1  % Check the "continuous" property of phase angle
    if (DeltaGFL_draw1(n)-DeltaGFL_draw1(n+1)) > 2*pi*4/5
        DeltaGFL_draw1(n+1:length(DeltaGFL_draw1)) = DeltaGFL_draw1(n+1:length(DeltaGFL_draw1)) + 2*pi;
    elseif (DeltaGFL_draw1(n)-DeltaGFL_draw1(n+1)) <  -2*pi*4/5
        DeltaGFL_draw1(n+1:length(DeltaGFL_draw1)) = DeltaGFL_draw1(n+1:length(DeltaGFL_draw1)) - 2*pi;
    end
end

for n = 1:length(DeltaGFL_draw2)-1  % Check the "continuous" property of phase angle
    if (DeltaGFL_draw2(n)-DeltaGFL_draw2(n+1)) > 2*pi*4/5
        DeltaGFL_draw2(n+1:length(DeltaGFL_draw2)) = DeltaGFL_draw2(n+1:length(DeltaGFL_draw2)) + 2*pi;
    elseif (DeltaGFL_draw2(n)-DeltaGFL_draw2(n+1)) <  -2*pi*4/5
        DeltaGFL_draw2(n+1:length(DeltaGFL_draw2)) = DeltaGFL_draw2(n+1:length(DeltaGFL_draw2)) - 2*pi;
    end
end

hold on
plot(DeltaGFL_draw1,OmegaGFL_draw1,'m','LineWidth',1.5);
hold on
plot(DeltaGFL_draw2,OmegaGFL_draw2,'g','LineWidth',1.5);




