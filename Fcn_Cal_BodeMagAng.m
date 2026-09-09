function [Mag,Ang]=Fcn_Cal_BodeMagAng(Input)
    Mag=20*log10(abs(Input));
    Ang=angle(Input)*180/pi;
end