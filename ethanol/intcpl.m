function cpT=intcpl(T0,T1)
%INTCPL(T0,T1)  Integrate CPL dT from T0 to T1 [J/kg].
%
%  Calls MOLM.
%  See also CPL.
%
%  Ethanol.
%  Valid for 159 < T < 390.
%  From Perry (1997), Table 2-196.

[R M]=molm;

c1=1.0264e5;
c2=-139.63;
c3=-3.0341e-2;
c4=2.0386e-3;

dT=T1-T0;
dT2=(T1.^2-T0.^2)/2;
dT3=(T1.^3-T0.^3)/3;
dT4=(T1.^4-T0.^4)/4;

%cp = ( c1 + c2*T + c3*T.^2 + c4*T.^3 )/M;
cpT = ( c1*dT + c2*dT2 + c3*dT3 + c4*dT4 )/M;

%%%
% T = 298 to 348 K.
% Cp(liq) = 98.39 + 0.5368(T/K-273.25) J/mol*K (298 to 348 K).; DH

%[R M]=molm;
%dT=T1-T0;
%dT2=T1.^2-T0.^2;
%cpTold = (98.39e3*dT + 0.5368e3*(dT2/2-273.15*dT))/M
