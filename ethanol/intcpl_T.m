function cp = intcpl_T(T0,T1)
%INTCPL_T(T0,T1) Integrate CPL/T dT from T0 to T1 [J/kgK].
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
%dT4=(T1.^4-T0.^4)/4;

%cp = ( c1 + c2*T + c3*T.^2 + c4*T.^3 )/M;
cp = ( c1*log(T1./T0) + c2*dT + c3*dT2 + c4*dT3 )/M;
