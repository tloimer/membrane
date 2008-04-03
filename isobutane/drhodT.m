function dr = drhodT(T)
%DRHODT(T)  Derivative drho/dT of the density of the liquid [kg/m3K].
%
%  See also RHO.
%
%  Ethanol.
%  Valid for 191K < T < 400K.
%  From Landolt-Börnstein: Group IV, vol. 8G (2000).
%  Perry (1997) gives an equation with slightly larger range of validity.

%t0= 1.16239e3;
t1=-2.25788;
t2=5.30621e-3;
t3=-6.63070e-6;

%r = t0 + t1.*T + t2.*T.^2 + t3.*T.^3;
dr =      t1    + 2*t2.*T  + 3*t3.*T.^2;
