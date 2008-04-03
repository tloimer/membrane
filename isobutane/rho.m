function r = rho(T)
%RHO(T)     Density of the liquid [kg/m3].
%
%  Isobutane. 2-Methylpropane. CAS 75-28-5.
%  Range: 188 K < T < 326K.
%
%  From Landolt-Börnstein: Group IV, vol. 8B (1996).  R. Wilhoit, K. Marsh, X.
%  Hong, N. Gadalla and M. Frenkel. An expression with larger range of validity
%  is also available.

t0= 870.93;
t1=-1.36494;
t2=2.56419e-3;
t3=-5.32743e-6;

r = t0 + t1.*T + t2.*T.^2 + t3.*T.^3;
