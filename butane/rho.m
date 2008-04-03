function r = rho(T)
%RHO(T)     Density of the liquid [kg/m3].
%
%  Butane.  CAS 106-97-8.
%  Range: 135 K < T < 340K.
%
%  From Landolt-Börnstein: Group IV, vol. 8B (1996).  R. Wilhoit, K. Marsh, X.
%  Hong, N. Gadalla and M. Frenkel. An expression with larger range of validity
%  is also available.

t0= 892.907;
t1=-1.45679;
t2=2.87931e-3;
t3=-5.35281e-6;

r = t0 + t1.*T + t2.*T.^2 + t3.*T.^3;
