function kc = kappac(T)
%KAPPAC(T)  Critical permeability [m2].
%
%  Calls NU, K, DPSDT, R.

kc = nu(T,0).*k(T,0)./(dpsdT(T).*r(T));
