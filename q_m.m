function q = q_m(T,a)
%Q_M(T,A)   Heat flux over mass flux, q/m, in the 2ph-region [J/kg].
%
%  Calls K, NU, DPSDT, KAPPA.

q = k(T,a).*nu(T,a)./(dpkdT(T)*kappa);
