function nu=nu(T,a)
%NU(T,A)    Effective kinematic viscosity of the 2ph-mixture, plug [m2/s].
%  NU = ( A*mug + (1-A)*mul ) / rho .
%
%  Calls MUG, MUL, RHO, V, PS.

nu = (a.*mug(T)+(1-a).*mul(T)) ./ ((1-a).*rho(T)+a./v(T,ps(T)));
