function [t,q] = temp(z,m,T0,q0,z0,cp,k)
%TEMP       Solve the energy equation.
%  T = TEMP(Z,M,T0,Q0,Z0,CP,K) solves the energy equation for constant
%  material properties CP and K, constant mass flux M and initial
%  conditions T0 and Q0 at position Z0. Returns the temperature T at
%  position Z.
%
%  [T,Q] = TEMP(Z,M,T0,Q0,Z0,CP,K) solves as above and returns the
%  temperature T and heat flux Q at position Z.
%
%  T = T0 + (Q0/(M*CP))[ 1 - exp((M*CP/K)(Z-Z0)) ],
%  Q = Q0 * exp((M*CP/K)(Z-Z0)).
%
%  Called from FLOWBACK.

t = T0 + (q0./(m.*cp)).*( 1 - exp((m.*cp./k).*(z-z0)) );
if (nargout==2)
  q = q0*exp((m.*cp./k).*(z-z0));
end
