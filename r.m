function r = r(T)
%R(T)       Evaporation enthalpy [J/kg].
%  Calculated by Clausius-Clapeyron equation, dp/dT = r/(T(v''- v')).
%
%  Calls V, PS, RHO, DPSDT.

% NIST: Majer and Svoboda, 1985.  T=298-469K.
% r1 = 50.43e6*exp(+0.4475*(T/513.9)).*((1-T/513.9).^0.4989)/46.07;
% Tabulated values also in Perry (1997), Table 2-256, p. 2-235.

r = T.*(v(T,ps(T))-1./rho(T)).*dpsdT(T);
