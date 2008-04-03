function A = Antoineeq(T)
%ANTOINEEQ  Provides the coefficients in Antoine's equation for isobutane.
%
% A = ANTOINEEQ(T) returns the coefficients for the classical and the extended
% Antoine's equation in the 1-by-8 Matrix A, A = [a b c T0 Tc n E F]. If T <= T0
% then A(4) = 0, only a, b, c are set and the saturation pressure [Pa] is given
% by
%   p = 10^(a-b/(c+T)).
%
% If T>T0 then
%
%   p = 10^(a - b/(c+T) + 0.4329*chi^n + E*chi^8 + F*chi^12).
%
% See Landolt-Börnstein, New Series, Group IV: Physical Chemistry.
% Vapor Pressure of Chemicals, vol.  20A: J. Dykyj, J. Svoboda, R.C. Wilhoit,
% M.  Frenkel, K.R. Hall (1999).
%
% Isobutane. 2-Methylpropane. CAS 75-28-5.
% Range:  188 < T < 407.1
%
% Called by PS, TS.
% See also SUBSTANCE.

A = zeros(1,8);
% A B C T0 Tc n E F
A(1:3) = [6.00272+3 947.54 -24.28];
T0 = 268;
if T > T0
  A(4:8) = [T0 407.1 2.6705 -19.64 2792];
end
