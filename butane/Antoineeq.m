function A = Antoineeq(T)
%ANTOINEEQ  Provides the coefficients in Antoine's equation for butane.
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
% Butane. n-butane. CAS 106-97-8.
% Range:  196 < T < 425.1
%
% Called by PS, TS.
% See also SUBSTANCE.

A = zeros(1,8);
% A B C T0 Tc n E F
A(1:3) = [5.93266+3 935.773 -34.361];
T0 = 288;
if T > T0
  A(4:8) = [T0 425.1 2.14767 -175.62 12204];
end
