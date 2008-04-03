function A = Antoineeq(T)
%ANTOINEEQ  Provides the coefficients in Antoine's equation for ethanol.
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
% Vapor Pressure of Chemicals, vol.  20B: J. Dykyj, J. Svoboda, R.C. Wilhoit,
% M.  Frenkel, K.R. Hall (2000).
%
% Ethanol. CAS 64-17-5.
% Range:  269 < T < 514
%
% Called by PS, TS.
% See also SUBSTANCE.

A = zeros(1,8);
% A B C T0 Tc n E F
% First Antoine eq., 269 < T < 341.2
if T <= 341.2
  A(1:3) = [7.33675+3 1648.22 -42.232];
  return
end
% Second Antoine eq., 341.2 < T < 358
T0 = 358;
A(1:3) = [6.92365+3 1410.46 -64.636];
if T > T0
  A(4:8) = [T0 513.9 0.434294 -255.71 300056];
end
