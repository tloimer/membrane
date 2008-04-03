function ps = ps(T)
%PS(T)      Vapor pressure [Pa].
%
%  See also TS, DPSDT.
%
%  Ethanol.
%  Antoine equation.
%  Valid for 341 < T < 358, error < 0.1% for 341.2K < T < 358K.
%  From Landolt-Börnstein: Group IV, vol. 20A (2000).
%  An equation with larger range of validity is given by Perry (1997),
%  see PSNEW(T).

Aa= 6.92365 + 3; % convert to Pa, not kPa
Ab=1410.46;
Ac=-64.636;

ps = 10.^(Aa-Ab./(Ac+T));
