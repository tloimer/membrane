function dps = dpsdT(T)
%DPSDT(T)   Derivative of the vapor pressure [Pa/K].
%
%  Calls PS.
%  See PS, TS.
%
%  Ethanol.
%  From an Antoine equation.
%  From Landolt-Börnstein: Group IV, vol. 20A (2000).
%  An equation with larger range of validity is given by Perry (1997).

%Aa= 6.923365 + 3; % convert to Pa, not kPa
Aa= 9.923365;
Ab=1410.46;
Ac=-64.636;

ps = 10.^(Aa-Ab./(Ac+T));
dps = ps*Ab*2.302585092994./(Ac+T).^2;
