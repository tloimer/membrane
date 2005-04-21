function dps = dpsdT(T)
%DPSDT(T)   Derivative of the vapor pressure [Pa/K].
%
%  See PS, TS.
%
%  Nitrogen.
%  Antoine equation.
%  Range 63.14K < T < 126K accuracy not given
%  From http://webbook.nist.gov.

Aa= 8.7362; % 5 added, for Pa, not bar
Ab=264.651;
Ac=-6.788;

ps = 10.^(Aa-Ab./(Ac+T));
dps = ps*Ab*2.302585092994./(Ac+T).^2;
% ln(10) = 2.30258..
