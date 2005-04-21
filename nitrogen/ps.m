function ps = ps(T)
%PS(T)      Vapor pressure [Pa].
%
%  See also TS, DPSDT.
%
%  Nitrogen.
%  Antoine equation.
%  Range 63.14K < T < 126K accuracy not given
%  From http://webbook.nist.gov.

Aa= 8.7362; % 5 added, for Pa, not bar
Ab=264.651;
Ac=-6.788;

ps = 10.^(Aa-Ab./(Ac+T));
