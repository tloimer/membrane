function Ts = Ts(p)
%TS(P)      Saturation temperature [K].
%
%  See also PS.
%
%  Ethanol.
%  Antoine equation. Valid for 269K < T < 341K.
%  From Landolt-Börnstein: Group IV, vol. 20A (2000).
%  An equation with larger range of validity is given by Perry (1997).

%log(10)=2.302585093
Aa= 6.923365 + 3; % convert to Pa, not kPa
Ab=1410.46;
Ac=-64.636;

%ps = 10.^(Aa-Ab./(Ac+T));
Ts = Ab./(Aa-log(p)/2.302585093)-Ac;
