function Ts = Ts(p)
%TS(P)      Saturation temperature [K].
%
%  See PS. See also DPSDT.
%
%  Ethanol.
%  Antoine equation.
%  Range 341K < T < 514K, error < 0.1% for 341.2K < T < 358K.
%  From Landolt-Börnstein: Group IV, vol. 20A (2000).
%  An equation with larger range of validity is given by Perry (1997).

%log(10)=2.302585093
Aa= 6.92365 + 3; % convert to Pa, not kPa
Ab=1410.46;
Ac=-64.636;

%ps = 10.^(Aa-Ab./(Ac+T));
Ts = Ab./(Aa-log(p)/2.302585093)-Ac;
