function Td2B = Td2bdT(T)
%TD2BDT(T)  Derivative of virial coefficient T*d^2b/dT^2 [cm3/molK].
%
%  Calls VIRIAL.
%  See also B, TDBDT, T2D2B, V.

% R=8314.4;
[b0 b1 b2 b3]=virial;

%B = b0 + b1./T + b2./T.^2 + b3./T.^3;
% dBdT= -b1/T^2 - 2*b2/T^3 - 3*b3/T^4;
% d2B = 2b1/T^3 + 6 b2/T^4 + 12b3/T^5;
%TdBdT= -b1./T - 2*b2./T.^2 - 3*b3./T.^3;
Td2B = 2*b1./T.^2 + 6*b2./T.^3 + 12*b3./T.^4;
