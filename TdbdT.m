function [B,TdBdT] = TdbdT(T)
%TDBDT(T)   Virial coefficient B and derivative T*dB/dT [cm3/mol].
%  [B TDBDT] = TDBDT(T).
%
%  Calls VIRIAL.
%  See also V, B.

% R=8314.4;
[b0 b1 b2 b3]=virial;

B = b0 + b1./T + b2./T.^2 + b3./T.^3;
% dBdT= -b1./T.^2 - 2*b2./T.^3 - 3*b3./T.^4;
TdBdT= -b1./T - 2*b2./T.^2 - 3*b3./T.^3;
