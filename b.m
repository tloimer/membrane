function B = b(T)
%B(T)       Second virial coefficient [cm3/mol].
%
%  Calls VIRIAL.
%  Called from V.
%  See also V.

[b0 b1 b2 b3]=virial;

B = b0 + b1./T + b2./T.^2 + b3./T.^3;
