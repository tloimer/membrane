function dcpdp = dcpdp(T)
%DCPDP(T)   Derivative dcp/dp for T constant [m3/kgK].
%  DCPDP(T) calculates the derivative dcp/dp, based on a virial
%  equation truncated after the first term.
%
%  Calls MOLM, TD2BDT.
%  See also V, VIRIAL.

%[B TdBdT T2d2B] = T2d2b(T);
Td2B = Td2bdT(T);
[R M] = molm;

%TB=TdBdT./B;
%T2B=T2d2B./B;

%a = 250*R*T;
%b = a + p.*B;
% 20*sqrt(10)=63.24555320
%dcpdpold = -(sqrt(R.*b)/(63.24555320*M.*sqrt(T))) * ( (a.*(1-TB).^2/b)...
%  + 2*T2B + TB.*(2-TB) -1 )

dcpdp = -Td2B./(1000*M);
