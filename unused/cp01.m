function dcp = cp01(T,p0,p1)
%CP01       Difference of cp for different pressures [J/kgK].
%  CP01(T,P0,P1) gives the difference of cp between P0 and P1,
%  cp(p1)-cp(p0). CP01 calculates the difference based on the virial
%  equation truncated after the first term.
%  (Equivalent to DCPDP(T)*(P1-P0).)
%
%  Calls MOLM, TD2BDT.
%  See also V, DCPDP, DHDP.

% see dcpdp.nb

%[B TdBdT T2d2B] = T2d2b(T);
Td2B = Td2bdT(T);
[R M] = molm;

%TB=TdBdT./B;
%T2B=T2d2B./B;

%a = 250*R*T;
%b = a + p1.*B;
% 20*sqrt(10)=63.24555320
%cp1= -(sqrt(R.*b)/(63.24555320*M.*sqrt(T))) * ( (a.*(1-TB).^2/b)...
%  + 2*T2B + TB.*(2-TB) -1 );
%b = a + p0.*B;
%cp0= -(sqrt(R.*b)/(63.24555320*M.*sqrt(T))) * ( (a.*(1-TB).^2/b)...
%  + 2*T2B + TB.*(2-TB) -1 );
%dcpold = cp1 - cp0

dcp = -Td2B.*(p1-p0)./(1000*M);
