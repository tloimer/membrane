function cp=cpl(T)
%CPL(T)     Constant pressure heat capacity of the liquid [J/kgK].
%
%  Calls MOLM.
%
%  Ethanol.
%  Valid for 159K < T < 390K.
%  From Perry (1997), Table 2-196.

[R M]=molm;

c1=1.0264e5;
c2=-139.63;
c3=-3.0341e-2;
c4=2.0386e-3;

cp = ( c1 + c2*T + c3*T.^2 + c4*T.^3 )/M;

%%%
%  Data taken from NIST: Pedersen, Kay, et al., 1975
%
% T = 298 to 348 K.
% Cp(liq) = 98.39 + 0.5368(T/K-273.25) J/mol*K (298 to 348 K).; DH

%cpold = (98.39e3 + 0.5368e3*(T-273.15))/M
