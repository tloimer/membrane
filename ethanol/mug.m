function mu=mug(T)
%MUG(T)     Dynamic viscosity of the vapor [Pas].
%
%  Ethanol.
%  Valid for T/Tc < 1.5 (0.6?).
%  Prediction method from Perry (1997), p 2-363.

%  Perry (1997), Table 2-365: mu=8.35e-6 at 20C;
%  muold(293) = 9.9e-6, mu=7.4e-6 -> 10% deviation.

%muold = sqrt(T/400)*11.6e-6

%[R M]=molm;

Tc=513.92;
%pc=6.12e6;

Tr=T/Tc;

% 4.6e-4*1e-3*0.00034 = 1.564e-10
% Ethanol: M=46.07.
% 4.6e-4*1e-3*0.00034*sqrt(M)*(pc^4/Tc)^(1/6) = 1.2550e-5.
mu=1.255e-5*Tr.^0.94;
