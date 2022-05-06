function [m,ms] = mnumTdiff(fac,T1,p1,p2,s,ms,accuracy)
%MNUMTDIFF  Mass flux for flow with given temperature difference.
%  MNUMTDIFF(FAC,T1,P1,P2,SUBSTANCE,MS) returns the mass flux [kg/m2s] for
%  adiabatic flow of SUBSTANCE through a stack of layered membranes MS. The
%  membrane struct MS is constructed with MSTACKSTRUCT. The temperature at the
%  downstream side of the struct is the temperature expected for a Joule-Thomson
%  process, T2AD, plus FAC times the temperature difference T1 - T2AD. Hence,
%  T1 - T2 = (1 - FAC) * (T1 - T2AD).
%  MNUMTDIFF(0,...) is equivalent to MNUMADIABAT(...).
%  MNUMTDIFF(1,...) is equivalent to MT1EQT2(...).
%
%  The boundary condition
%
%    T2 = T2AD + FAC * (T1 - T2AD)
%
%  when using the approximation dh = cp dT, because of
%
%    m * h1 + q1 = m * h2 + q2,  q1 = 0,  h2 = h1 + fac * cp * (T1 - T2ad)
%
%  is equivalent to
%
%   q2 = m * (h1 - h2) = -m * cp * fac * (T1 - T2ad)
%
%  where m is the mass flux denstiy, m = rho v.
%
%  MNUMTDIFF(FAC,T1,P1,P2,SUBSTANCE,MS,'crude') sets crude solver tolerances.
%
%  [M,MS] = MNUMTDIFF(FAC,T1,P1,P2,SUBSTANCE,MS) writes the solution to MS.
%
%  Calls ASYM.
%
%  See also ASYM, MNUMADIABAT, MNUMALPHA, MNUMHEATFLUX, MNUMRADIATION,
%  MNUMHTRAD, MNUMISO, MSTACKSTRUCT, MT1EQT2, SUBSTANCE.

% Some input sanitizing.
if s.ps(T1) < p1
    error([upper(mfilename)...
	': The upstream state is a liquid. This is not implemented.']);
elseif p2 < 0.
    error([upper(mfilename)...
	': The downstream pressure is negative. That is not possible.']);
elseif T1 < 0.
    error([upper(mfilename)...
	': The upstream temperature is below absolute zero. Impossible.']);
end

if nargin < 7
    accuracy = 'accurate';
end

% Downstream temperature for adiabatic downstream boundary condition
T2ad = s.intjt(T1,p1,p2);
if fac == 0
    T2 = T2ad;
elseif fac == 1
    T2 = T1;
else
    T2 = T2ad + fac * (T1 - T2ad);
end

% Enthalpy difference h12 = h1 - h2,
if fac == 0
    h12 = 0;
else
    h12 = -integral(@(T) s.cpg(T,p2), T2ad, T2);
end

% Set the downstream state, downstreamstate(T2,p2,a2,q2,s,m)
state2 = downstreamstate(T2,p2,[],0,s); % thus, m is empty

% Copy some values to the membrane struct;
ms.T1 = T1;
ms.p1in = p1;
ms.T2 = T2;
ms.p2 = p2;
ms.a2 = state2.a;

% writeflowsetups() by itself adds a considerable range above T1
ms = ms.writeflowsetups(T1,T2ad,s,ms);

% initial guess
if isfield(ms,'mguess') && isscalar(ms.mguess) && ~isempty(ms.mguess)
    mguess = ms.mguess;
elseif isinf(s.ps(T1)) % above the critical temperature
    mguess = ms.mfluxknudsen(T1,p1,p2,s,ms) + ms.mfluxviscous(T1,p1,p2,s,ms);
else
    mguess = ms.mfluxliquid(T1,p1,p2,s,ms);
end

solver = solverstruct(accuracy);
% m = findzero(@presiduum, mguess, (p1-p2)/1000);
[minterval,pinterval] = findinterval(@presiduum, mguess, p2-p1);
m = findzero(@presiduum, [minterval; pinterval], (p1-p2)/1000);

% Now write the solution
ms.m = m;
ms.q2 = m * h12;
solver.writesolution = true;
solver.fullsolution = true;
[p1,ms] = asym(m,state2,ms,solver);
ms.p1sol = p1;


%--- nested functions ------------------------------------- nested functions ---
function pres = presiduum(m) %---------------------------------------- presiduum
	% m h1 + q1 = m h2 + q2, q1 = 0  ->  q2 = m (h1 - h2)
	state2.q = m * h12;
	pres = asym(m, state2, ms, solver) - p1;
end %------------------------------------------------------------- end presiduum

end %%% END MNUMTDIFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MNUMTDIFF %%%
