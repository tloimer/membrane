function [m,ms] = mT1eqT2(T1,p1,p2,s,ms,accuracy)
%MT1EQT2     Mass flux for adiabatic flow but q2 such that T1 equal to T2.
%  MT1EQT2(T1,P1,P2,SUBSTANCE,MS) returns the mass flux [kg/m2s] for
%  adiabatic flow of SUBSTANCE through a stack of layered membranes MS. The
%  membrane struct MS is constructed with MSTACKSTRUCT.
%
%  MT1EQT2(T1,P1,P2,SUBSTANCE,MS,'crude') sets crude solver tolerances.
%
%  [M,MS] = MT1EQT2(T1,P1,P2,SUBSTANCE,MS) writes the solution to MS.
%
%  Calls ASYM.
%
%  See also ASYM, MNUMADIABAT, MSTACKSTRUCT, SUBSTANCE.

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

if nargin < 6
  accuracy = 'accurate';
end

% Set the downstream state, downstreamstate(T2,p2,a2,q2,s,m)
state2 = downstreamstate(T1,p2,[],0,s); % thus, m is empty

% Copy some values to the membrane struct;
ms.T1 = T1;
ms.p1in = p1;
ms.T2 = T1;
ms.p2 = p2;
ms.a2 = state2.a;

% Calculate and set flowsetup-structure. To have substance properties being
% computed in a wide enough range of temperatures, calculate the downstream
% temperature for adiabatic flow.
T2 = s.intjt(T1,p1,p2);
% writeflowsetups() by itself adds a considerable range above T1
ms = ms.writeflowsetups(T1,T2,s,ms);

% Compute h1 - h2; For positive Joule-Thomson coefficient, h2 > h1.
% mu_JT = -(dh/dp)_T/cp, cp > 0, for mu_JT > 0 yields (dh/dp)_T < 0;
% With p2 < p1 results h2 > h1.
if p1 == p2
    h12 = 0;
else
    h12 = integral(@(p) s.dhdp(T1,p), p2, p1);
end

% Find an interval for m where the residual pressure, p1sol - p1, changes sign.
solver = solverstruct(accuracy);

% findinterval.m uses m = 0 as one point, hence rather err towards large mass
% flows - but then, the achieved temperatures might be too high, out of range
if isfield(ms,'mguess') && isscalar(ms.mguess) && ~isempty(ms.mguess)
  mguess = ms.mguess;
elseif isinf(s.ps(T1)) % above the critical temperature
  mguess = ms.mfluxknudsen(T1,p1,p2,s,ms) + ms.mfluxviscous(T1,p1,p2,s,ms);
else
  mguess = ms.mfluxliquid(T1,p1,p2,s,ms);
end

[minterval,pinterval] = findinterval(@presiduum,mguess,p2-p1);
m = findzero(@presiduum,[minterval; pinterval],(p1-p2)/1000);

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

end %%% END MT1EQT2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MT1EQT2 %%%
