function [r9,ms] = rnumadiabat(m,T1,p1,p2,s,ms,accuracy)
%RNUMADIABAT Mass flux for adiabatic flow through layered membranes.
%  RNUMADIABAT(T1,M,P1,P2,SUBSTANCE,MS) returns the radius R9 [m] within the
%  membrane for adiabatic flow of SUBSTANCE through a stack of layered membranes
%  MS. The membrane struct MS is constructed with MSTACKSTRUCT.
%
%  RNUMADIABAT(T1,P1,P2,SUBSTANCE,MS,'crude') sets crude solver tolerances.
%
%  [M,MS] = RNUMADIABAT(T1,P1,P2,SUBSTANCE,MS) writes the solution to MS.
%
%  Calls ASYM.
%
%  See also ASYM, MSTACKSTRUCT, MTCONST, SUBSTANCE.

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

% Calculate the downstream temperature
T2 = s.intjt(T1,p1,p2);

% Set the downstream state, downstreamstate(T2,p2,a2,q2,s,m)
state2 = downstreamstate(T2,p2,[],0,s); % thus, m is empty TODO - why 

% Copy some values to the membrane struct;
ms.m = m;
ms.T1 = T1;
ms.p1in = p1;
ms.T2 = T2;
ms.p2 = p2;
ms.a2 = state2.a;
ms.q2 = state2.q; % = 0

solver = solverstruct(accuracy);
function pres = presiduum(curv)
	ms = ms.writecurvsetups(curv, T1, T2, s, ms);
	pres = asym(m, state2, ms, solver) - p1;
end

curv = fzero(@presiduum, 4./ms.membrane(1).layer(1).matrix.dia);

% Now write the solution
solver.writesolution = true;
solver.fullsolution = true;
[p1,ms] = asym(m,state2,ms,solver);

r9 = 2./curv;	% for cylindrical pores
end
