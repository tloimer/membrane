function [m,ms] = mnumiso(T1,p1,p2,s,ms,accuracy)
%MNUMISO     Mass flux for isothermal flow through layered membranes.
%  MNUMISO(T1,P1,P2,SUBSTANCE,MS) returns the mass flux [kg/m2s] for
%  isothermal flow of SUBSTANCE through a stack of layered membranes MS. The
%  membrane struct MS is constructed with MSTACKSTRUCT.
%
%  MNUMISO(T1,P1,P2,SUBSTANCE,MS,'crude') sets crude solver tolerances.
%
%  [M,MS] = MNUMISO(T1,P1,P2,SUBSTANCE,MS) writes the solution to MS.
%
%  Calls ISO.
%
%  See also ISO, MNUMADIABAT, MSTACKSTRUCT, MT1EQT2, SUBSTANCE.

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
% TODO: Does this work?
T2 = T1;
% writeflowsetups() by itself adds a considerable range above T1
ms = ms.writeflowsetups(T1,T2,s,ms);

% Find an interval for m where the residual pressure, p1sol - p1, changes sign.
solver = solverstruct(accuracy);
presiduum = @(m) iso(m, state2, ms, solver) - p1;

% findinterval.m uses m = 0 as one point, hence rather err towards large mass
% flows - but then, the achieved temperatures might be too high, out of range
if isfield(ms,'mguess') && isscalar(ms.mguess) && ~isempty(ms.mguess)
  mguess = ms.mguess;
elseif isinf(s.ps(T1)) % above the critical temperature
  mguess = ms.mfluxknudsen(T1,p1,p2,s,ms) + ms.mfluxviscous(T1,p1,p2,s,ms);
else
  mguess = ms.mfluxliquid(T1,p1,p2,s,ms);
end

[minterval,pinterval] = findinterval(presiduum,mguess,p2-p1);
m = findzero(presiduum,[minterval; pinterval],(p1-p2)/1000);

% Now write the solution
ms.m = m;
solver.writesolution = true;
solver.fullsolution = true;
[p1,ms] = iso(m,state2,ms,solver);
ms.p1sol = p1;
