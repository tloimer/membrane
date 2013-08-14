function [m,ms] = mnumadiabat(T1,p1,p2,s,ms,accuracy)
%MNUMADIABAT Mass flux under adiabatic conditions.
%  MNUMADIABAT(T1,P1,P2,SUBSTANCE,MS) returns the mass flux. The membrane struct
%  MS is constructed with MSTACKSTRUCT.
%
%  [M,MS] = MNUMADIABAT(T1,P1,P2,SUBSTANCE,MS) writes the solution to MS.
%
%  Calls ASYM.
%
%  See also ASYM, MSTACKSTRUCT, SUBSTANCE.

% Check the upstream state
if s.ps(T1) < p1
  error('The upstream state is a liquid. This is not implemented.');
end
if nargin < 6
  accuracy = 'accurate';
end

% Calculate the downstream temperature
T2 = s.intjt(T1,p1,p2);

% Set the downstream state, downstreamstate(T2,p2,a2,q2,s,m)
state2 = downstreamstate(T2,p2,[],0,s); % thus, m is empty

% Copy some values to the membrane struct;
ms.T1 = T1;
ms.p1in = p1;
ms.T2 = T2;
ms.p2 = p2;
ms.a2 = state2.a;
ms.q2 = state2.q; % = 0

% Calculate and set flowsetup-structures
ms = ms.writeflowsetups(T1,T2,s,ms);

% Find an interval for m where the residual pressure, p1sol - p1, changes sign.
solver = solverstruct('crude');
presiduum = @(m) asym(m,state2,ms,solver) - p1;

% findinterval.m uses m = 0 as one point, hence rather err towards large mass
% flows - but then, the achieved temperatures might be too high, out of range
if isfield(ms,'mguess') && isscalar(ms.mguess) && ~isempty(ms.mguess)
  mguess = ms.mguess;
else
  mguess = ms.mfluxliquid(T1,p1,p2,s,ms);
end
[minterval,pinterval] = findinterval(presiduum,mguess,p2-p1);

% previously, used crude solver for findinterval(), then accurate solver
solver = solverstruct('accurate');
%presiduum = @(m) asym(m,state2,ms,solver) - p1;

%m = fzero(presiduum,minterval);

% findzero needs a different call to presiduum
presiduum = @(m,solver) asym(m,state2,ms,solver) - p1;
m = findzero(presiduum,[minterval; pinterval],(p1-p2)/1000,solver);

% Now write the solution
ms.m = m;
solver.writesolution = true;
solver.fullsolution = true;
[p1,ms] = asym(m,state2,ms,solver);
ms.p1sol = p1;
