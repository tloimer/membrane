function [r9,ms] = rnumalpha(alpha,m,T1,p1,p2,s,ms,accuracy)
%RNUMALPHA  Radius of curvature for flow through layered membranes.
%  RNUMADIABAT(ALPHA,M,T1,P1,P2,SUBSTANCE,MS) returns the radius R9 [m] within
%  the membrane for the flow with heat transfer coefficient ALPHA [W/m²K] of
%  SUBSTANCE through a stack of layered membranes MS.
%  The membrane struct MS is constructed with MSTACKSTRUCT.
%
%  RNUMALPHA(ALPHA,M,T1,P1,P2,SUBSTANCE,MS,'crude') sets crude tolerances.
%
%  [R9,MS] = RNUMALPHA(ALPHA,M,T1,P1,P2,SUBSTANCE,MS) writes the solution
%  to MS.
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

if nargin < 8
  accuracy = 'accurate';
end

% Calculate and set flowsetup-structure. To have substance properties being
% computed in a wide enough range of temperatures, calculate the downstream
% temperature for adiabatic flow.
% writeflowsetups() by itself adds a considerable range above T1

% Compute downstream temperature and heat flux.
%
% Consider a T-s diagram,
%        h(T1,p1)
%                        p = const
%      ^               .
%  T1 -|- +----------+ h(T1, p2)
%      |   `       .
%      |     `   .   anywhere between: h(T, p2)
%  T2 -|-      `
%      |        isenthalpic line
%      +--------------> s
%
%
% Use the following symbols,
%
%   h1_T1p2  = h(T1, p1) - h(T1, p2) = int_p2^p1 dhdp(T1,p') dp' ,
%
%   hT1p2_h2 = h(T1, p2) - h(T, p2) = -int_T1^T cp(T',p2) dT' ,
%
%   h12 = h1 - h2 = h1_T1p2 + hT1p2_h2 .
%
% Note, h1_T1p2 < 0, hT1p2_h2 > 0. Also, dhdp < 0, cp > 0.
%

switch alpha
case Inf
	% diabatic flow
	T2 = T1;
	if p1 == p2
		h1_T1p2 = 0;
	else
		[~, h1_T1p2] = ode45(@(p,h) s.dhdp(T1,p), [p2 p1], 0);
		h1_T1p2 = h1_T1p2(end);
	end
	q2 = m*h1_T1p2;
case 0
	% adiabatic flow
	T2 = s.intjt(T1,p1,p2);
	q2 = 0;
otherwise
	T2ad = s.intjt(T1,p1,p2);

	% Get a minimum temperature for the computation domain
	Tmin = T2ad - max((T1 - T2ad)/10., 0.1);

	if p1 == p2
		h1_T1p2 = 0;
	else
		[~, h1_T1p2] = ode45(@(p,h) s.dhdp(T1,p), [p2 p1], 0);
		h1_T1p2 = h1_T1p2(end);
	end
	hT1p2_h2 = ode45(@(T,h) -s.cpg(T,p2), [T1 Tmin], 0);

	h12 = @(T) h1_T1p2 + deval(hT1p2_h2, T);

	%  m h1 + q1 = m h2 + q2,
	%         q1 = 0   ->     q2 = m (h1 - h2),
	%  q2 = alpha(T2 - T1)  ->  alpha(T2 - T1) = m(h1 - h2)
	T2 = fzero(@(T) alpha*(T - T1) - m*h12(T), [T1 Tmin]);
	q2 = alpha*(T2 - T1);
end

% Set the downstream state, downstreamstate(T2,p2,a2,q2,s,m)
state2 = downstreamstate(T2,p2,[],q2,s,m); % thus, m is empty TODO - why 

% Copy some values to the membrane struct;
ms.m = m;
ms.T1 = T1;
ms.p1in = p1;
ms.T2 = T2;
ms.p2 = p2;
ms.a2 = state2.a;
ms.q2 = state2.q;

solver = solverstruct(accuracy);
function pres = presiduum(curv)
	ms = ms.writecurvsetups(curv, T1, T2, s, ms);
	pres = asym(m, state2, ms, solver) - p1;
end

curv = fzero(@presiduum, 4./ms.membrane(1).layer(1).matrix.dia);

% Now write the solution
solver.writesolution = true;
solver.fullsolution = true;
[p1sol,ms] = asym(m,state2,ms,solver);

% Do not use [ms.p1sol,ms] = asym(...) above, ms would overwrite ms.p1sol.
ms.p1sol = p1sol;
r9 = 2./curv;	% for cylindrical pores
end
