function [m,ms] = mnumalpha(alpha,T1,p1,p2,s,ms,accuracy)
%MNUMALPHA   Mass flux for flow with downstreams heat transfer coefficient.
%  MNUMALPHA(ALPHA,T1,P1,P2,SUBSTANCE,MS) returns the mass flux [kg/m2s] for
%  adiabatic flow of SUBSTANCE through a stack of layered membranes MS. At the
%  downstream side of the struct, a heat transfer coefficent applies,
%
%    q2 = -ALPHA * (T1 - T2),
%
%  where q2 is the heat flux density at the downstream location in the
%  downstream direction, T1 is the upstream and T2 the downstream temperature.
%  The membrane struct MS is constructed with MSTACKSTRUCT.
%
%  MNUMALPHA(ALPHA,T1,P1,P2,SUBSTANCE,MS,'crude') sets crude solver tolerances.
%
%  [M,MS] = MNUMALPHA(ALPHA,T1,P1,P2,SUBSTANCE,MS) writes the solution to MS.
%
%  Calls ASYM.
%
%  See also ASYM, MNUMADIABAT, MT1EQT2, MNUMISO, MSTACKSTRUCT, SUBSTANCE.

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

% Set the downstream state, downstreamstate(T2,p2,a2,q2,s,m)
state2 = downstreamstate(T1,p2,[],0,s); % thus, m is empty

% Copy some values to the membrane struct;
ms.T1 = T1;
ms.p1in = p1;
ms.p2 = p2;
ms.a2 = state2.a;

% Calculate and set flowsetup-structure. To have substance properties being
% computed in a wide enough range of temperatures, calculate the downstream
% temperature for adiabatic flow.
T2ad = s.intjt(T1,p1,p2);
% writeflowsetups() by itself adds a considerable range above T1
ms = ms.writeflowsetups(T1,T2ad,s,ms);

% From the global energy balance,
%
%   m h1 + q1 = m h2 + q2
%
% state 2 and the heat flux q2 are not independent,
%
%   h2 = f(m,h1,q2)   if q1 = 0.
%
% Therefore, given q2, or alpha, the downstream state must be determined,
% depending on the mass flux density m.
%
% Construct a function h12(T) that gives the enthalpy difference
% h(T1,p1) - h(T,p2).
%
% Consider a T-s diagram,
%     T
%      ^                 p = const
%      | h(T1,p1)      .
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

% Get a minimum temperature,
Tmin = T2ad - max((T1 - T2ad)/10., 0.1);

if p1 == p2
    h1_T1p2 = 0;
else
    [~, h1_T1p2] = ode45(@(p,h) s.dhdp(T1,p), [p2 p1], 0);
    h1_T1p2 = h1_T1p2(end);
end
hT1p2_h2 = ode45(@(T,h) -s.cpg(T,p2), [T1 Tmin], 0);

h12 = @(T) h1_T1p2 + deval(hT1p2_h2, T);


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
ms.T2 = state2.T;
ms.q2 = state2.q;
solver.writesolution = true;
solver.fullsolution = true;
[p1,ms] = asym(m,state2,ms,solver);
ms.p1sol = p1;

%--- nested functions ------------------------------------- nested functions ---

function pres = presiduum(m) %---------------------------------------- presiduum
	% m h1 + q1 = m h2 + q2, q1 = 0  ->  q2 = m (h1 - h2)
	T2 = fzero(@(T) alpha*(T - T1) - m*h12(T), [T1 Tmin]);
	state2.T = T2;
	state2.q = m * h12(T2);
	pres = asym(m, state2, ms, solver) - p1;
end %------------------------------------------------------------- end presiduum

end %%% END MNUMALPHA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MNUMALPHA %%%
