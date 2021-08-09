function [m,ms] = mnumradiation(T1,p1,p2,s,ms,accuracy)
%MNUMRADIATION  Mass flux for flow with downstreams radiation heat transfer.
%  MNUMRADIATION(T1,P1,P2,SUBSTANCE,MS) returns the mass flux [kg/m2s] for
%  adiabatic flow of SUBSTANCE through a stack of layered membranes MS. At the
%  downstream side of the struct, there is black body radiation between the
%  downstream front of the membrane and a wall having the upstream temperature,
%
%    q2 = -sigma * (T1^4 - T2^4),
%
%  where q2 is the heat flux density at the downstream location in the
%  downstream direction, sigma is the Stefan-Boltzmann constant, T1 the upstream
%  and T2 the downstream temperature. The membrane struct MS is constructed with
%  MSTACKSTRUCT.
%
%  MNUMRADIATION(T1,P1,P2,SUBSTANCE,MS,'crude') sets crude solver tolerances.
%
%  [M,MS] = MNUMRADIATION(T1,P1,P2,SUBSTANCE,MS) writes the solution to MS.
%
%  Calls ASYM.
%
%  See also ASYM, MNUMADIABAT, MNUMALPHA, MT1EQT2, MNUMISO, MSTACKSTRUCT,
%  SUBSTANCE.

% Radiation heat exchange
% From
%            a2 e1 sigma T1^4 - a1 e2 sigma T2^2
%    q_12 = -------------------------------------
%                     a1 + a2 - a1 a2
%
% with a1, a2 - absortivities of walls 1, 2; e1, e2 - emmisivities.
% For small temperature differences (T1 - T2)/T1, e1 = a1 and e2 = a2, hence
%
%            sigma (T1^4 - T2^4)
%    q_12 = ---------------------
%             1/a1 + 1/a2 - 1
%

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

% Get a minimum temperature
Tmin = T2ad - max((T1 - T2ad)/10., 0.1);

if p1 == p2
    h1_T1p2 = 0;
else
    h1_T1p2 = integral(@(p) s.dhdp(T1,p), p2, p1);
end
hT1p2_h2 = ode45(@(T,h) -s.cpg(T,p2), [T1 Tmin], 0);

h12 = @(T) h1_T1p2 + deval(hT1p2_h2, T);


% inital guess
if isfield(ms,'mguess') && isscalar(ms.mguess) && ~isempty(ms.mguess)
  mguess = ms.mguess;
elseif isinf(s.ps(T1)) % above the critical temperature
  mguess = ms.mfluxknudsen(T1,p1,p2,s,ms) + ms.mfluxviscous(T1,p1,p2,s,ms);
else
  mguess = ms.mfluxliquid(T1,p1,p2,s,ms);
end

solver = solverstruct(accuracy);
m = findzero(@presiduum,mguess,(p1-p2)/1000);

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
	% m h1 + q1 = m h2 + q2, q1 = 0  ->
	% q2 = m (h1 - h2) = sigma (T2^4 - T1^4)
	% I believe, numerically it is more accurate to write
	% T1^4 - T^4 = (T1 - T)(T1 + T)(T1^2 + T^2)
	T2 = fzero(@(T) 5.67e-8*(T1 - T)*(T1 + T)*(T1^2 + T^2) + m*h12(T),...
			[Tmin T1]);
	state2.T = T2;
	state2.q = m * h12(T2);
	pres = asym(m, state2, ms, solver) - p1;
end %------------------------------------------------------------- end presiduum
end %%% END MNUMRADIATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MNUMRADIATION %%%
