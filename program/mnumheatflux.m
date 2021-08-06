function [m,ms] = mnumheatflux(negq2,T1,p1,p2,s,ms,accuracy)
%MNUMHEATFLUX   Mass flux for flow with given downstream heat flux.
%  MNUMHEATFLUX(NEGQ2,T1,P1,P2,SUBSTANCE,MS) returns the mass flux [kg/m2s] for
%  the flow of SUBSTANCE through a stack of layered membranes MS. At the
%  downstream side of the stack, a heat flux NEGQ2, which is positive if the
%  heat flows in the upstream direction, is given. T1 is the upstream
%  temperature and T2 the downstream temperature. The membrane struct MS is
%  constructed with MSTACKSTRUCT.
%
%  MNUMHEATFLUX(NEGQ2,T1,P1,P2,SUBSTANCE,MS,'crude') sets crude solver
%  tolerances.
%
%  [M,MS] = MNUMHEATFLUX(NEGQ2,T1,P1,P2,SUBSTANCE,MS) writes the solution to MS.
%
%  Calls ASYM.
%
%  See also ASYM, MNUMADIABAT, MNUMALPHA, MT1EQT2, MNUMISO, MSTACKSTRUCT,
%  SUBSTANCE.

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
state2 = downstreamstate(T1,p2,[],-negq2,s); % thus, m is empty

% Copy some values to the membrane struct;
ms.T1 = T1;
ms.p1in = p1;
ms.p2 = p2;
ms.a2 = state2.a;

% An estimate of the lower bound of the mass flux
% For R142b and isobutane, 25 nm pores, the liquid mass flux is smaller than
% the gaseous mass flux. For larger pores, the reverse becomes true.
if isfield(ms,'mguess') && isscalar(ms.mguess) && ~isempty(ms.mguess)
    mguess = ms.mguess;
elseif isinf(s.ps(T1)) % above the critical temperature
    mguess = ms.mfluxknudsen(T1,p1,p2,s,ms) + ms.mfluxviscous(T1,p1,p2,s,ms);
else
    mvapor = ms.mfluxknudsen(T1,p1,p2,s,ms) + ms.mfluxviscous(T1,p1,p2,s,ms);
    mguess = min(ms.mfluxliquid(T1,p1,p2,s,ms), mvapor);
end

% The downstram temperature under adiabatic conditions
T2ad = s.intjt(T1,p1,p2);

% Calculate and set flowsetup-structure. That needs the full range of possible
% temperatures of the process, hence provide a minimum and maximum temperature.

% With an estimate for the mass flux, define qdia, the heat flux necessary to
% force T1 = T2, which is given by
%   qdia = m (h(T1,p1) - h(T1,p2)).
% Define h1_3 = h(T1,p1) - h(T1,p2), hence  qdia = m * h1_3.  Note, qdia < 0.
% Estimate the temperature difference by taking cpg constant,
%   q2 = m (h(T1,p1) - h(T2,p2)) = m (h(T1,p1) - h(T1,p2) + h(T1,p2) - h(T2,p2))
%   q2 = qdia + m cpg (T1 - T2).
% For q2 < qdia follows T2 > T1.
%   -negq2 = m h1_3 - m cpg (T2-T1).
% Otherwise, for q2 > 0, i.e., cooling from downstreams,
%   q2 = m cpg (T2ad - T2),
%   -negq2 = m cpg (T2ad - T2)

if p1 == p2
    h1_3 = 0;
    qdia = 0;
else
    h1_3 = integral(@(p) s.dhdp(T1,p), p2, p1);
    qdia = mguess * h1_3;   % qdia < 0
end

% Add some safety margin to the minimal temperature.
% Since writeflowsetups() anyhow doubles to requested temperature range,
% extending beyond the maximum temperature, Tmax does not need to be extended.
if negq2 < 0    % heat flux in downstream direction = cooling from downstreams
    Tmax = T1;
    Tmin = T2ad + negq2 / (mguess * s.cpg(T2ad,p2));
elseif negq2 < -qdia   % a bit of heating from downstreams
    Tmax = T1 + (T1 - T2ad);
    Tmin = T2ad - 0.1 * (T1 - T2ad);
else % negq2 > -qdia, a lot of heating from downstreams
    Tmax = T1 + (h1_3 + negq2/mguess) / s.cpg(T1,p2);
    Tmin = T2ad - 0.1 * (T1 - T2ad);
end

% Otherwise, a temperature difference given by heat conduction would be, for
% a membrane consisting of one layer,
%deltaTcond = negq2 * ms.membrane(1).layer(1).matrix.L / ...
%    ms.membrane(1).layer(1).fmodel.kmgas(...
%                ms.membrane(1).layer(1).matrix.epsilon,...
%                ms.membrane(1).layer(1).matrix.km,...
%                s.kg(T1));
% OR, even worse, compute (large) temperature difference for heat conduction
% through the gas, which presumably has the worst thermal conductivity,
%deltaTcond = negq2 * ms.membrane(1).layer(1).matrix.L / s.kg(T1);
% However, this corresponds to q1 = q2 and does not fulfil the boundary
% condition q1 = 0. Also, the energy equation, q1 + m h1 = q2 + m h2, is only
% fulfilled for  m = 0, or h1 = h2.

% writeflowsetups() by itself adds a considerable range above T1,
ms = ms.writeflowsetups(Tmax,Tmin,s,ms);

% For computing the enthalpy below, now extend the maximum temperature.
Tmax = Tmax + (Tmax - Tmin);

% From the global energy balance,
%
%   m h1 + q1 = m h2 + q2
%
% state 2 and the heat flux q2 are not independent,
%
%   h2 = f(m,h1,q2)   if q1 = 0.
%
% Therefore, given q2, the downstream state must be determined,
% depending on the mass flux density m.
%
% Construct a function h1_2(T2) that gives the enthalpy difference h1 - h2,
% h(T1,p1) - h(T2,p2).
%
% Consider a T-s diagram,
%     T                          . p = const
%      ^                       .
% Tmax-|-                    . -- hmax
%      |                   .
%      |                 .  anywhere between: h(T,p2) = h2
%      | h(T1,p1) = h1 .
%  T1 -|- +----------+ h(T1,p2) = h3
%      |   `       .
%      |     `   .
%  T2 -|-      `  h(T2,p2) = h4
%      |        isenthalpic line
%      +--------------> s
%
%
% Use the following symbols,
%
%   h1 = h(T1,p1),  h3 = h(T1,p2),  h4 = h(T2,p2),  h2 = h(T,p2)
%   h1_3  = h(T1,p1) - h(T1,p2) = int_p2^p1 dhdp(T1,p) dp,
%   h3_4 = h(T1,p2) - h(T2,p2) = int_T2^T1 cp(T,p2) dT,
%   h4_2 = h(T2,p2) - h(T,p2) = -int_T2^T cp(T',p2) dT',
%
% Note, h1_3 < 0, h3_4 > 0. Also, dhdp < 0, cp > 0.
% Really, below use Tmin instead of T2.
% h1_3 is evaluated further above.

h3_4 = integral(@(T) s.cpg(T,p2), Tmin, T1);

h3_2 = ode45(@(T,h) -s.cpg(T,p2), [Tmin Tmax], h3_4);

h1_2 = @(T) h1_3 + deval(h3_2, T);

solver = solverstruct(accuracy);
try
    m = findzero(@presiduum,mguess,(p1-p2)/1000);
catch ME
    if strcmp(ME.identifier, 'MATLAB:fzero:ValuesAtEndPtsSameSign')
        error(['No interval: mguess = %.3g, mvapor = %.3g, q = %.0f,' ...
                ' T = [%.2f %.2f], m*h1_2[max, min] = [%.2g %.2g]\n'], ...
                mguess, mvapor, negq2, Tmax-273.15, Tmin-273.15, ...
                mguess*h1_2(Tmax), mguess*h1_2(Tmin));
    else
        rethrow(ME);
    end
end

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
    T2 = fzero(@(T) negq2 + m * h1_2(T), [Tmax Tmin]);
    state2.T = T2;
    pres = asym(m, state2, ms, solver) - p1;
end %------------------------------------------------------------- end presiduum

end %%% END MNUMHEATFLUX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MNUMHEATFLUX %%%
