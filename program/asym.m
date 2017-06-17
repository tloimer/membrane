function [p1,ms] = asym(m,state,ms,solver)
%ASYM       Upstream pressure for adiabatic flow through a MSTACKSTRUCT.
%  ASYM(M,STATE,MSTACKSTRUCT,SOLVER) returns the upstream pressure for the
%  mass flux M of the fluid MSTACKSTRUCT.SUBSTANCE through an assemblage of
%  membranes described by the struct MSTACKSTRUCT. The downstream state of
%  the fluid is given by STATE, the solution tolerances are controlled via
%  the struct SOLVER.
%
%  [P1,MSTACKSTRUCT] = ASYM(M,STATE,MSTACKSTRUCT,SOLVER) returns the
%  upstream pressure P1 and a MSTACKSTRUCT describing the solution.
%
%  Restrictions: The upstream state must be a gaseous phase, because the
%  liquid film and the temperature boundary layer integrators upstream of
%  the membrane (ifreevapor, ifreeliquid) use termination conditions that
%  work only with an gaseous upstream state.
%  The vapor phase integrator (integratevapor) always starts at the
%  downstream end of a layer.
%  For a two-phase upstream state, integration direction would have to be
%  downstream.
%
%  ASYM>FRONT should probably be a subfunction or nested function in STATE,
%  STATE1 = STATE2.FRONT(STATE2,...). However, algorithms in FRONT and
%  INTEGRATE are partially similar, therefore they are kept in one place.
%
%  See also DOWNSTREAMSTATE, MSTACKSTRUCT, FLOWSETUP, SOLVERSTRUCT.

if solver.writesolution
  % construct empty (0x0) struct flow.
  % Use of [] instead of {} would create a 1x1 struct flow.
  zeroflowstruct = struct('z',{},'T',{},'p',{},'a',{},'q',{},'Kn',{},'color',{});
else
  zeroflowstruct = [];
end

% This part of the code must know the mstackstruct-structure.
nmembranes = length(ms.membrane);
s = ms.substance;
freesetup = ms.freesetup; % free space flow setup

% integrate in upstream direction
% Cycle over all membranes
for i = nmembranes:-1:1
  nlayers = length(ms.membrane(i).layer);

  % Cycle over all layers in one membrane
  for j = nlayers:-1:1
    % Start is the downstream end of the last layer of the last membrane;
    % Then, at the end of each layer.
    z = ms.membrane(i).layer(j).matrix.L;
    % flow grows in integrate(), if solver.writesolution is true
    flow = zeroflowstruct;
    % Stay in one layer, as long as it is not ready.
    while z > 0
      state = front(state, ms.membrane(i).layer(j).flsetup, m, s);
      % [state, z, flow] = integrate(state, z, matrix, flsetup, m, s, solver);
      [state,z,flow] = integrate(state,z,flow,ms.membrane(i).layer(j).matrix,...
				 ms.membrane(i).layer(j).flsetup, m,s,solver);
    end
    % Arrived at the upstream end of a layer
    ms.membrane(i).layer(j).flow = flow;
  % Go into the next layer
  end

  % Upstream front of a membrane; Calculate the state in the intermediate space.
  state = front(state,freesetup,m,s);
  % calculate the temperature far upstream
  % but not necessarily for the upstream-most layer of the upstream-most membrane
  % if i == 1 && solver.partialsolution, else  ... integratefree(..);  end
  % TODO: here, rewrite by checking for q == 0
  if i > 1 || solver.fullsolution
    flow = zeroflowstruct;
    % this might either be a liquid film, the gaseous temperature boundary
    % layer, or two-phase flow;
    [state,z,flow] = integratefree(state,z,flow,m,s,solver);
    while state.q ~= 0
      state = front(state,freesetup,m,s);
      [state,z,flow] = integratefree(state,z,flow,m,s,solver);
    end
    % if the last call to integratefree called ifreevapor, zscale is returned
    if state.phase == 'g'
      ms.membrane(i).zscale = z;
    end
    ms.membrane(i).flow = flow;
  end

% Next membrane
end
if isempty(state.p)	% equivalent to: state.phase == '2'
  p1 = state.pk;
else
  p1 = state.p;
end

end %%% END ASYM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END ASYM %%%


%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%

function state1 = front(state2,fs,m,s) %---------------------------------- front
%FRONT      Return the state of the fluid at the upstream side of a front.
%
%  STATE1 = FRONT(STATE2,MEMBRANE,FLSETUP,M,SUBSTANCE) returns the state STATE1
%  at the upstream side of a front. STATE1 is calculated from the downstream
%  state STATE2, the membrane layer at the upstream side of the front, the mass
%  flux M and a SUBSTANCE. The properties of the membrane layer are given by
%  MEMBRANE and FLSETUP. STATE2 conveys sufficient information of the downstream
%  side of the front, such that the properties of the downstream layer do not
%  need to be given. The calculation runs in upstream direction, against the
%  flow direction.
%
%                        front
%                          |
%      MEMBRANE, FLSETUP   |   STATE2
%            STATE1        |              ---->  flow direction
%                          |              <----  direction of calculation
%
%Input
%  STATE2 must contain
%    STATE2.phase (one of 'g', 'l' or '2'),
%    STATE2.T and the utility functions
%    STATE2.avapor, STATE2.aliquid and STATE2.atwophase.
%  In addition, if STATE2.phase is 'g' or 'l', STATE2 must contain
%    STATE2.p,
%    STATE2.a,
%    STATE2.q.
%  If STATE2.phase is '2', STATE2 must in addition contain
%    STATE2.q_mh,
%    STATE2.hvapK,
%    STATE2.pk,
%    STATE2.pliq.
%  Note, that STATE2.a is not needed, STATE2.q_mh and STATE2.hvapk is used
%  instead.
%
%Output
%  STATE1 contains,
%    STATE1.phase,
%    STATE1.T,
%    STATE1.avapor, STATE1.aliquid and STATE1.atwophase.
%  If STATE1.phase is 'g' or 'l', STATE1 contains in addition
%    STATE1.p,
%    STATE1.a,
%    STATE1.q.
%  If STATE1.phase is '2', STATE1 contains in addition
%    STATE1.q_mh,
%    STATE1.hvapK,
%    STATE1.pk,
%    STATE1.dpk,
%    STATE1.dpcap.
%
%  See also DOWNSTREAMSTATE, ASYM.

% Some abbreviations, for simplicity.
T2 = state2.T;
p2 = state2.p;
a2 = state2.a;
q2 = state2.q;

% A possibility, to define mustbecomegaseous and mustbecomeliquid for any phase
%paux1 = pk1 - (1-a2)*pcap1;
%	    pk1		for state2 = gaseous
%  puax1 =  pk1 - pcap1	for state2 = liquid
%	    in between	for state2 = two-phase mixture
%mustbecomegaseous = p2 < paux1;
%mustbecomeliquid = p2 > paux1;
% not(mustbecomegaseous) && not(mustbecomeliquid) = p2 == paux1

[pk1, pcap1] = fs.pkpcap(T2);

switch state2.phase
  case 'g' % gaseous
    % First deal with the signal p2 == pk1 (within computer accuracy). This
    % signal is only triggered within a membrane layer.
    if state2.vapliqequilibrium % p2 == pk1, to within computer accuracy
      [canbecomevapor,canbecomeliquid,hvapK1,dpk1,dpcap1] = heatfluxcriterion;
      if canbecomeliquid
	interface_liqvap;
      else
	interface_2phvap;
      end
    else
      % the logic based on numerically correct comparisons follows
      if p2 < pk1 % mustbecomegaseous
	% nointerface
	state1 = state2;
      elseif p2 > pk1 % mustbecomeliquid
	interface_liqvap;
      else % p2 == pk1, gaseous, liquid or two-phase possible
	% pk1 and pcap1 are calculated again - this does not happen too often, anyway
	[canbecomevapor,canbecomeliquid,hvapK1,dpk1,dpcap1] = heatfluxcriterion;
	if canbecomevapor
	  % nointerface
	  state1 = state2;
	  if canbecomeliquid
	    warning('Double solution possible. Here, vapor remains vapor.');
	  end
	elseif canbecomeliquid
	  interface_liqvap;
	else
	  interface_2phvap;
	end
      end
    end

  case 'l' % liquid
    % Deal with the signal p2 == pk1 -pcap1. This may be true only within
    % computer accuracy, and it is raised only within a membrane layer.
    if state2.vapliqequilibrium % p2 == pk1 - pcap1, to within computer accuracy
      [canbecomevapor,canbecomeliquid,hvapK1,dpk1,dpcap1] = heatfluxcriterion;
      if canbecomevapor
	interface_vapliq;
      else
	interface_2phliq;
      end
    else
      % Signal exception is treated, continue with correct logic
      if p2 > pk1 - pcap1 % mustbecomeliquid
	% nointerface
	state1 = state2;
      elseif p2 < pk1 - pcap1 % mustbecomegaseous
	interface_vapliq;
      else % gaseous, liquid or two-phase possible
	% pk1 and pcap1 are calculated again - this does not happen too often, anyway
	[canbecomevapor,canbecomeliquid,hvapK1,dpk1,dpcap1] = heatfluxcriterion;
	if canbecomeliquid
	  % nointerface
	  state1 = state2;
	  if canbecomevapor
	    warning('Double solution possible. Here, liquid remains liquid.');
	  end
	elseif canbecomevapor
	  interface_vapliq;
	else
	  interface_2phliq;
	end
      end
    end

  case '2' % twophase
    if state2.pk > pk1 % mustbecomeliquid
      interface_liq2ph;
    elseif state2.pk < pk1 % mustbecomevapor
      interface_vap2ph;
    else % state2.pk == pk1
      % The two-phase heatfluxcriterion
      [qmin,qmax,hvapK1,dpk1,dpcap1] = fs.qminqmax(m,T2);
      % qvapor would be state2.q_mh - m*hvapK1, because qvap + m*hvapK = q_mh
      % similarly: qliq = state2.q_mh ( hvapK1 == state2.hvapK )
      % Allow a tad negative q1 for vapor flow in free space. The negative q1
      % corresponds to 1-dotx < 0.001. See also case '2' in integratefree below
      % dotx = state.q_mh/(m*state.hvapK);
      % q_mh = q + m*dotx*hvapK, must be smaller or equal to m*hvapK
      if state2.q_mh - m*state2.hvapK >= qmin || state2.q_mh/(m*state2.hvapK) > 0.999
	interface_vap2ph;
      elseif state2.q_mh <= qmax
	interface_liq2ph;
      else % two-phase flow
	% nointerface
	state1 = state2;
	state1.dpk = dpk1;
	state1.dpcap = dpcap1;
	% if, during iteration, there is two-phase flow upstream of the membrane
	% (which is possible for theta = 90), do not leave state.p empty - it
	% will be passed to presiduum, eventually, and findzero fails.
	% For theta = 90, p = pk is even correct.
	state1.p = state2.pk;
      end
    end

  otherwise
    error('How can we come here?');
end

%--- nested functions ----------------------------------------- nested functions

function [canbecomevapor,canbecomeliquid,hvapK1,dpk1,dpcap1] = heatfluxcriterion
%HEATFLUXCRITERION Return whether a fluid can become vapor or liquid upstreams.
  [qmin,qmax,hvapK1,dpk1,dpcap1,pcap1] = fs.qminqmax(m,T2);
  % canbecomevapor = q1 >= qmin    ---> flow   vapor 1 | state 2 <--- compute
  % canbecomeliquid = q1 <= qmax   ---> flow  liquid 1 | any state 2
  % The heat flux q1 is, depending on the phase change tabulated below,
  % vap  vap   q1 = q2
  % vap  liq   q1 = q2 - m*hvapK
  % liq  vap   q1 = q2 + m*hvapK
  % liq  liq   q1 = q2
  % allow a tad negative q1, which is necessary for the flow in free space,
  % with the  -0.001*m*hvapK1  term!
  canbecomevapor = q2 - (1-a2)*m*hvapK1 >= qmin - 0.001*m*hvapK1;
  % A liquid is possible for q1 < qmax, the last two lines above give
  canbecomeliquid = q2 + a2*m*hvapK1 <= qmax;
  % For complete phase change, xdot is not necessary, a2 can be used instead.
end %---------------------------------------------------------------------------

% Convention: pcap = p_vapor - p_liquid.
% heat balance: q_mh = q + m h = const.

function interface_liqvap %-----------------------------------------------------
  %
  % -------------
  %     liq.  1 ( 2  vap.
  % -------------
  %
  % The radius of curvature can be explicitly calculated because
  % p2 = pk(T,r). Hence, r is immediately known. See flow12>front62.
  %   p2 = pk(T1,r),
  %   pcap = p2 - p1.  (pcap > 0)
  % With Kelvin's equation,
  %   ln(p2/psat) = -pcap/(rho*R*T),
  %   pcap = ln(psat/p2)*rho*R*T.
  psat = s.ps(T2);
  [drho2, rho2] = s.drho(T2);
  pcap = log(psat/p2)*rho2*s.R*T2;
  % arguments: hvapKraw(T,pk,psat,pcap,rho,drho)
  hvapK12 = fs.hvapKraw(T2,p2,psat,pcap,rho2,drho2);
  p1 = p2 - pcap;
  q1 = q2 + m*hvapK12;
  state1 = state2.aliquid(T2,p1,q1);
end %---------------------------------------------------------------------------

function interface_vapliq %-----------------------------------------------------
  %
  % -----------------
  %	    vap.  1 ) 2  liq.
  % -----------------
  %
  % calculate r, radius of curvature - but not directly, only pk
  % copied and modified from flow12>front35
  %   p1 = pk(T1,r),
  %   p1 = p2 + pcap  (pcap < 0, hence p1 < p2),
  % Evaluate Kelvin's eq.,
  %   ln(p1/psat) = -pcap/(rho*R*T),
  %   ln(p1/psat) = (p2-p1)/(rho*R*T).
  % Solve the eq. above by Newton iteration,
  %   F(x) = ln(x) + (x*psat - p2)/rho*R*T,
  %   F'(x) = 1/x + psat/rho*R*T.
  % A good first guess is, with ln(p1/psat) ~ p1/psat - 1,
  %   p1/psat - 1 = p2/(rho*R*T) - p1/(rho*R*T),
  %   ( with R = rho*R*T, p1*(1/psat + 1/R) = p2/R + 1, )
  %   ( p1 = (p2/R + 1)*psat*R/(R + psat)		   )
  %   p1/psat = (p2 + rho*R*T) / (psat + rho*R*T).
  psat = s.ps(T2);
  [drho2, rho2] = s.drho(T2);
  rhoRT = rho2*s.R*T2;
  % initial guess; here, p1 = p1/psat
  p1 = (p2 + rhoRT)/(psat + rhoRT);
  % setup newton
  ps_rhoRT = psat/rhoRT;
  function [F, dF] = sol35(x)
    F = log(x) + x*ps_rhoRT - p2/rhoRT;
    dF = 1/x + ps_rhoRT;
  end
  p1 = newton(@sol35,p1,1e-15) * psat;
  hvapK12 = fs.hvapKraw(T2,p1,psat,p1-p2,rho2,drho2);
  % hvapKraw(T,pk,psat,pcap,rho,drho)
  q1 = q2 - m*hvapK12;
  state1 = state2.avapor(T2,p1,q1);
end %---------------------------------------------------------------------------

function interface_2phliq %-----------------------------------------------------
  % Solve for q1 and a1 in the integrator, not here. Also, pass a few variables
  % (dpk, dpcap, ...), so they need not calculated twice.
  q_mh12 = q2;
  state1 = state2.atwophase(T2,q_mh12,hvapK1,pk1,pcap1,dpk1,dpcap1);
  % set a (dummy) pressure, if the liquid integrator just stops at the upstream
  % membrane front, and theta = 90; Then, this pressure is even correct.
  % state1.p = pk1;
end %---------------------------------------------------------------------------

function interface_2phvap %-----------------------------------------------------
  % Solve for q1 and a1 in the integrator. Two-phase flow is fully determined by
  % the temperature and flux of enthalpy.
  q_mh12 = q2 + m*hvapK1; % = q2 + m*hvapK2, because p2 = pk1.
  state1 = state2.atwophase(T2,q_mh12,hvapK1,pk1,pcap1,dpk1,dpcap1);
  % set a (dummy) pressure, if the vapor integrator just stops at the upstream
  % membrane front, and theta = 90; Then, this pressure is even correct.
  % state1.p = pk1;
end %---------------------------------------------------------------------------

function interface_liq2ph %-----------------------------------------------------
  % in general, here p1 = p2liq; because, there are pores where a liquid part
  % of the 2ph-region might border to the liquid in part1. liq-liq, so, no
  % pressure difference
  p1 = state2.pliq;
  q1 = state2.q_mh;
  state1 = state2.aliquid(T2,p1,q1);
end %---------------------------------------------------------------------------

function interface_vap2ph %-----------------------------------------------------
  % in general, here p1 = p2vap; for reason, see interface_liq2ph.
  % For memfree: pk2 = p2, pcap2 = 0, hvapk2 = s.hvap(T2).
  p1 = state2.pk;
  q1 = state2.q_mh - m*state2.hvapK;
  state1 = state2.avapor(T2,p1,q1);
end %---------------------------------------------------------------------------

end %----------------------------------------------------------------- end front

%--------------------------------------------------------------------- integrate
function [state,z,flow] = integrate(state,z,flow,matrix,flsetup,m,s,solver)
%INTEGRATE  Integrate the flow within the layer.

% integrate knows the state-struct. The specific integrators do not need to know
% about the state-struct.
switch state.phase
  case 'g' % gaseous
    [T,p,q,z,flow,state.vapliqequilibrium] = integratevapor(...
	m, state.T, state.p, state.q, state.vapliqequilibrium, flow, matrix,flsetup,s,solver);
    % Above, state.vapliqequilibrium is constructed and, possibly, set; Another
    % possibility is:   if z > 0, state.vapliqequilibrium = true; end
    % Update the state.
    % Another possibility is to construct a new state, state=state.avapor(T,p,q)
    state.T = T;
    state.p = p;
    state.q = q;
  case 'l' % liquid
    [T,p,q,z,flow,state.vapliqequilibrium] = integrateliquid(...
	m, state.T, state.p, state.q, z, state.vapliqequilibrium, flow, matrix,flsetup,s,solver);
    % Update the state.
    state.T = T;
    state.p = p;
    state.q = q;
  case '2' % two-phase
    [T,q_mh,hvapK,pk,pliq,z,flow] = integratetwophase(m,state.T,state.q_mh, ...
     state.hvapK,state.pk,state.dpk,state.dpcap,z,flow,matrix,flsetup,solver);
    state.T = T;
    state.q_mh = q_mh;
    state.hvapK = hvapK;
    state.pk = pk;
    state.pliq = pliq;
    state.dpk = [];
    state.dpcap = [];
  otherwise
    error('No phase letter?');
end
end %------------------------------------------------------------- end integrate

function [state,z,flow] = integratefree(state,z,flow,m,s,solver) % integratefree
switch state.phase
  case 'g' % gaseous
    % the thermal boundary layer.
    [T,p,q,z,flow] = ifreevapor(m,state.T,state.p,state.q,z,flow,s,solver);
    % Update the state.
    state.T = T;
    state.p = p;
    state.q = q;
  case 'l' % liquid
    % liquid film. z is not given, must be zero.
    [T,p,q,z,flow] = ifreeliquid(m,state.T,state.p,state.q,flow,s,solver);
    % Update the state.
    state.T = T;
    state.p = p;
    state.q = q;
  case '2' % two-phase
    % two-phase in front - occured, as long as the heatfluxcriterion did not
    % allow   q1 >= qmin - minute correction   (instead of q1 >= qmin.)
    dotx = state.q_mh/(m*state.hvapK);
    % 1 - dotx < 0.001 (dotx > 0.999) should have been caught in front, case '2'
    warning('Two-phase flow in free space, 1 - dotx = %.3g, should be zero.',...
	    1 - dotx);
    state.q = 0; % This is the signal for the upstream front code to terminate.
    % state.p was set in the front, state.p = state.pk.
    % This should be:
    % f = fmodel('plug'); a = f.a(dotx,s.v(state.T,..),1/s.rho(...))
    % homogeneous flow model: inverse function to f.xdot
    a = dotx./(1 + dotx.*(1-1./(s.rho(state.T).*s.v(state.T,state.pk))));
    flow = writeflow(flow,{'z','T','p','a','q','color'},...
			  {0,state.T,state.p,a,state.q,'g'});
  otherwise
    error('Can not integrate phase %c in free space.', state.phase);
end
end %--------------------------------------------------------- end integratefree

%---------------------------------------------------------------- integratevapor
function [T9,p9,q9,z9,flow,vapliqequilibrium] = integratevapor(m,T2,p2,q2,...
					vapliqequilibrium,flow,mem,fs,s,solver)
%INTEGRATEVAPOR Vapor flow within the membrane - a copy of FLOW12>FLOW92

% TODO: z is not referred to! It is assumed z = mem.L.

% Termination condition: p = pk(T), valid: p < pk(T).
% integrate from z = L to z = 0, initial conditions T2, p2, q2.
% eqs.: Darcy's law, Fourier, energy: mh + q = const => m dh/dz = -dq/dz,
% dh/dz = (dh/dT)_p dT/dz + (dh/dp)_T dp/dz, dh/dT = s.cpg, dh/dp = s.dhdp
%   dp/dz = -m nu/kappa
%   dT/dz = -q/k
%   dq/dz = -m cpg dT/dz - m dh/dp dp/dz
% Note: jt = -(dh/dp)/cpg
% dimensionless (w, instead of star; with zscale = mem.L):
%   zscale = mem.L, pscale = ps(T2)-p2,  qscale = -m dhdp2 pscale.
% But really, pw varies to O(1) not on zscale but on a fraction
%   deltazw = pscale*mem.kappa/(m*nu2*L).
% We enforce zscale = mem.L, which is not appropriate to the problem, are
% content with dpw/dzw >> 1 but care for the temperature variation to be O(1):
%   Tscale = qscale*zscale*deltazw/km.
%   z = L*zw;  p = pscale*pw + p2;  T = Tscale*Tw + T2;  q = qscale*qw;
% dimensionless eqs:
% pscale/zscale dpw/dzw = -m nu/kap,
%   dpw/dzw = -m*zscale/kappa*pscale nu,  O(pw) = 1, O(dpw/dzw) = m*zscale/...
% Tscale/zscale dTw/dzw = -qw*qscale/k,
%   dTw/dzw = -km/deltazw qw/k
% qscale dqw/dzw = -m (cpg Tscale dTw/dzw + dh/dp pscale dpw/dzw)  |  eq./zscale
%   dqw/dzw = - cpg m*zscale*deltazw/km dTw/dzw + dhdp/dhdp2 dpw/dzw

%  CHARACTERISTIC SCALES
% scales to make dimensionless
%pdarcy = m*s.nug(T2,p2)*mem.L/mem.kappa;
pdarcy = m*fs.nuapp(T2,p2)*mem.L/mem.kappa;
pscale = min(s.ps(T2) - p2,pdarcy);
[dhdp2 cpg2] = s.dhcpg(T2,p2);
%dhdp2 = s.dhdp(T2,p2);
qscale = -m*dhdp2*pscale;
deltazw = pscale/pdarcy; % = 1 if pscale = pdarcy
Tscale = qscale*mem.L*deltazw/mem.km;
step92 = -min(solver.odemaxstep(pscale,solver.maxpperstep),solver.maxpperstep/pdarcy);
% The length scale on which the temperature varies is k/(m cp),
% This was necessary for large mass fluxes, e.g., the solution for the flow of
% nitrogen through the cotton plug used by Thomson (1853).
init92 = -min(-step92, fs.kmgas(T2)/(mem.L*m*cpg2));
% functions to re-calculate dimensional values
mkTdim = @(Tw) Tscale.*Tw + T2;
mkpdim = @(pw) pscale.*pw + p2;
function [T,p,q,z] = mkdimensional(Tw,pw,qw,zw)
  T = mkTdim(Tw); p = mkpdim(pw);
  q = qscale.*qw; z = mem.L.*zw;
end
% coefficients for the eqs. in int92w:
coeff1 = -m*mem.L/(mem.kappa*pscale);
coeff2 = -mem.km/deltazw;
coeff3 = m*mem.L/coeff2;

%  INTEGRATE
% make initial conditions dimensionless; integrate
%T2w = 0; p2w = 0; z2w = 1;
q2w = q2/qscale; % q2w = 0;
options=odeset('Events',@term92w,...%'AbsTol',atol*[1 1 -coeff1*cpg2],...
  'RelTol',solver.rtol,'Refine',1,'InitialStep',init92,'MaxStep',step92);
%sol92 = ode45(@int92w,[z2w 0],[T2w p2w q2w],options);
sol92 = ode45(@int92w,[1 0],[0 0 q2w],options);

%function dy = int92(z,y)
%dimensional. Equations see above.
%T = y(1);  p = y(2);  q = y(3);
%dp = -m*s.nug(T,p)/mem.kappa;  dT = -q/kmgas(T);
%[dhdp cpg] = s.dhcpg(T,p);
%dq = -m*( cpg*dT + dhdp*dp );
%dy = [dT;dp;dq];
%end

function dy = int92w(z,y)
  % dimensionless eqs., cf. above.
  Tw = y(1);  pw = y(2);  qw = y(3);
  T = mkTdim(Tw); p = mkpdim(pw);
  dpw = fs.nuapp(T,p)*coeff1;
  dTw = coeff2*qw/fs.kmgas(T);
  [dhdp, cpg] = s.dhcpg(T,p);
  dqw = coeff3*cpg*dTw + dpw*dhdp/dhdp2;
  dy = [dTw; dpw; dqw];
end

function [val,isterm,direction] = term92w(z,y)
  % Terminate integration when pk(T)-p = 0, only when falling. OK: pk(T)-p > 0.
  % without direction, integration already terminates at start
  isterm = 1; direction = 0; %-1;
  % T = y(1), p = y(2)
  % dimensional: val=pk(y(1))-y(2);
  % dimensionless:
  val = fs.pkelv(mkTdim(y(1))) - mkpdim(y(2));
end

%  ASSIGN LAST POINT
% dimensionalize and assign actual state, here T9, p9, q9, z9
last = size(sol92.x,2);
[T9, p9, q9, z9] = mkdimensional(sol92.y(1,last),sol92.y(2,last),...
   sol92.y(3,last),sol92.x(last));
% not necessary: .ye, .xe are also in .y(last), .x(last).
%  [T9 p9 q9 z9] = mkdimensional(sol92.ye(1,end),sol92.ye(2,end),...
%    sol92.ye(3,end),sol92.xe(end));

% but necessary: RAISE A SIGNAL IF INTEGRATION TERMINATED PREMATURELY
% This means, p9 == pK, but because of mkdimensional() within computer accuracy
% In contrast to what the documentation makes believe, the field .ie is always
% present, but empty when no termination event occured.
if ~isempty(sol92.ie) % integration terminated by event function term92w
  vapliqequilibrium = true;
end

%  WRITE SOLUTION
% if wanted, write the solution
if solver.writesolution
  % allocate space for all points; assign last point
  T92(last) = T9; p92(last) = p9; q92(last) = q9; z92(last) = z9;
  Kn92(last) = fs.knudsen(T9,p9);
  % and write the solution but the last point
  last = last - 1;
  [T92(1:last), p92(1:last), q92(1:last), z92(1:last)] = mkdimensional(...
    sol92.y(1,1:last),sol92.y(2,1:last),sol92.y(3,1:last),sol92.x(1:last));
  for i = 1:last
    Kn92(i) = fs.knudsen(T92(i),p92(i));
  end
  flow = writeflow(flow,{'z','T','p','a','q','Kn','color'},...
			{z92,T92,p92,ones(1,last+1),q92,Kn92,'r'});
end

end %-------------------------------------------------------- end integratevapor

%--------------------------------------------------------------- integrateliquid
function [T5,p5,q5,z5,flow,vapliqequilibrium] = integrateliquid(m,T6,p6,q6,z6,...
					vapliqequilibrium,flow,mem,fs,s,solver)
%INTEGRATELIQUID Liquid flow within the membrane - a copy of FLOW12>FLOW56.
% Termination condition: p = pk(T) - pcap, valid: p > pk(T) - pcap,

% integrate from z = z6 to z = 0, initial conditions T6, p6, q6.
% eqs.: Darcy's law, Fourier, energy: mh + q = const => m dh/dz = -dq/dz,
% Because cpl is mostly given as independent of p, but we would like to satisfy
% the integrability condition on d^2h/dTdp at least within the liquid region, we
% set dh/dT(T,p) = cpl(T,psat) + psat_int^p d^2h/dTdp' dp'. Then, with d^2h/dTdp
% = -T d^2v/dT^2, which is a function of T only, we can integrate dh, instead of
% using the second, presumably quite inaccurate, derivative of rho(T). Hence,
%   m h + q = m h6 + q6, q = q6 - m (h - h6)
% integrate dh; abbreviation: g(T) = dh/dp = v - T*dv/dT.
%   dh/dT(T,p) =  cpl(T,psat) + psat_int^p' dg/dT dp' = cpl(T) + dg/dT*(p-psat)
%   dh = cpl(T) dT + dg/dT*(p-psat(T)) dT + g dp
%     = cpl dT + dg p - dg/dT psat dT  + g dp
%     = cpl dT + d(g*p) - d(g*psat) + g dpsat/dT dT. Integration yields
%   h - h6
%     = T6_int^T cpl dT + g*(p-psat) - g6*(p6-psat6) - T6_int^T g dpsat/dT dT
% (Verify by taking the definite integrals, probably along different paths.)
% The last term above must be evaluated. This can be done once in the preamble.
% Equations:
%   q = q6 - m (h - h6)
%   dp/dz = -m nul/kappa
%   dT/dz = -q/kmliq
% dimensionless:
%   z = z6*zw;  p = (z6*m*nul6/kappa)*pw + p6; nul6 = s.nul(T6);
%   q = q6*qw;  T = (q6*z6/kmliq6)*Tw + T6;
% dimensionless eqs:
%   qw = 1 - m (h - h6)/q6;
%   dpw/dzw = -s.nul(T)/nul6;
%   dTw/dzw = -qw*kmliq6/kmliq;

%  CHARACTERISTIC SCALES
% scales to make dimensionless
nul6 = s.nul(T6);  kmliq6 = fs.kmliq(T6);
pscale = z6*m*nul6/mem.kappa;
Tscale = q6*z6/kmliq6;
% zscale = z6; qscale = q6;
% integration steps
% liquid properties only depend on temperature, not on pressure
% step56 = -min(odemaxstep(Tscale,maxTperstep),odemaxstep(pscale,maxpperstep))
step56 = -solver.odemaxstep(Tscale,solver.maxTperstep);
% functions to re-calculate dimensional values
mkTdim = @(Tw) Tscale.*Tw + T6;
mkpdim = @(pw) pscale.*pw + p6;
% coefficients for the eqs. in int56w:
% primarily, for the invariable part of (q6 + m*h6)/q6:
[drho6, rho6] = s.drho(T6);
coeff1 = q6 + m*( (p6-s.ps(T6))*(1+T6*drho6)/rho6 - fs.intdhdpdpsatdT(T6) );
coeff2 = -kmliq6/q6;

% calculate q. This function is needed twice, once during integration and later
% for evaluation of the solution.
function q = calcq(T,p)
  % do not vectorize! s.ps is scalar, gives wrong values on vectors.
  [drho, rho] = s.drho(T);
  q = coeff1 - m ...
    * (s.intcpl(T6,T) + (p-s.ps(T))*(1+T*drho)/rho - fs.intdhdpdpsatdT(T));
end

%  INTEGRATE
% integrate; dimensionless initial conditions
%T6w = 0; p6w = 0; q6w = 1; z6w = 1;  % O(dTw/dzw) = 1; O(dpw/dzw) = 1;
options=odeset('Events',@term56w,'RelTol',solver.rtol,...
  'Refine',1,'InitialStep',step56,'MaxStep',step56);
%sol56 = ode45(@int56w,[z6w 0],[T6w p6w],options);
sol56 = ode45(@int56w,[1 0],[0 0],options);

function dy = int56w(z,y)
  T = mkTdim(y(1));  p = mkpdim(y(2));
  dTw = calcq(T,p)*coeff2/fs.kmliq(T);
  dpw = -s.nul(T)/nul6;
  dy = [dTw; dpw];
end

function [val,isterm,direction] = term56w(z,y)
  % Termination condition: p = pk(T) - pcap, valid: p > pk(T) - pcap, pcap =
  % curv*sigma.
  % Terminate integration when p + pcap - pk = 0, falling.
  isterm = 1; direction = -1;
  % Tw = y(1);  pw = y(2);
  T = mkTdim(y(1));
  [pk, pcap] = fs.pkpcap(T);
  val = mkpdim(y(2)) + pcap - pk;
end

%  ASSIGN LAST POINT
% dimensionalize and assign actual state, here T5, p5, q5, z5
last = size(sol56.x,2);
T5 = mkTdim(sol56.y(1,last));  p5 = mkpdim(sol56.y(2,last));
q5 = calcq(T5,p5); z5 = sol56.x(last)*z6; % z5 = 0;

% RAISE A SIGNAL IF INTEGRATION TERMINATED PREMATURELY, see integratevapor
if ~isempty(sol56.ie) % integration terminated by event function term56w
  vapliqequilibrium = true;
end

%  WRITE SOLUTION
% if wanted, write the solution
if solver.writesolution
  % allocate space for all points; assign last point
  T56(last) = T5; p56(last) = p5; q56(last) = q5; z56(last) = z5;
  % and write remaining solution
  last = last - 1;
  T56(1:last) = mkTdim(sol56.y(1,1:last));
  p56(1:last) = mkpdim(sol56.y(2,1:last));
  for i = 1:last
    q56(i) = calcq(T56(i),p56(i));
  end
  z56(1:last) = z6*sol56.x(1:last);
  flow = writeflow(flow,{'z','T','p','a','q','color'},...
			{z56,T56,p56,zeros(1,last+1),q56,'b'});
end

end %------------------------------------------------------- end integrateliquid

%------------------------------------------------------------- integratetwophase
function [T7,q_mh7,hvapK7,pk7,pliq7,z7,flow] = integratetwophase(m,T8,q_mh8,...
				hvapK8,pk8,dpk8,dpcap8,z8,flow,mem,fs,solver)
%INTEGRATETWOPHASE Two-phase flow - copy of FLOW12>FRONT89 and FLOW12>FLOW78.

% The part below is from FLOW12>FRONT89

% Solve
%   (1) xdot9*hvapK9 + q8/m = q_mh8/m,
%   (2) q8 = -k2ph*dT/dz.
% With dp/dz = -m*nu/kappa, dp/dz = dp2ph/dz = dp2ph/dT*dT/dz + dp2ph/da*da/dz.
%  UNDER THE ASSUMPTION  da/dz = 0 (somewhat arbitrarily, i believe)
% follows q/m = k2ph*nu2ph/(kappa*dp2ph/dT)
% Remark: Eq. (1) above sets the enthalpy of the liquid at T8 and pliq8 to zero.
% Eq. (1) follows from m*h + q = const.
sol89 = @(a) q_mh8/m - fs.xdot(T8,pk8,a)*hvapK8 ...
  - fs.k2ph(T8,a)*fs.nu2ph(T8,pk8,a)/((dpk8-(1-a)*dpcap8)*mem.kappa);
a8 = fzero(sol89,[0 1],optimset('TolX',solver.tola));
dp2ph8 = dpk8 - (1-a8)*dpcap8;

% Below follows the part from FLOW12>FLOW78

% Below, hgK returns the enthalpy of the vapor phase at T, pk(T) with respect to
% the enthalpy of the vapor phase at T2, pk(T2). Hence,
% m*h = m(hgk - (1-xdot)hvapK). (Compare to above, m*h = m*xdot*hvapK.)
% Solve
%   (1) m(hgK - (1-xdot)hvapK) + q = m*hgK8 + q_mh8 - m*hvapK8
%   (2) q = -k dT/dz
%   (3) dp/dz = -m nu/kappa = dp2phdT*dT/dz + dp2phda*da/dz.
% HOWEVER, CHECK THIS OUT - TODO!
% Presently, da/dz is set to 0, hence dp2ph/dT = -m nu2ph/kappa, and
% q = -k2ph dT/dz = (-k2ph/(dp2ph/dT)) dp/dz = m nu2ph k2ph/((dp2ph/dT) kappa).

%  CHARACTERISTIC SCALES
% scales to make dimensionless
nu8 = fs.nu2ph(T8,pk8,a8);
pscale = m*z8*nu8/mem.kappa; % zscale = z8;
Tscale = pscale/dp2ph8;
mkTdim = @(Tw) Tscale.*Tw + T8;
coeff1 = -dp2ph8/nu8;
coeff2 = fs.hgK(T8) + q_mh8/m - hvapK8;
step78 = -min(solver.odemaxstep(pscale,solver.maxpperstep), ...
	      solver.odemaxstep(Tscale,solver.maxTperstep));

%  INTEGRATE
% integrate the dimensionless equations; step78/16 is necessary for large
% pressure drops on the order of 2 bar.
options=odeset('InitialSlope',[-1;0],'RelTol',solver.rtol,...
  'Mass',[1 0;0 0],'MStateDependence','none','MassSingular','yes',...
  'Refine',1,'InitialStep',step78/16,'MaxStep',step78); %,'OutputFcn',@odeplot);
sol78 = ode23t(@int78w,[1 0],[0 a8],options);

function dy = int78w(z,y)
  Tw = y(1); a = y(2); T = mkTdim(Tw);
  [pk, dpk, hvK, dpcap] = fs.hvapK(T);
  nu = fs.nu2ph(T,pk,a);
  dp2ph = dpk - (1-a)*dpcap;
  dTw = coeff1*nu/dp2ph;
  da = fs.hgK(T) - (1-fs.xdot(T,pk,a))*hvK ...
       + nu*fs.k2ph(T,a)/(dp2ph*mem.kappa) - coeff2;
  dy = [dTw; da];
end

%  ASSIGN LAST POINT
% dimensionalize last point, here T, a, z
last = size(sol78.x,2);
T7 = mkTdim(sol78.y(1,last));  a7 = sol78.y(2,last);  z7 = z8*sol78.x(last);
if a7 > 1 || a7 < 0
  error(['Two-phase flow vanished, a = %.1f. ' ...
	 'Implement a termination condition.'], a7);
end

[pk7, dpk7, hvapK7, dpcap7, pcap7] = fs.hvapK(T7);
% This is repeated in flsetup.q2ph; Here, though, hvapK is also needed.
q7 = m*fs.nu2ph(T7,pk7,a7)*fs.k2ph(T7,a7)/((dpk7-(1-a7)*dpcap7)*mem.kappa);
q_mh7 = q7 + m*fs.xdot(T7,pk7,a7)*hvapK7;
pliq7 = pk7 - pcap7;

%  WRITE SOLUTION
% if wanted, write the solution
if solver.writesolution
  % allocate space for all points; assign last point
  T78(last) = T7; a78(last) = a7; z78(last) = z7; q78(last) = q7;
  % here the 2ph-pressure, p2ph = pK - (1-a)*pcap, not p2ph = pK
  p78(last) = pk7 - (1-a7)*pcap7;
  last = last - 1;
  T78(1:last) = mkTdim(sol78.y(1,1:last)); a78(1:last) = sol78.y(2,1:last);
  z78(1:last) = z8*sol78.x(1:last);
  % pk, q2ph, is not vektorizable; Gives a result, but probably wrong numbers.
  for i = 1:last
    [q78(i), p78(i)] = fs.q2ph(m,T78(i),a78(i));
  end
  flow = writeflow(flow,{'z','T','p','a','q','color'},{z78,T78,p78,a78,q78,'g'});
end

end %----------------------------------------------------- end integratetwophase

%-------------------------------------------------------------------- ifreevapor
function [T1,p1,q1,zscale,flow] = ifreevapor(m,T3,p3,q3,z3,flow,s,solver)
%IFREEVAPOR Integrate vapor flow in free space - a copy of FLOW12>FLOW13.
% Yields the upstream state. Only the temperature changes.

% m h + q = const => m cp dT + dq = 0,  with dh = cp dT + (dh/dp) dp, dp = 0.
% Then, dT/dq = -1/m*cp. This would be sufficient to solve for T and p, but we
% also want to know a little bit about z:
% dz = -k/q dT = -k/q (-dq/m*cp).
% Scaling:  q = qw*q3; T = T3 + Tw q3/m*cp3. Hence, dq = q3 dqw, dT = q3/m*cp3
% dTw and we have
%   dTw/dqw = -cp3/cp.
% The z-coordinate:  zw = exp(m*cp3/k3 (z-z3)), dz = dzw/(m*cp/k3 zw), therefore
%   dzw/dqw = cp3/cp k/k3 zw/qw.
% Integrate from qw = 1 to qw = 0, initial conditions:
%   Tw(qw=1) = 0, zw(qw=1) = 1.

%  ASSIGN WHAT IS ALREADY KNOWN
p1 = p3;

%  CHARACTERISTIC SCALES
% scales to make dimensionless
cp3 = s.cpg(T3,p3);
Tscale =  q3/(m*cp3); % qscale = q3;

%  SKIP CALCULATION OF THE THERMAL BOUNDARY LAYER
% solver.writesolution is used as a signal that the final solution is computed
% ( Could also do: if abs(...) && ~solver.writesolution - to skip computation
%   only during iteration. )
if abs(Tscale) < solver.rtol
  %  THE END
  q1 = 0;
  T1 = T3;
  zscale = [];
  if solver.writesolution
    warning(['Did not calculate the thermal boundary layer, set q to zero '...
	     'instead of q = %.3g W/m2K'], q3);
  end
  return
end
k3 = s.kg(T3);
zscale =  k3/(m*cp3); % plot z3 + (zw-1)*k3/m*cp3; zw runs from 0 -- 1.
mkTdim = @(Tw) Tscale.*Tw + T3;

%  INTEGRATE
% make initial conditions dimensionless; integrate
% step size depends on the temperature variation
% integration steps
step13 = -solver.odemaxstep(Tscale,solver.maxTperstep);
%T6w = 0; p6w = 0; q6w = 1; z6w = 1;  % O(dTw/dzw) = 1; O(dpw/dzw) = 1;
% Tw(zw=1) = 0, qw(zw=1) = 1; Analytical solution: Tw(zw=0) = 1, dTw/dzw = -1.
options = odeset('RelTol',solver.rtol,'AbsTol',solver.rtol,'Refine',1,...
  'InitialStep',step13,'MaxStep',step13);
sol13 = ode45(@int13q,[1 0],[0 1],options);

function dy = int13q(qw,y)
  Tw = y(1); zw = y(2);
  T = mkTdim(Tw);
  dTw = -cp3/s.cpg(T,p3);
  dzw = -(dTw*s.kg(T)/k3)*(zw/qw);
  if qw == 0
    % in the singularity, zw = 0, qw = 0, zw/qw is replaced by dzw/dqw.
    dzw = sqrt(-(dTw*s.kg(T)/k3)); % should always be positive!
  end
  dy = [dTw; dzw];
end

%  ASSIGN LAST POINT
% dimensionalize last point, here T1
last = size(sol13.x,2);
T1 = mkTdim(sol13.y(1,last));
q1 = 0; % we integrate with q as independent variable, until q = 0.

%  WRITE SOLUTION
% if wanted, write the solution
if solver.writesolution
  T13(last) = T1; p13(1:last) = p1; q13(last) = q1;
  last = last - 1;
  T13(1:last) = mkTdim(sol13.y(1,1:last)); q13(1:last) = q3*sol13.x(1:last);
  % plotted is the characteristic scale, (exp((z-z3)/zscale) - 1)zscale + z3
  z13 = z3 + (sol13.y(2,:)-1)*zscale;
  % To plot the physical z-coordinate:
  %  z3 = z13(1);
  %  (1)  z = z3 + zscale * ln( (z13-z3)/zscale + 1 )
  % Because zw runs from 0 to 1 in flow direction, zw = sol13.y(2,:)),
  % z = ln(zw)*zscale+z3 (see above), zw = exp((z-z3)/zscale), hence we plot
  % z13 = (exp((z-z3)/zscale)-1)*zscale+z3. To recover z, use eq. (1) above.
  flow = writeflow(flow,{'z','T','p','a','q','color'},...
			{z13,T13,p13,ones(1,last+1),q13,'r'});
end

end %------------------------------------------------------------ end ifreevapor

function [T4,p4,q4,z4,flow] = ifreeliquid(m,T5,p5,q5,flow,s,solver)% ifreeliquid
%IFREELIQUID Integrate liquid flow in free space - a copy of FLOW12>FLOW45.
% Only valid for z5 = 0.

% Only the temperature changes, see flow13:
%   dq/dz = - m cp dT/dz, or   dq/dz = m*cp/k q,
%   dT/dz = - q/k.
% dimensionless:
%   T = (T4-T5)*Tw + T5;
%   q = q5*qw;
%   z = (T4-T5)*k5/q5 zw;
% dimensionless Eqs.:
%   dTw/dzw = -k5/k qw;
%   dqw/dzw = -cp m*(T4-T5)/q5 dTw/dzw
% Because m*cp*(T4-T5)/q5, with q5 approx hvap, dqw/dzw is small.
% Initial conditions: At z5 = 0: Tw5 = 0, qw5 = 1, to be integrated
% approximately to zw4 = -1.
% Instead of zw, use T as independent variable and integrate between T5 and T4.
% Advantage: no termination condition necessary and no need to guess a large
% enough range of z.
%   dq/dT = -m cp,   dqw/dTw = -cp m*(T4-T5)/q5
%   dz/dT = -k/q,  dzw/dTw = -k/k5 1/qw.
% TODO, Disadvantage: Only valid for a vapor upstream of the membrane, to
% generalize use termination conditions q = 0 and T < Ts(p).

%  ASSIGN LAST POINT
% we already know part of the solution; needed for characteristic scales.
T4 = s.Ts(p5);

if isempty(T4), warning('Empty Ts in liquid film calculation?'); return; end

p4 = s.ps(T4); % not p4 = p5, instead tolerate a minute pressure glitch

%  CHARACTERISTIC SCALES
% scales to make dimensionless
k5 = s.kl(T5); % cp5 = s.cpl(T5); probably not needed
Tscale = T4 - T5; % qscale = q5; % must be q5 > 0
zscale = Tscale*k5/q5;
mkTdim = @(Tw) Tscale.*Tw + T5;
coeff1 = -m*Tscale/q5;
step45 = -solver.odemaxstep(Tscale,solver.maxTperstep);

%  INTEGRATE
% integrate from Tw5 to Tw4, i.e. from 0 to 1.
options = odeset('RelTol',solver.rtol,'Refine',1,...
		 'InitialStep',step45,'MaxStep',step45);
%sol45 = ode45(@int45,[Tw5 Tw4],[qw5 zw5],options);
sol45 = ode45(@int45,[0 1],[1 0],options);

function dy = int45(Tw,y)
  T = mkTdim(Tw);
  dqw = coeff1*s.cpl(T);
  dzw = -s.kl(T)/(k5*y(1));
  dy = [dqw; dzw];
end

%  ASSIGN LAST POINT
% dimensionalize and assign actual state, here T, p, q, z
% T, p assigned already above
last = size(sol45.x,2);
q4 = q5*sol45.y(1,last); z4 = zscale*sol45.y(2,last);

%  WRITE SOLUTION
% if wanted, write the solution
if solver.writesolution
  T45(last) = T4; p45(1:last) = p4; q45(last) = q4; z45(last) = z4;
% i could do T45(1) = T5; T45(2:last-1) = ...
  last = last - 1;
  T45(1:last) = mkTdim(sol45.x(1:last)); q45(1:last) = q5*sol45.y(1,1:last);
  z45(1:last) = zscale*sol45.y(2,1:last);
  flow = writeflow(flow,{'z','T','p','a','q','color'},...
			{z45,T45,p45,zeros(1,last+1),q45,'b'});
end

end %----------------------------------------------------------- end ifreeliquid

function flow = writeflow(flow,vars,values) %------------------------- writeflow
%WRITEFLOW  Write the solution to the flowstruct.

j = length(flow) + 1;
last = size(vars,2);
for i = 1:last
  flow(j).(vars{i}) = values{i};
end
end %------------------------------------------------------------- end writeflow

function x = newton(fun,x0,res,iter) %----------------------------------- newton
%NEWTON     Newton iteration. Controls tolerance in the solution.
%
%  NEWTON(FUN,X0) finds a solution to F(X) = 0 by Newton iteration, starting
%  from X = X0. NEWTON expects FUN to be the handle to a function that returns a
%  vector [F(X) DF(X)], where DF is the derivative of F at X.
%
%  NEWTON(FUN,X0,RES,ITER) finds a solution to F(X) = 0. NEWTON iterates until
%  X changes less than RES in one step or if more than ITER iterations are done.
%  Default values are RES = 1e-12 and ITER = 100.
%
%  Example
%    Given
%      function [s ds] = senus(x)
%      s = sin(x); ds = cos(x);
%    should give pi, using a function handle to senus.
%      probablypi = newton(@senus,3.1,0)
%
%  See also the MATLAB-functions function_handle, feval.

if nargin < 4, iter = 100; end
if nargin < 3, res = 1e-12; end

x = x0;
for i = 1:iter
  [y, dy] = fun(x);
  dx = y/dy;
  x = x - dx;
  if abs(dx) < res, return; end
end

% Be verbose if no solution is found.
funname = func2str(fun);
error(['NEWTON not succesful. Increase RES or ITER. Type help newton.\n' ...
  '  Function %s, initial guess x0: %g, last found value x: %g,\n' ...
  '  function value %s(x) = %g. Allowed residual: %g, Iterations: %u.'], ...
  funname,x0,x,funname,y,res,iter)
end %---------------------------------------------------------------- end newton
