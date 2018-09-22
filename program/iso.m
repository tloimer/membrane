function [p1,ms] = iso(m,state,ms,solver)
%ISO        Upstream pressure for isothermal flow through a MSTACKSTRUCT.
%  ISO(M,STATE,MSTACKSTRUCT,SOLVER) returns the upstream pressure for the
%  mass flux M of the fluid MSTACKSTRUCT.SUBSTANCE through an assemblage of
%  membranes described by the struct MSTACKSTRUCT. The downstream state of
%  the fluid is given by STATE, the solution tolerances are controlled via
%  the struct SOLVER.
%
%  [P1,MSTACKSTRUCT] = ISO(M,STATE,MSTACKSTRUCT,SOLVER) returns the
%  upstream pressure P1 and a MSTACKSTRUCT describing the solution.
%
%  This is a crude copy of ASYM, hacked to isothermal flow.
%
%  See also ASYM, DOWNSTREAMSTATE, MSTACKSTRUCT, FLOWSETUP, SOLVERSTRUCT.

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

% Next membrane
end
if isempty(state.p)	% equivalent to: state.phase == '2'
  p1 = state.pk;
else
  p1 = state.p;
end

end %%% END ISO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END ISO %%%


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
      interface_liqvap;
    else
      % the logic based on numerically correct comparisons follows
      if p2 < pk1 % mustbecomegaseous
	% nointerface
	state1 = state2;
      elseif p2 > pk1 % mustbecomeliquid
	interface_liqvap;
      else % p2 == pk1, gaseous, liquid or two-phase possible
	% pk1 and pcap1 are calculated again - this does not happen too often, anyway
        warning('Double solution possible. Here, evaporation front.');
	interface_liqvap;
      end
    end

  case 'l' % liquid
    % Deal with the signal p2 == pk1 -pcap1. This may be true only within
    % computer accuracy, and it is raised only within a membrane layer.
    if state2.vapliqequilibrium % p2 == pk1 - pcap1, to within computer accuracy
      interface_vapliq;
    else
      % Signal exception is treated, continue with correct logic
      if p2 > pk1 - pcap1 % mustbecomeliquid
	% nointerface
	state1 = state2;
      elseif p2 < pk1 - pcap1 % mustbecomegaseous
	interface_vapliq;
      else % gaseous, liquid or two-phase possible
	% pk1 and pcap1 are calculated again - this does not happen too often, anyway
	    warning('Double solution possible. Here, condensation front.');
	  interface_vapliq;
      end
    end

  otherwise
    error('How can we come here?');
end

%--- nested functions ----------------------------------------- nested functions

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
  p1 = p2 - pcap;
  q1 = q2;
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
  % Pade approximation, ln(p1/psat) ~ 2*(p1-psat)/(p1+psat),
  % something wrong, did not work
  %fprintf('Initial guess, Taylor: %.4f, Pade: %.4f, ', p1,...
  %  (-rhoRT-(psat-p2)*.5 + sqrt((rhoRT + (psat-p2)*.5)^2 + psat*(2.*rhoRT-p2)))/psat);
  % setup newton
  ps_rhoRT = psat/rhoRT;
  function [F, dF] = sol35(x)
    F = log(x) + x*ps_rhoRT - p2/rhoRT;
    dF = 1/x + ps_rhoRT;
  end
  p1 = newton(@sol35,p1,1e-15) * psat;
  %fprintf('Result: %.4f\n', p1/psat);
  q1 = q2;
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
deltazw = pscale/pdarcy; % = 1 if pscale = pdarcy
step92 = -min(solver.odemaxstep(pscale,solver.maxpperstep),solver.maxpperstep/pdarcy);
% functions to re-calculate dimensional values
mkpdim = @(pw) pscale.*pw + p2;
function [p, z] = mkdimensional(pw, zw)
  p = mkpdim(pw);
  z = mem.L.*zw;
end
% coefficients for the eqs. in int92w:
coeff1 = -m*mem.L/(mem.kappa*pscale);
pk2 = fs.pkelv(T2);

%  INTEGRATE
% make initial conditions dimensionless; integrate
options=odeset('Events',@term92w,...%'AbsTol',atol*[1 1 -coeff1*cpg2],...
  'RelTol',solver.rtol,'Refine',1,'InitialStep',step92,'MaxStep',step92);
sol92 = ode45(@int92w, [1 0], 0, options);

function dy = int92w(z,y)
  p = mkpdim(y);
  dy = fs.nuapp(T2, p)*coeff1;
end

function [val,isterm,direction] = term92w(z,y)
  % Terminate integration when pk(T)-p = 0, only when falling. OK: pk(T)-p > 0.
  % without direction, integration already terminates at start
  isterm = 1; direction = 0; %-1;
  % T = y(1), p = y(2)
  % dimensional: val=pk(y(1))-y(2);
  % dimensionless:
  val = pk2 - mkpdim(y);
end

%  ASSIGN LAST POINT
% dimensionalize and assign actual state, here T9, p9, q9, z9
last = size(sol92.x,2);
T9 = T2;
q9 = q2;
[p9, z9] = mkdimensional(sol92.y(last), sol92.x(last));
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
  T92(1:last) = T2;
  q92(1:last) = q2;
  % allocate space for all points; assign last point
  p92(last) = p9; z92(last) = z9;
  Kn92(last) = fs.knudsen(T9, p9);
  % and write the solution but the last point
  last = last - 1;
  [p92(1:last), z92(1:last)] = mkdimensional(sol92.y(1:last), sol92.x(1:last));
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
% integration steps
% step56 = -min(odemaxstep(Tscale,maxTperstep),odemaxstep(pscale,maxpperstep))
step56 = -solver.odemaxstep(pscale,solver.maxpperstep);
% functions to re-calculate dimensional values
mkpdim = @(pw) pscale.*pw + p6;
[pk6, pcap6] = fs.pkpcap(T6);

%  INTEGRATE
% integrate; dimensionless initial conditions
%T6w = 0; p6w = 0; q6w = 1; z6w = 1;  % O(dTw/dzw) = 1; O(dpw/dzw) = 1;
options=odeset('Events',@term56w,'RelTol',solver.rtol,...
  'Refine',1,'InitialStep',step56,'MaxStep',step56);
%sol56 = ode45(@int56w,[z6w 0],[T6w p6w],options);
sol56 = ode45(@int56w, [1 0], 0, options);

function dy = int56w(z,y)
  dy = -s.nul(T6)/nul6;
end

function [val,isterm,direction] = term56w(z,y)
  % Termination condition: p = pk(T) - pcap, valid: p > pk(T) - pcap, pcap =
  % curv*sigma.
  % Terminate integration when p + pcap - pk = 0, falling.
  isterm = 1; direction = -1;
  val = mkpdim(y) + pcap6 - pk6;
end

%  ASSIGN LAST POINT
% dimensionalize and assign actual state, here T5, p5, q5, z5
last = size(sol56.x,2);
T5 = T6;
q5 = q6;
p5 = mkpdim(sol56.y(last));
z5 = sol56.x(last)*z6;

% RAISE A SIGNAL IF INTEGRATION TERMINATED PREMATURELY, see integratevapor
if ~isempty(sol56.ie) % integration terminated by event function term56w
  vapliqequilibrium = true;
end

%  WRITE SOLUTION
% if wanted, write the solution
if solver.writesolution
  T56(1:last) = T6;
  q56(1:last) = q6;
  % allocate space for all points; assign last point
  p56(last) = p5;  z56(last) = z5;
  % and write remaining solution
  last = last - 1;
  p56(1:last) = mkpdim(sol56.y(1:last));
  z56(1:last) = z6*sol56.x(1:last);
  flow = writeflow(flow,{'z','T','p','a','q','Kn','color'},...
			{z56,T56,p56,zeros(1,last+1),q56,zeros(1,last+1),'b'});
end

end %------------------------------------------------------- end integrateliquid

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
