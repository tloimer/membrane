function fl = flow12(m,T2,p2,fl,flsetup,solver)
%FLOW12     Flow through the membrane.
%  FLOW12(M,T2,P2,FLOWSTRUCT,FLSETUP,SOLVER) returns a structure FLOWSTRUCT that
%  contains the solution for the mass flux M [kg/m2s] at a downstream
%  temperature T2 [K] and a downstream pressure P2 [Pa]. FLSETUP is a struct
%  that contains input to FLOW12, calculated from combining SUSTANCE, MEMBRANE
%  and FMODEL properties and the maximally expected Temperature. The struct
%  FLOWSTRUCT is also used to pass input to FLOW12. Solution tolerances are
%  controlled via the struct SOLVER.
%
%  FLSETUP is used because if MNUM iteratively calls FLOW12 to find a solution,
%  the setup only needs to be calculated once.
%
%  Fields in FLOWSTRUCT (FL) used for input to FLOW12:
%    FL.info.theta      Contact angle (degree).
%    FL.info.substance  Struct SUBSTANCE.
%    FL.info.membrane   Struct MEMBRANE.
%    FL.info.fmodel     Struct FMODEL.
%
%  Fields in FLOWSTRUCT used for output from FLOW12:
%    FL.sol.m           Mass flux [kg/m2s].
%    FL.sol.T1          Upstream temperature [K].
%    FL.sol.p1          Upstream pressure [Pa].
%    FL.sol.q1          Upstream heat flux [W/m2].
%    FL.sol.T3          Temperature at the liquid film or the membrane front.
%    FL.sol.states      Flow states, see below.
%    FL.sol.len         Length of the struct FL.flow. Negative number.
%    FL.flow            Struct of length -FL.sol.len.
%  If SOLVER.writesolution is true, the following fields are written: 
%    FL.flow.z          Coordinate [m].
%    FL.flow.T          Temperature [K].
%    FL.flow.p          Pressure [p].
%    FL.flow.q          Heat flux [W/m2].
%    FL.flow.Kn         Knudsen number, for vapor flow.
%    FL.flow.a          Vapor volume fraction.
%    FL.flow.color      Plotting specification.
%
%  FLSETUP must contain the fields
%     .curv             Curvature of the meniscus in a pore.
%     .kelv             Function to calculate the ratio pkelv/psat.
%     .pkelv            Function to calculate the pressure pkelv.
%     .hgK              Function to return the vapor enthalpy hvap(T,pkelv(T)).
%     .intdhdpdpsatdT
%     .nuapp
%     .knudsen
%     .nu2ph
%     .kmgas
%     .kmliq
%     .k2ph
%     .xdot
%     .odemaxstep        Function to calculate integration step width.
%
%  SOLVER must contain the fields
%    .rtol              Relative error tolerance in odeset ('RelTol').
%    .atol              Absolute error tolerance in odeset ('AbsTol').
%    .tola              Tolerance for FZERO. Default in FZERO: 1e-16.
%    .maxTperstep       Maximum temperature difference per integration step.
%    .maxpperstep       Maximum pressure difference per integration step.
%    .writesolution     If true, write FL.flow.
%    .partialsolution   Only write a partial solution: p1, T3.
%
%  See also FMODEL, MEMBRANE, MNUM, SUBSTANCE.
%
% Nested functions: hvapK, hvapKraw, from2, flow92, front62, flow56,
%    check69or89, front89, flow78, front37, front57, flow13, flow45, front35,
%    front85, writetostruct.
%  Subfunction: newton.
%  Try, e.g., help flow12>newton.
%
%  States of the fluid (see Table 2 in Loimer, 2007):
%    1 - upstream state
%    2 - downstream state
%    3 - unsaturated vapor at upstream membrane front
%    4 - saturated liquid at upstream end of liquid film
%    5 - liquid in the liquid film at the upstream membrane front and/or
%    5 - liquid at upstream end of liquid flow within the membrane
%    6 - liquid at downstream end of liquid flow within the membrane
%    7 - two-phase mixture at upstream end of 2ph-flow within the membrane
%    8 - two-phase mixture at downstream end of 2ph-flow within the membrane
%    9 - vapor at upstream end of vapor flow within the membrane
%
%  Note: State 5 may denote two different states in the non-linear
%  calculation, at contrast to the linear approximation. For conformity, an
%  additional name is not introduced.

%  TODO
% rewrite intcpl(T1,T2) to icpl(T), See substance. Would have to implement
% (rewrite) icpleq2, ipoly4, ipoly3.
% probably use q2ph, to calculate two-phase heat flux

%  COMMON VARIABLES
s = fl.info.substance;  mem = fl.info.membrane;  f = fl.info.fmodel;
% definition of colors, i.e., LineSpecs
liqcolor = 'b';  vapcolor = 'r';  twophcolor = 'g-+';

% Combination of theta and membrane.
% fl.info.theta is accessed directly
curv = flsetup.curv; %curv = fl.info.curv;

%  COMMON FUNCTIONS
% Thermodynamic properties from combination of substance, theta and membrane.
kelv = flsetup.kelv;  pkelv = flsetup.pkelv;
hgK = flsetup.hgK;  intdhdpdpsatdT = flsetup.intdhdpdpsatdT;

% Effective single flow - membrane and two-phase properties.
nuapp = flsetup.nuapp;  knudsen = flsetup.knudsen;  nu2ph = flsetup.nu2ph;
kmgas = flsetup.kmgas;  kmliq = flsetup.kmliq;  k2ph = flsetup.k2ph;
xdot = flsetup.xdot;  %x = flsetup.x;

% Step size calculation for ode45, ode23t.
odemaxstep = flsetup.odemaxstep;

% Solver values.
[rtol atol tola maxTperstep maxpperstep writesolution partialsolution] ...
  = deal(solver.rtol, solver.atol, solver.tola, solver.maxTperstep,...
  solver.maxpperstep, solver.writesolution, solver.partialsolution);


%  INITIALIZATION
fl.sol.states = [];
if writesolution
  % initialize writestruct
  fl.sol.len = 0;
  % make sure, fl.flow is empty
  if ~isempty(fl.flow)
    fl.flow = fl.flow(1);
    [fl.flow.z fl.flow.T fl.flow.p fl.flow.q fl.flow.a fl.flow.color] ...
      = deal([]);
    % in 92, we assign a 'Kn'-field.
    if isfield(fl.flow,'Kn'), fl.flow = rmfield(fl.flow,'Kn'); end
  end
end

%  START
from2(m,T2,p2);
% The rest plays all in the nested functions.

%%% NESTED FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NESTED FUNCTIONS %%%

function [pk dpk hvapK] = hvapK(T) %-------------------------------------- hvapK
%HVAPK      Enthalpy of vaporization at a curved interface.
%
% [PK DPK HVK] = HVAPK(T) returns the enthalpy of vaporization as well as the
% vapor pressure and the derivative of the vapor pressure with respect to
% temperature at a curved interface. [PK] = Pa, [DPK] = Pa/K, [HVK] = J/kg.
%
% DPK: See Eq. (11) in Loimer (2007).
%
% HVK: See Eq. (14) in [Loimer, Proc. STAMM 2004; Wang, Hutter (eds.), Trends in
% Applications of Mathematics to Mechanics, 247-256, Shaker (2005)].
% With dh/dp = v - T dv/dT.
[psat dps] = s.ps(T);
[dsig sigma] = s.dsig(T);
[drho rho] = s.drho(T);
pk_ps = kelv(T,sigma,rho);
pk = pk_ps*psat; pcap = curv*sigma;
dpk = pk_ps * (dps + psat*pcap*(1/T-dsig+drho)/(s.R*rho*T));
%hvapK = s.hvap(T) + (pk-psat)*(s.dhdp(T,pk)-dhldp) + dhldp*curv*sigma;
hvapK = hvapKraw(T,pk,psat,pcap,rho,drho);
end %----------------------------------------------------------------- end hvapK

function hvK = hvapKraw(T,prad,psat,pcap,rho,drho) %------------------- hvapKraw
  dhldp = (1 + T*drho)/rho;
  hvK = s.hvap(T) + (prad-psat)*(s.dhdp(T,prad)-dhldp) + dhldp*pcap;
end %-------------------------------------------------------------- end hvapKraw

function from2(m,T2,p2) %------------------------------------------------- from2
%FROM2      Start flow through the membrane.
% p2 >= pk(T2)  Are we in the capillary condensation regime?
%  |- no:  flow92 (unsaturated vapor flow)
%  |- yes: front62 -> flow56 (condensation and liquid flow)

% FROM2 used in order to encapsulate variables (T6,p6,q6).
q2 = 0;
if p2 >= pkelv(T2)
  [T6 p6 q6] = front62(m,T2,p2,q2);
  flow56(m,T6,p6,q6,mem.L);
else
  flow92(m,T2,p2,q2); % In flow92: z2 = mem.L.
end
end %----------------------------------------------------------------- end from2

function flow92(m,T2,p2,q2) %-------------------------------------------- flow92
%FLOW92     Vapor flow within the membrane.
% Termination condition: p = pk(T), valid: p < pk(T).
% If FLOW92 terminates before z = 0
%  |- yes: full condensation possible ?
%  |   |- yes: front69 -> flow56
%  |   |- no:  front89 -> flow78
%  |- no:  p9 > psat(T9)  Condensation in front of the membrane possible?
%  |   |- yes: front59 -> flow45  (p9 > psat only for non-wetting possible)
%  |   |- no:  ccond39 -> flow13  (unsaturated vapor in front of membrane;
%  |   |                           limit: saturated vapor at z = 0 possible)

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
pdarcy = m*nuapp(T2,p2)*mem.L/mem.kappa;
pscale = min(s.ps(T2) - p2,pdarcy);
%[dhdp2 cpg2] = s.dhcpg(T2,p2); % s.dhdp
dhdp2 = s.dhdp(T2,p2);
qscale = -m*dhdp2*pscale;
deltazw = pscale/pdarcy; % = 1 if pscale = pdarcy
Tscale = qscale*mem.L*deltazw/mem.km;
step92 = -min(odemaxstep(pscale,maxpperstep),maxpperstep/pdarcy);
% functions to re-calculate dimensional values
mkTdim = @(Tw) Tscale.*Tw + T2;
mkpdim = @(pw) pscale.*pw + p2;
function [T p q z] = mkdimensional(Tw,pw,qw,zw)
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
  'RelTol',rtol,'Refine',1,'InitialStep',step92,'MaxStep',step92);
%sol92 = ode45(@int92w,[z2w 0],[T2w p2w q2w],options);
sol92 = ode45(@int92w,[1 0],[0 0 q2w],options);

%function dy = int92(z,y)
%% dimensional. Equations see above.
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
  dpw = nuapp(T,p)*coeff1;
  dTw = coeff2*qw/kmgas(T);
  [dhdp cpg] = s.dhcpg(T,p);
  dqw = coeff3*cpg*dTw + dpw*dhdp/dhdp2;
  dy = [dTw;dpw;dqw];
end

function [val,isterm,direction] = term92w(z,y)
  % Terminate integration when pk(T)-p = 0, only when falling. OK: pk(T)-p > 0.
  % without direction, integration already terminates at start
  isterm = 1; direction = 0; %-1;
  % T = y(1), p = y(2)
  % dimensional: val=pk(y(1))-y(2);
  % dimensionless:
  val = pkelv(mkTdim(y(1))) - mkpdim(y(2));
end

%  ASSIGN LAST POINT
% dimensionalize and assign actual state, here T9, p9, q9, z9
last = size(sol92.x,2);
[T9 p9 q9 z9] = mkdimensional(sol92.y(1,last),sol92.y(2,last),...
   sol92.y(3,last),sol92.x(last));
% not necessary: .ye, .xe are also in .y(last), .x(last).
%  [T9 p9 q9 z9] = mkdimensional(sol92.ye(1,end),sol92.ye(2,end),...
%    sol92.ye(3,end),sol92.xe(end));

%  WRITE SOLUTION
% if wanted, write the solution
if writesolution
  % allocate space for all points; assign last point
  T92(last) = T9; p92(last) = p9; q92(last) = q9; z92(last) = z9;
  Kn92(last) = knudsen(T9,p9);
  pk92(last) = pkelv(T9); %DEBUG
  % and write the solution but the last point
  last = last - 1;
  [T92(1:last) p92(1:last) q92(1:last) z92(1:last)] = mkdimensional(...
    sol92.y(1,1:last),sol92.y(2,1:last),sol92.y(3,1:last),sol92.x(1:last));
  for i = 1:last
    Kn92(i) = knudsen(T92(i),p92(i));
    pk92(i) = pkelv(T92(i));
  end
%DEBUG  writetostruct('92',{'T','p','q','Kn','z','color'},...
%DEBUG    {T92,p92,q92,Kn92,z92,vapcolor});
  writetostruct('92',{'T','p','q','Kn','z','color','pk'},... %DEBUG
    {T92,p92,q92,Kn92,z92,vapcolor,pk92}); %DEBUG
end

%  DECIDE AND CALL NEXT
if not(isempty(sol92.ie))
  % terminated before z = 0
  % full condensation possible?
  % p9 = pk(T9), because this is the termination condition
  [condensation pk9 dpkdT9 hvapK9 q6] = check69or89(m,T9,q9);% q6 for front69
  if condensation
    % yes: front69 -> flow56
    T6 = T9; p6 = p9 - curv*s.sigma(T9); %============================== front69
    if writesolution, writetostruct('-',{},{}); end
    flow56(m,T6,p6,q6,z9);
  else
    % no:  front89 -> flow78
    [T8 a8 doth8] = front89(m,T9,p9,q9,hvapK9,dpkdT9);
    flow78(m,T8,a8,z9,doth8,p9,dpkdT9,hvapK9);
  end
else % not(isempty(sol92.ie))
  % vapor flow all through to z = 0
  % liq. film or vapor flow?
  if p9 > s.ps(T9)
    % yes: front59 -> flow45  (p9 >= psat only for non-wetting possible)
    % front59 is the same as front62
    [T5 p5 q5] = front62(m,T9,p9,q9); %================================= front59
    flow45(m,T5,p5,q5);
  else
    % no:  ccond39 -> flow13  (unsaturated vapor in front of membrane)
    T3 = T9; p3 = p9; q3 = q9; %======================================== ccond39
    flow13(m,T3,p3,q3,0);
  end
end % end terminated
end %---------------------------------------------------------------- end flow92

function [T6 p6 q6] = front62(m,T2,p2,q2) %----------------------------- front62
%FRONT62    Evaporation of liquid at a meniscus at downstream end of membrane.
% Kelvin's eq. + Young-Laplace eq.; radius of curvature must be found.
% ln(p2/psat(T2)) = -pcap/(s.rho*R*T); pcap = p2 - p6;
% -pcap = ln(p2/psat)*rho*R*T;
psat2 = s.ps(T2); [drho2 rho2] = s.drho(T2);
pcap = log(psat2/p2)*rho2*s.R*T2;
hvapK2 = hvapKraw(T2,p2,psat2,pcap,rho2,drho2);
T6 = T2;
p6 = p2 - pcap;
q6 = m*hvapK2 + q2;
end %--------------------------------------------------------------- end front62

function flow56(m,T6,p6,q6,z6) %----------------------------------------- flow56
%FLOW56     Liquid flow within the membrane.
% Termination condition: p = pk(T) - pcap, valid: p > pk(T) - pcap,
%   pcap =  curv*sigma.
% If FLOW56 goes through to z = 0
%  |- yes: p5 > psat(T5)
%  |   |- yes: ccond55 -> flow45  (liquid film in front of membrane)
%  |   |- no:  front35 -> flow13  (unsaturated or saturated vapor)
%  |- no: front85 -> flow78  (z > 0, partial condensation within the membrane)

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
nul6 = s.nul(T6);  kmliq6 = kmliq(T6);
pscale = z6*m*nul6/mem.kappa;
Tscale = q6*z6/kmliq6;
% zscale = z6; qscale = q6;
% integration steps
% liquid properties only depend on temperature, not on pressure
% step56 = -min(odemaxstep(Tscale,maxTperstep),odemaxstep(pscale,maxpperstep))
step56 = -odemaxstep(Tscale,maxTperstep);
% functions to re-calculate dimensional values
mkTdim = @(Tw) Tscale.*Tw + T6;
mkpdim = @(pw) pscale.*pw + p6;
% coefficients for the eqs. in int56w:
% primarily, for the invariable part of (q6 + m*h6)/q6:
[drho6 rho6] = s.drho(T6);
coeff1 = q6 + m*( (p6-s.ps(T6))*(1+drho6)/rho6 - intdhdpdpsatdT(T6) );
coeff2 = -kmliq6/q6;

% calculate q. This function is needed twice, once during integration and later
% for evaluation of the solution.
function q = calcq(T,p)
  % do not vectorize! s.ps is scalar, gives wrong values on vectors.
  [drho rho] = s.drho(T);
  q = coeff1 - m ...
    * (s.intcpl(T6,T) + (p-s.ps(T))*(1+drho)/rho - intdhdpdpsatdT(T));
end

%  INTEGRATE
% integrate; dimensionless initial conditions
%T6w = 0; p6w = 0; q6w = 1; z6w = 1;  % O(dTw/dzw) = 1; O(dpw/dzw) = 1;
options=odeset('Events',@term56w,'RelTol',rtol,...
  'Refine',1,'InitialStep',step56,'MaxStep',step56);
%sol56 = ode45(@int56w,[z6w 0],[T6w p6w],options);
sol56 = ode45(@int56w,[1 0],[0 0],options);

function dy = int56w(z,y)
  T = mkTdim(y(1));  p = mkpdim(y(2));
  dTw = calcq(T,p)*coeff2/kmliq(T);
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
  val = mkpdim(y(2)) + curv*s.sigma(T) - pkelv(T);
%  disp(sprintf('T = %g, p = %g: val = %g',T-273.15,mkpdim(y(2)),val));%DEBUG
end

%  ASSIGN LAST POINT
% dimensionalize and assign actual state, here T5, p5, q5, z5
last = size(sol56.x,2);
T5 = mkTdim(sol56.y(1,last));  p5 = mkpdim(sol56.y(2,last));
q5 = calcq(T5,p5); z5 = sol56.x(last)*z6; % z5 = 0;

%  WRITE SOLUTION
% if wanted, write the solution
if writesolution
  % allocate space for all points; assign last point
  T56(last) = T5; p56(last) = p5; q56(last) = q5; z56(last) = z5;
  pp56(last) = p5 + curv*s.sigma(T5); pk56(last) = pkelv(T5); %DEBUG
  % and write remaining solution
  last = last - 1;
  T56(1:last) = mkTdim(sol56.y(1,1:last));
  p56(1:last) = mkpdim(sol56.y(2,1:last));
  for i = 1:last
    q56(i) = calcq(T56(i),p56(i));
    pp56(i) = p56(i) + curv*s.sigma(T56(i)); %DEBUG
    pk56(i) = pkelv(T56(i)); %DEBUG
  end
  z56(1:last) = z6*sol56.x(1:last);
%DEBUG  writetostruct('56-',{'T','p','q','z','color'},{T56,p56,q56,z56,liqcolor});
  writetostruct('56',{'T','p','q','z','color','pp','pk'},... %DEBUG
  {T56,p56,q56,z56,liqcolor,pp56,pk56}); %DEBUG
end

%  DECIDE AND CALL NEXT
if isempty(sol56.ie)
  % got through to z = 0
  if p5 > s.ps(T5)
    % liquid film in front of membrane
    % ccond55 %========================================================= ccond55
    flow45(m,T5,p5,q5);
  else
    % saturated or unsaturated vapor in front of membrane
    [T3 p3 q3] = front35(m,T5,p5,q5);
    flow13(m,T3,p3,q3,0);
  end
else % ~isempty(sol56.ie), terminated before z = 0, at z > 0
  [T8 a8 doth8 pk8 dpkdT8 hvapK8] = front85(m,T5,p5,q5);
  flow78(m,T8,a8,z5,doth8,pk8,dpkdT8,hvapK8);
end % end gottrough
end %---------------------------------------------------------------- end flow56

function [condensation pk9 dpkdT9 hvapK9 q6] = check69or89(m,T9,q9)% check69or89
%CHECK69OR89 Evaporation of liquid (69) or two-phase mixture (89)?
% Check if full condensation is possible by taking the 2ph-flow to the limit of
% full condensation,  a -> 0, k -> kmliq, nu -> nuliq.
%   (liq)  m*h6 + q6 = m*h9 + q9,
%   (2ph)  m*xdot*hvapK + q8 = m*hvapK + q9,   if h = 0 for sat. liquid at T8.
% In 2ph-flow the heat flux is a function of a,
%   q8 = -k*dT/dz = -k*(dT/dp)*dp/dz = k*nu*m/((dpk/dT)*kap)
% For full condensation q6 must be small enough such that the temperature falls
% below the saturation temperature. q = q8,liq gives the temperature gradient
% such that T follows Tk(p), hence q6 <=! q8,liq, with q6 = m*hvapK + q9 and q8
% from above
% full condensation possible:
%   m*kmliq*nuliq/((dpk/dT)*kap) >= m*hvapK + q9.
% However, since min(q8) not= q8,liq in the case of the homogeneous 2ph-flow
% model (but close to it, see Fig. 6 in Loimer, 2007), eq. (2ph) may still have
% a solution,
%   q8 = (1-xdot) mhvapK + q9, although ~min(q8) >= max(rhs).
% The additional solutions are not tracked, instead full condensation is used.
% Full condensation is necessary if eq. (2ph) does not have a solution, which my
% be estimated with min(q8) >= m*hvapk + q9.

% three vars we calculate once (and pass to front89 and flow78)
[pk9 dpkdT9 hvapK9] = hvapK(T9);
q6 = q9 + m*hvapK9; %========================================== part for front69
if m*kmliq(T9)*s.nul(T9)/(dpkdT9*mem.kappa) >= q6
  condensation = true;
else
  condensation = false;
end
end %----------------------------------------------------------- end check69or89

function [T8 a8 doth8] = front89(m,T9,p9,q9,hvapK9,dpkdT9) %------------ front89
%FRONT89    Full evaporation of two-phase mixture within the membrane.
% Find a8 that solves 
%   (1) m*xdot*hvapK + q8 = m*hvapK + q9,   if h = 0 for sat. liquid at T = T8.
% with q8 = q8(a8) = -k*dT/dz = -k*(dT/dp)*dp/dz = k*nu*m/((dpk/dT)*kap)
% Eq. (1) may have two solutions, cf. Fig. 6 in (Loimer, 2007) for homogeneous
% flow. However, we only come here for q8,liq < m*hvapK + q9 where only one
% solution is possible, cf. Fig. 6.
T8 = T9; doth8 = q9 + m*hvapK9; %pk8 = p9;
% (1-xdot)*hvapK + q9/m - k*nu/((dpk/dT)*kap) =! 0;
sol89 = @(a) (1-xdot(T9,p9,a))*hvapK9 + q9/m ...
  - k2ph(T9,a)*nu2ph(T9,p9,a)/(dpkdT9*mem.kappa);
a8 = fzero(sol89,[0 1],optimset('TolX',tola));
end %--------------------------------------------------------------- end front89

function flow78(m,T8,a8,z8,doth8,pk8,dpkdT8,hvapK8) %-------------------- flow78
%FLOW78     Two-phase flow within the membrane.
% Does not terminate, goes through to z = 0.
% p7 > psat(T7) ?
%  |- yes: front57 -> flow45  (non-wetting, liq. film in front of membrane)
%  |- no:  front37 -> flow13  (saturated or unsaturated vapor in front)

% integrate two-phase flow (differential-algebraic eqs.)
% eqs.: m*h + q = const; m*h = m(hgk - (1-xdot)hvapK); hgK ist the enthalpy of
% the vapor phase at pk, hg(T,pk(T));
% with q = -k dT/dz = (-k/(dpk/dT)) dp/dz = m nu k /((dpk/dT) kappa)
% abbreviations: nu, k = nu2ph, k2ph
% dT/dz = -m nu/(dpk/dT kappa)
% m(hgK-(1-xdot)hvapK) + m nu k/((dpk/dT) kappa)
%  = m(hgK8-(1-xdot8)hvapK8) + m nu8 k8/((dpkdT8 kappa)
% ( = m*hgK8 + doth8 - m*hvapK8,  see front89: doth8 = m*xdot*hvapK8 + q8 ).
% dimensionless
% z = z8*zw; T = (z8 m nu8/(dpkdT8 kappa)*Tw + T8;
% dimensionless eqs:
% dTw/dzw = -(nu/(dpk/dT) (dpk/dT)8/nu8;  nu8 = nu2ph(T8,a8);

%  CHARACTERISTIC SCALES
% scales to make dimensionless
nu8 = nu2ph(T8,pk8,a8); % tmp1 = nu8/(dpkdT8*mem.kappa);
pscale = m*z8*nu8/mem.kappa; % zscale = z8;
Tscale = pscale/dpkdT8;
mkTdim = @(Tw) Tscale.*Tw + T8;
coeff1 = -dpkdT8/nu8;
%coeff2 = hgK(T8) - (1-xdot(T8,pk8,a8))*hvapK8 + k2ph(T8,a8)*tmp1;
coeff2 = hgK(T8) + doth8/m - hvapK8; %+ q9/m;
step78 = -min(odemaxstep(pscale,maxpperstep),odemaxstep(Tscale,maxTperstep));

%  INTEGRATE
% integrate the dimensionless equations; step78/16 is necessary for large
% pressure drops on the order of 2 bar.
options=odeset('InitialSlope',[-1;0],'RelTol',rtol,...
  'Mass',[1 0;0 0],'MStateDependence','none','MassSingular','yes',...
  'Refine',1,'InitialStep',step78/16,'MaxStep',step78);
sol78 = ode23t(@int78w,[1 0],[0 a8],options);

function dy = int78w(z,y)
  Tw = y(1); a = y(2); T = mkTdim(Tw);
  [p dpk hvK] = hvapK(T);
  nu = nu2ph(T,p,a);
  dTw = coeff1*nu/dpk;
  da = hgK(T) - (1-xdot(T,p,a))*hvK + nu*k2ph(T,a)/(dpk*mem.kappa) - coeff2;
  dy = [dTw;da];
end

%  ASSIGN LAST POINT
% dimensionalize last point, here T7, a7
last = size(sol78.x,2);
T7 = mkTdim(sol78.y(1,last));  a7 = sol78.y(2,last);  z7 = z8*sol78.x(last);
% z7 = 0;  % p7 = pk(T7);

%  WRITE SOLUTION
% if wanted, write the solution
if writesolution
  % allocate space for all points; assign last point
  T78(last) = T7; a78(last) = a7; z78(last) = z7; p78(last) = pkelv(T7);
  last = last - 1;
  T78(1:last) = mkTdim(sol78.y(1,1:last)); a78(1:last) = sol78.y(2,1:last);
  z78(1:last) = z8*sol78.x(1:last);
  for i = 1:last
    % pkelv is not vektorizable; Gives a result, but probably wrong numbers.
    p78(i) = pkelv(T78(i));
  end
  % vielleicht q78 berechnen? q2ph?
  writetostruct('78-',{'T','p','a','z','color'},{T78,p78,a78,z78,twophcolor});
end

%  DECIDE AND CALL NEXT
if fl.info.theta > 90 % equivalent: theta > 90 or p7 > s.ps(T7), costheta < 0
  % yes: front57 -> flow45 (non-wetting, liq. flow in front of membrane)
  [T5 p5 q5] = front57(m,T7,a7);
  if partialsolution
    % THE END
    fl.sol.p1 = p5;
  else
    flow45(m,T5,p5,q5);
  end
else % theta > 90 % curv <=0, p7 <=s.ps(T7)
  % no:  front37 -> flow13  (saturated or unsaturated vapor in front)
  [T3 p3 q3] = front37(m,T7,a7);
  if partialsolution
    % THE END
    fl.sol.p1 = p3;
  else
    flow13(m,T3,p3,q3,0);
  end
end
end %---------------------------------------------------------------- end flow78

function [T3 p3 q3] = front37(m,T7,a7) %-------------------------------- front37
%FRONT37    Partial condensation at z = 0.
% With mh + q = const,
%   m hvapK + q3 = m xdot7 hvapK + q7,
%   q3 = q7 + m*(xdot-1)*hvapK, for q7 see front89
T3 = T7;
%p3 = pkelv(T7); % p3 = p7;
[p3 dpkdT7 hvapK7] = hvapK(T7);
q7 = m*nu2ph(T7,p3,a7)*k2ph(T7,a7)/(dpkdT7*mem.kappa); %q2ph?
q3 = q7 + m*(xdot(T7,p3,a7)-1)*hvapK7;
end  %-------------------------------------------------------------- end front37

function [T5 p5 q5] = front57(m,T7,a7) %-------------------------------- front57
%FRONT57    Partial evaporation at z = 0.
% Liquid film in front, occurs for non-wetting.
% With mh + q = const,
%   q5 = m xdot7 hvapK + q7,
T5 = T7;
%p7 = pk(T7);
[p7 dpkdT7 hvapK7] = hvapK(T7);
p5 = p7 - curv*s.sigma(T7);
q7 = m*nu2ph(T7,p7,a7)*k2ph(T7,a7)/(dpkdT7*mem.kappa); %q2ph?
q5 = q7 + m*xdot(T7,p7,a7)*hvapK7;
end  %-------------------------------------------------------------- end front57

function flow13(m,T3,p3,q3,z3) %----------------------------------------- flow13
%FLOW13     Vapor flow in front of the membrane.
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

%  CHARACTERISTIC SCALES
% scales to make dimensionless
cp3 = s.cpg(T3,p3);
Tscale =  q3/(m*cp3); % qscale = q3;

% DONE ALREADY?
if abs(Tscale) < rtol 
  %  THE END
  % repeated below
  fl.sol.p1 = p3;  fl.sol.T1 = T3;  fl.sol.T3 = T3; fl.sol.q1 = q3;
  return
end
k3 = s.kg(T3);
zscale =  k3/(m*cp3); % plot z3 + (zw-1)*k3/m*cp3; zw runs from 0 -- 1.
mkTdim = @(Tw) Tscale.*Tw + T3;

%  INTEGRATE
% make initial conditions dimensionless; integrate
% step size depends on the temperature variation
% integration steps
step13 = -odemaxstep(Tscale,maxTperstep);
%T6w = 0; p6w = 0; q6w = 1; z6w = 1;  % O(dTw/dzw) = 1; O(dpw/dzw) = 1;
% Tw(zw=1) = 0, qw(zw=1) = 1; Analytical solution: Tw(zw=0) = 1, dTw/dzw = -1.
options = odeset('RelTol',rtol,'AbsTol',rtol,'Refine',1,'InitialStep',step13,...
  'MaxStep',step13);
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
p1 = p3;
% dimensionalize last point, here T1
last = size(sol13.x,2);
T1 = mkTdim(sol13.y(1,last));
q1 = 0; % we integrate with q as independent variable, until q = 0.

%  WRITE SOLUTION
% if wanted, write the solution
if writesolution
  T13(last) = T1; p13(1:last) = p3; q13(last) = q1;
  last = last - 1;
  T13(1:last) = mkTdim(sol13.y(1,1:last)); q13(1:last) = q3*sol13.x(1:last);
  % the exponential of z - z3 is plotted;
  z13 = z3 + (sol13.y(2,:)-1)*zscale;
  writetostruct('13-',{'T','p','q','z','color'},{T13,p13,q13,z13,vapcolor});
end

%  THE END
fl.sol.p1 = p1;  fl.sol.T1 = T1;  fl.sol.T3 = T3;
fl.sol.q1 = q1; % here, but not above, q1 is always zero

% OBSOLETE FLOW13
% m h + q = const => m dh/dz + dq/dz = 0; dh/dz = cp dT/dz ( + dh/dp dp/dz = 0).
% Therefore, with q = -k dT/dz we have
%   dq/dz = - m cp dT/dz, or   dq/dz = m*cp/k q,
%   dT/dz = - q/k.
% The analytical solution for k, cp = const gives an exponential decay to
% z -> -infty:   m cp dT/dz - k d^2T/dz^2 = 0;   lambda^2 - m*cp/k lambda = 0;
% lambda1 = 0, lambda2 = m*cp/k. Hence T = Ta + Tb*exp(m*cp/k z). With boundary
% conditions T(z=0) = T3, dT/dz(z=0) = -q3/k yields Tb = -q3/m*cp and the
% solution is
%    T = T3 + q3/m*cp3 ( 1 - exp(m*cp/k z) ).
% To get a decent numerical solution, we scale the z-coordinate (for z3 = 0):
%   zw = exp(m*cp3/k3 z),  dzw = m*cp3/k3 exp(m*cp3/k3 z) dz,
%                          dz = dzw/(m*cp/k3 zw)
%   qw = q/q3,  dq = q3 dqw;  (q3 may be negative!)
%   Tw = (T - T3) / (q3/m*cp3),  dT = q3/m*cp3 dTw.
% Hence, the differential equations become in dimensionless form
%   dqw/dzw = cp/cp3 k3/k qw/zw,  (dqw/dzw = -cp/cp3 dTw/dzw),
%   dTw/dzw = -k3/k qw/zw 
% with initial conditions
%   Tw(zw=1) = 0,  qw(zw=1) = 1.
% Test this: For k = k3, c = cp3 the analytical solution must be recovered. The
% dimensionless eqs. are solved by qw = zw, then dqw/dzw = 1 and dTw/dzw = -1.
% Integrate Tw to get T = -zw + c1, bc. gives: Tw = 1 - zw. Made dimensional
% this correctly yields the analytical solution above.
% For z3 arbitrary the scaling is   z = ln(zw) k3/m*cp3 + z3.
%
%options=odeset('RelTol',rtol,'Refine',1,'InitialStep',step13,'MaxStep',step13);
%%sol13 = ode45(@int13w,[zw3 0],[Tw3 qw3],options);
%%sol13 = ode45(@int13w,[1 0],[0 1],options); % 0 gives singularity;
%sol13 = ode45(@int13w,[1 0.01],[0 1],options);
%
%function dy = int13w(z,y)
%% dimensionless eqs. see above
%  Tw = y(1);  qw = y(2);
%  T = mkTdim(Tw);
%  dTw = -k3*qw/(z*s.kg(T));
%  dqw = -s.cpg(T,p3)*dTw/cp3;
%  dy = [dTw; dqw];
%end
%
% END OBSOLETE FLOW13
end %---------------------------------------------------------------- end flow13

function flow45(m,T5,p5,q5) %-------------------------------------------- flow45
%FLOW45     Liquid film in front of the membrane.
%  | -> front34 -> flow13

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
% Instead of zw, i could use T as independent variable and integrate between T5
% and T4. Advantage: no termination condition necessary and i do not need to
% guess a large enough range of z.
%   dq/dT = -m cp,   dqw/dTw = -cp m*(T4-T5)/q5
%   dz/dT = -k/q,  dzw/dTw = -k/k5 1/qw.

%  ASSIGN LAST POINT
% we already know part of the solution; needed for characteristic scales.
p4 = p5; T4 = s.Ts(p5);

%  CHARACTERISTIC SCALES
% scales to make dimensionless
k5 = s.kl(T5); % cp5 = s.cpl(T5); probably not needed
Tscale = T4 - T5; % qscale = q5; % must be q5 > 0
zscale = Tscale*k5/q5;
mkTdim = @(Tw) Tscale.*Tw + T5;
coeff1 = -m*Tscale/q5;
step45 = -odemaxstep(Tscale,maxTperstep);

%  INTEGRATE
% integrate from Tw5 to Tw4, i.e. from 0 to 1.
options = odeset('RelTol',rtol,...
  'Refine',1,'InitialStep',step45,'MaxStep',step45);
%sol45 = ode45(@int45,[Tw5 Tw4],[qw5 zw5],options);
sol45 = ode45(@int45,[0 1],[1 0],options);

function dy = int45(Tw,y)
  T = mkTdim(Tw);
  dqw = coeff1*s.cpl(T);
  dzw = -s.kl(T)/(k5*y(1));
  dy = [dqw; dzw];
end

%  ASSIGN LAST POINT
% dimensionalize and assign actual state, here T4, p4, q4, z4
% T4, p4 assigned already above
last = size(sol45.x,2);
q4 = q5*sol45.y(1,last); z4 = zscale*sol45.y(2,last);

%  WRITE SOLUTION
% if wanted, write the solution
if writesolution
  T45(last) = T4; p45(1:last) = p4; q45(last) = q4; z45(last) = z4;
% i could do T45(1) = T5; T45(2:last-1) = ...
  last = last - 1;
  T45(1:last) = mkTdim(sol45.x(1:last)); q45(1:last) = q5*sol45.y(1,1:last);
  z45(1:last) = zscale*sol45.y(2,1:last);
  writetostruct('45-',{'T','p','q','z','color'},{T45,p45,q45,z45,liqcolor});
end

%  CALL NEXT
T3 = T4; p3 = p4; z3 = z4; q3 = q4 - m*s.hvap(T3); %==================== front34
flow13(m,T3,p3,q3,z3);

% OBSOLETE
%%  INTEGRATE
%% integrate; dimensionless initial conditions
%% qw stays close to 1, no AbsTol (as in flow92) needed.
%options=odeset('Events',@term45w,'RelTol',rtol);
%%sol56 = ode45(@int56w,[0 filmthickness?],[T5w q5w],options);
%sol45 = ode45(@int45w,[0 -2],[0 1],options);
%%
%function dy = int45w(z,y)
%  T = mkTdim(y(1)); qw = y(2);
%  dTw = -qw*k5/s.kl(T);
%  dqw = s.cpl(T)*coeff1*dTw;
%  dy = [dTw; dqw]
%end
%%
%function [val,isterm,direction] = term45w(z,y)
%isterm = 1; direction = 0; % direction of zero-crossing does not matter
%% dimensionless, Ts(p5) = T4, hence T4w = 1
%val = y(1) - 1;
%end
% END OBSOLETE
end %---------------------------------------------------------------- end flow45

function [T3 p3 q3] = front35(m,T5,p5,q5) %----------------------------- front35
%FRONT35    Condensation at upstream end of membrane.
% Kelvin's eq. + Young-Laplace eq.; radius of curvature must be found.
%   ln(p3/psat(T3)) = -pcap/(s.rho*R*T); pcap = p3 - p5;
% Hence, solve the implicit eq.
%   ln(p3/psat(T3)) = -(p3-p5)/(s.rho*R*T);
% Initial guess: With the first order approx. ln(p3/psat3) = p3/psat3 - 1,
% p3/psat3 -1 = (p5-p3)/(s.rho*R*T), hence
%   p3/psat3 =  (p5 + rho*R*T) / (psat3 + rho*R*T).
% Newton iteration:
%   F(x) = ln(p3/psat3) + (p3-p5)/rho*R*T,   F'(x) = 1/x + psat3/rho*R*T.
% hvapK .. see rkelv; Eq. (14) in [Loimer, 2005];
% See also front62.
T3 = T5;
psat3 = s.ps(T3);  [drho3 rho3] = s.drho(T3); 
rhoRT = rho3*s.R*T3;
% inital guess 
p3 = (p5 + rhoRT)/(psat3+rhoRT); % p3 here really is p3/psat3
% setup newton
ps_rhoRT = psat3/rhoRT; p5_rhoRT = p5/rhoRT;
function [F dF] = sol35(x)
  F = log(x) + x*ps_rhoRT - p5_rhoRT; % F = log(x) + (x*psat3-p5)/rhoRT;
  dF = 1/x + ps_rhoRT;                % dF = 1/x + psat3/rhoRT;
end
p3 = newton(@sol35,p3,1e-12)*psat3; % ln of x approx. 1 - res must be very small
% evaporation enthalpy
hvapK3 = hvapKraw(T3,p3,psat3,p3-p5,rho3,drho3);
q3 = q5 - m*hvapK3;
end %--------------------------------------------------------------- end front35

function [T8 a8 doth8 pk8 dpkdT8 hvapK8] = front85(m,T5,p5,q5) %-------- front85
%FRONT85    Full condensation of two-phase mixture within the membrane.
% Find a8 that solves 
%    m*xdot*hvapK + q8 =  q5.
% With q8 =  k*nu*m/((dpk/dT)*kap), (see front89), we have
%  0 =! q5/m - xdot*hvapK - k*nu/((dpk/dT)*kap).
% doth8 = q8 + m*xdot*hvapK8 = q5
T8 = T5; doth8 = q5;
[p8 dpkdT8 hvapK8] = hvapK(T8);
sol85 = @(a) q5/m - xdot(T8,p8,a)*hvapK8 ...
   - k2ph(T8,a)*nu2ph(T8,p8,a)/(dpkdT8*mem.kappa);
% DEBUG
%try
a8 = fzero(sol85,[0 1],optimset('TolX',tola));
%catch
%a8 = [0:0.01:1]; sola8 = a8;
%for i = 1:length(a8), sola8(i) = sol85(a8(i)); end
%disp(sprintf('min sola8: %g, max sola8: %g', min(sola8), max(sola8)));
%plot(a8,sola8,'k-');
%end
end %--------------------------------------------------------------- end front85

function writetostruct(name,vars,values) %------------------------ writetostruct
%WRITETOSTRUCT Write the solution to the flowstruct.

% COMMON VARIABLE FL
% some information on the integration path
fl.sol.states = [name fl.sol.states];

fl.sol.len = fl.sol.len - 1;
j = -fl.sol.len;
last = size(vars,2);
for i = 1:last
  fl.flow(j).(vars{i}) = values{i};
end
end %--------------------------------------------------------- end writesolution

end %%% END FLOW12 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END FLOW12 %%%


%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%

function x = newton(fun,x0,res,iter) %----------------------------------- newton
%NEWTON     Newton iteration. Controls tolerance in the solution.
%
%  NEWTON(FUN,X0) finds a solution to F(X) = 0 by Newton iteration, starting
%  from X = X0. NEWTON expects FUN to be the handle to a function that returns a
%  vector [F(X) DF(X)], where DF is the derivative of F at X.
%
%  NEWTON(FUN,X0,RES,ITER) finds a solution to F(X) = 0. NEWTON iterates until
%  X changes less than RES in one step or if more than ITER iterations are done.
%  Default values are RES = 1e-6 and ITER = 100.
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
if nargin < 3, res = 1e-6; end

x = x0;
for i = 1:iter
  [y dy] = fun(x);
  dx = y/dy;
  if abs(dx) < res, return
  else
    x = x - dx;
  end
end

% Be verbose if no solution is found.
funname = func2str(fun);
error(['NEWTON not succesful. Increase RES or ITER. Type help newton.\n' ...
  '  Function %s, initial guess x0: %g, last found value x: %g,\n' ...
  '  function value %s(x) = %g. Allowed residual: %g, Iterations: %u.'], ...
  funname,x0,x,funname,y,res,iter)
end %-----------------------------------------------------------------end newton
