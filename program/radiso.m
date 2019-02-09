function [r9,r3,ms] = radiso(m,T1,p1,p2,s,ms)
%RADIUS     Radii of curvature at upstream and downstream meniscus.
%  RADIUS(M,STATE,MSTACKSTRUCT,SOLVER) returns the radius of curvature of
%  the downstream meniscus, given a mass flux M of MSTACKSTRUCT.SUBSTANCE
%  through an single homogeneous membrane described by the the struct
%  MSTACKSTRUCT. The downstream state of the fluid is given by STATE, the
%  solution tolerances are controlled via the struct SOLVER.
%
% TODO complete help / description
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
%  See also DOWNSTREAMSTATE, MSTACKSTRUCT, FLOWSETUP, SOLVERSTRUCT.

if length(ms.membrane) > 1
    error([upper(mfilename) ': Only one homogeneous membrane possible.']);
end

% States:
%       _______________________
%  1   3)5 liquid   6(9        2
%       ^^^^^^^^^^^^^^^^^^^^^^^
% Here: isothermal, T = constant, T1 = T3, because of p1 = p3
% state 1 = state 3.

% some material properties
psat = s.ps(T1);
T2 = T1;
rhol = s.rho(T1);
nuliq = s.mul(T1)/rhol;;
sigma = s.sigma(T1);

% a dimension
L =  ms.membrane(1).layer(1).matrix.L;

if p1 > psat
    error([upper(mfilename)...
	': P1 is larger than the saturation pressure, must be smaller.']);
end

zeroflowstruct = struct('z',{},'T',{},'p',{},'a',{},'q',{},'Kn',{},'color',{});
solver = solverstruct('accurate');
solver.writesolution = true;

% integrate the vapor flow in upstream direction
[p9,z9,flow,sol92,mkp9dim] = integratevapor(m,T2,p2,psat,zeroflowstruct,...
	ms.membrane(1).layer(1).matrix, ms.membrane(1).layer(1).flsetup, s,...
	solver);

if z9 < L/1000. && p9 <= p1
    % vapor flow all through, and the upstream pressure is not reached - no
    % possibility for the flow of condensate, the flow remains gaseous
    ms.membrane(1).layer(1).flow = flow;
    r3 = 0;
    r9 = 0;
    fprintf('gaseous flow only: m = %.1f kg/m²s\n', m*1000);
    return
end

% get the state of the liquid at the upstream front
% Kelvin's eq.,
% ln(pk/psat) = -(2 sigma/r) (v/RT),  pk = p2, ...
pcap3 = log(psat/p1).*rhol.*s.R.*T1;
r3 = 2.*sigma./pcap3;
p5 = p1 - pcap3;


% liquid flow is given by
%    dp/dz = -m nul/kappa,
% with the right hand side being constant.
% the liquid pressure in the liquid flow region
pliqliq = @(z) p5 - z.*m.*nuliq./ms.membrane(1).layer(1).matrix.kappa;

function p96 = pres(z)
    p9 = mkp9dim(deval(sol92,z/L));
    % the liquid pressure, from the vapor flow region,  is p9 - pcap,
    % pcap = ln(psat/p9)*rho*R*T, see above
    % pliq (from the vapor side) minus pliq (from the liquid side)
    p96 = p9 - log(psat/p9)*rhol*s.R.*T1 - pliqliq(z);
end

% +/- one thousandth - once, an error occured whereupon sol92 was evaluated
% outside the interval where it is defined; Therefore, reduce the interval
% by just a tiny margin.
z69 = fzero(@pres, [z9+L/1000 L-L/1000]);

r9 = 2*sigma/(mkp9dim(deval(sol92,z69/L)) - pliqliq(z69));

% TODO: would now have to write the vapor flow solution to the flowstruct

% write the solution for the liquid flow to the flowstruct
% gaseous and liquid flow now overlap
flow = solveliquid(m, T1, p5, z69, flow, ms.membrane(1).layer(1).matrix,...
		ms.membrane(1).layer(1).flsetup, s,solver);
ms.membrane(1).layer(1).flow = flow;

% write the value in front - only one value
ms.membrane(1).flow = zeroflowstruct;
ms.membrane(1).flow(1).z = 0;
ms.membrane(1).flow(1).T = T1;
ms.membrane(1).flow(1).p = p1;
ms.membrane(1).flow(1).a = 0;
ms.membrane(1).flow(1).q = 0;
ms.membrane(1).flow(1).Kn = 0;
ms.membrane(1).flow(1).color = 'r';

zeroflowstruct = struct('z',{},'T',{},'p',{},'a',{},'q',{},'Kn',{},'color',{});

% write the solution to the membranestruct
ms.m = m;
ms.T1 = T1;
ms.p1in = p1;
ms.p1sol = p1;
ms.a1 = 1;
ms.q1 = 0;
ms.T2 = T1;
ms.p2 = p2;
ms.a2 = 0;
ms.q2 = 0;


end %%% END RADISO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END RADISO %%%


%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%

%---------------------------------------------------------------- integratevapor
function [p9,z9,flow,sol92,mkpdim] = integratevapor(m,T2,p2,psat2,flow,...
							mem,fs,s,solver)
%INTEGRATEVAPOR Unconditionally integrate vapor flow within the membrane.
%
%  A copy of ISO->INTEGRATEVAPOR without termination condition.
%  ISO->INTEGRATEVAPOR is itself a copy of ASYM->INTEGRATEVAPOR.
%
%  See also ASYM>INTEGRATEVAPOR, ISO->INTEGRATEVAPOR.

% TODO: z is not referred to! It is assumed z = mem.L.

% integrate from z = L to z = 0, initial conditions T2, p2, q2.
% eq.: Darcy's law,
%   dp/dz = -m nu/kappa,
%   zscale = mem.L,
% Dimensionless,
%   z = L*zw;  p = pscale*pw + p2;
% dimensionless eqs:
% pscale/zscale dpw/dzw = -m nu/kap,
%   dpw/dzw = -m*zscale/kappa*pscale nu,  O(pw) = 1, O(dpw/dzw) = m*zscale/...

%  CHARACTERISTIC SCALES
% scales to make dimensionless
%pdarcy = m*s.nug(T2,p2)*mem.L/mem.kappa;
pdarcy = m*fs.nuapp(T2,p2)*mem.L/mem.kappa;
pscale = min(psat2-p2, pdarcy);
step92 = -min(solver.odemaxstep(pscale,solver.maxpperstep),...
	solver.maxpperstep/pscale);
init92 = step92;

% functions to re-calculate dimensional values
mkpdim = @(pw) pscale.*pw + p2;

% dimensionless saturation pressure
ps2w = (psat2 - p2)/pscale;

function [p, z] = mkdimensional(pw, zw)
  p = mkpdim(pw);
  z = mem.L.*zw;
end

% coefficient for the eq. in int92w:
coeff = -m*mem.L/(mem.kappa*pscale);

%  INTEGRATE
% make initial conditions dimensionless; integrate
options=odeset('Events',@term92w,...	%'AbsTol',atol*[1 1 -coeff1*cpg2],...
  'RelTol',solver.rtol,'Refine',1,'InitialStep',step92,'MaxStep',step92);
sol92 = ode45(@int92w, [1 0], 0, options);

function dy = int92w(z,y)
  % dimensionless eqs., cf. above.
  p = mkpdim(y);
  dy = fs.nuapp(T2,p) * coeff;
end

% termination condition
function [val,isterm,direction] = term92w(z,y)
  isterm = 1;
  direction = 0; % 0 .. all zero crossings, -1 .. event function decreasing
  val = ps2w - y;	% psat2 - p, must be > 0
end

%  ASSIGN LAST POINT
% dimensionalize and assign actual state, here T9, p9, q9, z9
last = size(sol92.x,2);
T9 = T2;
q9 = 0;
[p9, z9] = mkdimensional(sol92.y(last), sol92.x(last));

%  WRITE SOLUTION
% if wanted, write the solution
if solver.writesolution
  % allocate space for all points; assign last point
  T92(1:last) = T2;
  q92(1:last) = q9;
  p92(last) = p9; z92(last) = z9;
  Kn92(last) = fs.knudsen(T9,p9);
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

%------------------------------------------------------------------- solveliquid
function flow = solveliquid(m,T5,p5,z6,flow,mem,fs,s,solver)
%SOLVELIQUID Solve for liquid flow within the membrane.
%  The liquid viscosity does not change with pressure, hence create a
%  straight line and fill the solution struct.

%   dp/dz = -m nul/kappa

% Assume, that z5 = 0.
p6 = p5 - z6*m*s.nul(T5)/mem.kappa;

%  WRITE SOLUTION
% if wanted, write the solution
if solver.writesolution
  z56 = [z6 0];
  T56 = [1 1] * T5;
  p56 = [p6 p5];	% the flow struct is reverse
  q56 = [0 0];
  Kn56 = [0 0];
  flow = writeflow(flow,{'z','T','p','a','q','Kn','color'},...
			{z56,T56,p56,q56,q56,Kn56,'b'});
end

end %----------------------------------------------------------- end solveliquid

function flow = writeflow(flow,vars,values) %------------------------- writeflow
%WRITEFLOW  Write the solution to the flowstruct.

j = length(flow) + 1;
last = size(vars,2);
for i = 1:last
  flow(j).(vars{i}) = values{i};
end
end %------------------------------------------------------------- end writeflow
