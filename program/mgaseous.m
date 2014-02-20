function [m,ms] = mgaseous(T,p1,p2,s,ms,type,accuracy)
%MGASEOUS    Mass flux for isothermal gaseous flow.
%  MGASEOUS(T1,P1,P2,SUBSTANCE,MS) returns the isothermal, gaseous mass flux.
%  Viscous and free molecular flow are taken into account. MS must contain a
%  field MS.SUBSTANCE, apart from being constructed with MSTACKSTRUCT. Accurate
%  solver settings are used.
%
%  [M,MS] = MGASEOUS(T1,P1,P2,SUBSTANCE,MS) writes the solution to MS.
%
%  [M,MS] = MGASEOUS(T1,P1,P2,SUBSTANCE,MS,TYPE,'coarse') calculates the
%  isothermal gaseous mass flux according to TYPE. Type is either 'viscous',
%  'gaseous' or 'knudsen'. With 'coarse', coarser solver settings are used.
%  Accurate solver settings are set with 'accurate'.
%
%  See also MSTACKSTRUCT, SUBSTANCE.

if s.ps(T) < p1
  error('The usptream state is a liquid. Not implemented.');
end

if nargin < 7
  accuracy = 'accurate';
  if nargin < 6
    type = 'gaseous';
  end
end

% Copy some values to the membrane struct
ms.T1 = T;
ms.p1in = p1;
ms.T2 = T;
ms.p2 = p2;
ms.a2 = 1;
ms.q2 = 0;

ms.substance = s;

mguess = ms.mfluxviscous(T,p1,p2,s,ms);
if strcmp(type,'gaseous')
  % T2 = 1000 - surely supercritical
  ms = ms.writeflowsetups(T,1000,s,ms);
end

% Set up the solver and solution iteration
solver = solverstruct(accuracy);
solver.gasflow = type;
presiduum = @(m) viscous(m,T,p2,ms,solver) - p1;
[minterval,pinterval] = findinterval(presiduum, mguess, p2-p1);

m = findzero(presiduum,[minterval; pinterval],(p1-p2)/10000);
%fprintf('Mass flux guessed mguess = %g, calculated m = %g, mguess - m = %g%%\n',...
%  mguess, m, 100*(mguess-m));

end %%% END MGASCOUS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MGASCOUS %%%

%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%

function [p1,ms] = viscous(m,T,p2,ms,solver) %-------------------------- viscous
%VISCOUS    Viscous gaseous flow through a porous medium with several layers.
%  VISCOUS(M,T,P2,MSTACKSTRUCT,SOLVER) returns the upstream pressure for the
%  gaseous, viscous flow of a fluid with temperature T and downstream pressure
%  P2 through an assemblage of membranes described by
%  at a downstream state STATE through an assemblage of membranes described by
%  the membrane struct MSTACKSTRUCT. The solution tolerances are controlled via
%  the struct SOLVER.
%
%  [P1,MSTACKSTRUCT] = VISCOUS(M,T,P2,MSTACKSTRUCT,SOLVER)) returns the upstream
%  pressure P1 and a MSTACKSTRUCT describing the solution.

%  TODO - rewrite below here, after first runs
%  Restrictions: The upstream state must be a gaseous phase, because the liquid
%  film and the temperature boundary layer integrators upstream of the membrane
%  (ifreevapor, ifreeliquid) use termination conditions that work only with an
%  gaseous upstream state.
%  The vapor phase integrator (integratevapor) always starts at the downstream
%  end of a layer.
%  For a two-phase upstream state, integration direction would have to be
%  downstream.
%
%  ASYM>FRONT should probably be a subfunction or nested function in STATE,
%  STATE1 = STATE2.FRONT(STATE2,...). However, algorithms in FRONT and INTEGRATE
%  are partially similar, therefore they are kept in one place.
%
%  See also DOWNSTREAMSTATE, MSTACKSTRUCT, SOLVERSTRUCT.

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
%freesetup = ms.freesetup; % free space flow setup

% integrate in upstream direction
% Cycle over all membranes
for i = nmembranes:-1:1
  nlayers = length(ms.membrane(i).layer);

  % Cycle over all layers in one membrane
  for j = nlayers:-1:1
    % Start is the downstream end of the last layer of the last membrane;
    % flow grows in integrate(), if solver.writesolution is true
    flow = zeroflowstruct;
    [p2,flow] = integratevapor(m,T,p2,flow,ms.membrane(i).layer(j).matrix,...
				ms.membrane(i).layer(j).flsetup,s,solver);
    % Arrived at the upstream end of a layer
    ms.membrane(i).layer(j).flow = flow;
  % Go into the next layer
  end

% Next membrane
end
p1 = p2;
end %--------------------------------------------------------------- end viscous

function [p9,flow] = integratevapor(m,T,p2,flow,mem,fs,s,solver) %integratevapor
%INTEGRATEVAPOR Vapor flow within the membrane - a copy of FLOW12>FLOW92

% TODO: z is not referred to! It is assumed z = mem.L.

% integrate from z = L to z = 0, initial condition T, p2.
% eqs.: Darcy's law,
%   dp/dz = -m nu/kappa.
% Integrate dimensionless eqs. (w, instead of star), with
%   zscale = mem.L,
% hence
% pscale/zscale dpw/dzw = -m nu/kap,
%   dpw/dzw = -m*mem.L/kappa*pscale nu,  O(pw) = 1, O(dpw/dzw) = m*zscale/...

if strcmp(solver.gasflow,'gaseous')
  nu = fs.nuapp;
else % 'viscous'
  nu = s.nug;
end
%  CHARACTERISTIC SCALES
% scales to make dimensionless
pscale = m*nu(T,p2)*mem.L/mem.kappa;
step92 = -min(solver.odemaxstep(pscale,solver.maxpperstep),solver.maxpperstep/pscale);
% functions to re-calculate dimensional values
mkpdim = @(pw) pscale.*pw + p2;
% coefficient for the eq. in int92w:
coeff = -m*mem.L/(mem.kappa*pscale);

%  INTEGRATE
% make initial conditions dimensionless; integrate
% p2w = 0; z2w = 1;
options=odeset('RelTol',solver.rtol,'Refine',1,'InitialStep',step92,'MaxStep',step92);
%sol92 = ode45(@int92w,[z2w 0],[T2w p2w q2w],options);
sol92 = ode45(@int92w,[1 0],0,options);

function dy = int92w(z,y)
  % dimensionless eqs., cf. above.
  % pw = y;
  p = mkpdim(y);
  % dy = dpw;
  dy = nu(T,p)*coeff;
end

%  ASSIGN LAST POINT
% dimensionalize and assign actual state, here T9, p9, q9, z9
last = size(sol92.x,2);
p9 = mkpdim(sol92.y(last));
if sol92.x(last) ~= 0
  error('Integratevapor did not reach upstream end of layer, zend/L = %g',...
	sol92.x(last));
end
z9 = sol92.y(last);

%  WRITE SOLUTION
% if wanted, write the solution
if solver.writesolution
  % allocate space for all points; assign last point
  p92(last) = p9;
  p92(1:last-1) = mkpdim(sol92.y(1:last-1));
  z92(1:last) = sol92.y(1:last);
  T92(1:last) = T; q92(1:last) = 0;
  % and write the solution but the last point
  flow = writeflow(flow,{'z','T','p','a','q','color'},...
			{z92,T92,p92,ones(1,last),q92,'r'});
end

end %-------------------------------------------------------- end integratevapor

function flow = writeflow(flow,vars,values) %------------------------- writeflow
%WRITEFLOW  Write the solution to the flowstruct.

j = length(flow) + 1;
last = size(vars,2);
for i = 1:last
  flow(j).(vars{i}) = values{i};
end
end %------------------------------------------------------------- end writeflow
