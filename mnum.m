function [m fl] = mnum(T1,p1,p2,theta,s,mem,f,Tmax,mguess)
%MNUM       Mass flux by numerical integration.
%  MNUM(T1,P1,P2,THETA,S,MEM,F) calculates the mass flux [kg/m2s] for a
%  substance S through a membrane MEM with a two-phase flow model F. The
%  upstream conditions are T1 [K] and P1 [Pa], the downstream pressure is P2
%  [Pa] and the contact angle is THETA, in degrees.
%
%  [M FL] = MNUM(T1,P1,P2,THETA,S,MEM,F) returrns the mass flux M and a
%  flowstruct FL containing the solution.
%
%  MNUM(T1,P1,P2,THETA,S,MEM,F,'m',MGUESS) uses MGUESS as the inital guess for
%  the mass flux.
%
%  MNUM(M,T2,P2,THETA,S,MEM,F,TMAX) constructs a struct SOLVER and invokes
%  FLOW12. Returns a FLOWSTRUCT. See FLOW12.
%
%  Calls FLOW12.
%
%  See also FLOW12, FMODEL, MEMBRANE, SUBSTANCE.
%
%  Subfunctions:  FLOWSETUP, FLOWSTRUCT.
%  Try, e.g., help mnum>flowsetup.
%
%  Nested functions: accurate, crude, findzero, presiduum.

%  INITIALIZATION
% Construct and initialize the flow struct.
fl = flowstruct(theta,s,mem,f);
fl.info.p2 = p2;
% Construct the solver struct.
solver = struct('rtol',[],'atol',[],'tola',[],...
  'maxTperstep',[],'maxpperstep',[],'writesolution',[],'partialsolution',[]);

%  SOLVER SETUPS
% Set solver to accurate or crude tolerances.
function solver = accurate(solver) %----------------------------------- accurate
  solver.rtol = 1e-3;  solver.atol = 1e-6;  solver.tola = 1e-6;
  solver.maxTperstep = .4;  solver.maxpperstep = 4e4;
end %-------------------------------------------------------------- end accurate
function solver = crude(solver) %----------------------------------------- crude
%  solver.rtol = 5e-2  <--  Das geht nicht für flow78, non-wetting.
  solver.rtol = 1e-3;  solver.atol = 1e-5;  solver.tola = 1e-6;
  solver.maxTperstep = 2;  solver.maxpperstep = 1e5;
end %----------------------------------------------------------------- end crude

%  CALL FLOW12 ONLY
% Either setup solver, get thermodynamic properties, call flow12(m,..) and
% return. Or setup solver, get TD properties and continue.
if nargin == 8 % Tmax defined
  % mnum called as flow12, mnum(m,T2,p2,theta,s,mem,f,Tmax)
  m = T1; T2 = p1;
  fl.info.m = m; fl.info.T2 = T2;
  flsetup = flowsetup(T2,Tmax,theta,s,mem,f);
  % Assign values for rather accurate calculation.
  solver = accurate(solver);
  solver.writesolution = true;  solver.partialsolution = ~solver.writesolution;
  fl = flow12(m,T2,p2,fl,flsetup,solver);
  % only one output if mnum is called as flow12
  m = fl;
  return
end % end nargin == 8

%  INITIALIZE FOR MNUM PROPER
% Assign T2, Temperature integration range, solver and call tdproperties.
%[ps dps] = s.ps(T1);
if  p1 > s.ps(T1)
  error(['MNUM: Fluid is in liquid state at upstream location:'...
    ' p1 > psat(T1), p1 = %g bar, psat = %g bar'],p1*1e-5,s.ps(T1)*1e-5);
end
p12 = p1 - p2;
T2 = s.intjt(T1,p1,p2);
fl.info.T1 = T1; fl.sol.T2 = T2;
fl.info.p1 = p1;
flsetup = flowsetup(T2,T1,theta,s,mem,f);
% Here, one could ...
% Calculate the information already possible: 
%[dpk pk] = dpkdT(T);
%fl.info.kappac = s.nul(T1)*flsetup.kmliq(T1)/(dps*s.hvap(T1));
%if theta ~= 90
%  if theta < 90
%    fl.info.Ccc = 
%end
% How much output is wanted.
if nargout > 1 % mnum called with [m fl] = mnum(...
  solver.writesolution = true;
else
  solver.writesolution = false;
end
solver.partialsolution = ~solver.writesolution;

%  FIRST GUESS
% Guess the mass flow either from linear theory or assuming pure vapor flow.
if nargin < 9
  if p1 >= flsetup.pkelv(T1)
    % for nearly saturated vapor use linear theory
    flin = mlinear(p1,p2,T1,theta,s,mem,f);
    mguess = flin.lin.m;
  else
    % for unsaturated vapor or gas, use pure vapor flow as a guess
    % mliq = p12*mem.kappa/(s.nul(T1)*mem.L);
    mguess = p12*mem.kappa/(flsetup.nuapp(T2,(p1+p2)/2)*mem.L);
  end
end

%  SHOOT AT P1
p1tol = 1e-3*p12;
m = findzero(mguess,p1tol);
fl.sol.m = m;

%% -- garbage 1

%%% NESTED FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NESTED FUNCTIONS %%%

function pres = presiduum(m,solver) %--------------------------------- presiduum
%PRESIDUUM  Calculate the residual pressure for a given mass flux M.
  fl = flow12(m,T2,p2,fl,flsetup,solver);
  pres = fl.sol.p1 - p1;
end %------------------------------------------------------------- end presiduum

function b = findzero(mguess,p1tol) %---------------------------------- findzero
%FINDZERO  Find a solution to PRESIDUUM(M) = 0.

%  INITIALIZE
% Information is plotted for trace > 1. Iteration results for trace > 2.
trace = 1;
fcount = 0; procedure = ' ';
savewritesolution = solver.writesolution;

%  ASSIGN CRUDE SOLVER TOLERANCES
% During iteration, use crude tolerances.
solver = crude(solver); solver.writesolution = false;
solver.partialsolution = ~solver.writesolution;

a(1) = mguess; fa(1) = presiduum(a(1),solver); fcount = 1;

%  SEARCH AN INTERVAL
% Now look, in which direction we have to search for pres =  0.
if fa(1) > 0, fac = 0.5; else fac = 1.4; end

for i = 2:6
  old = i - 1;
  a(i) = a(old)*fac;
  fa(i) = presiduum(a(i),solver);
  fcount = fcount + 1;
  if (fa(old) < 0) == (fa(i) > 0) || fa(i) == 0
    % found an interval.
    a = [a(i) a(old)];
    fa = [fa(i) fa(old)];
%DEBUG disp(sprintf('Intervall: %8f < m < %8f, %8f < pres < %8f',...
%DEBUG   a(1),a(2),fa(1),fa(2)));
    break
  end
end

%  PLOT SEARCH, IF NO INTERVAL FOUND
if size(a,2) ~= 2 %DEBUG || true(); %DEBUG
  % No interval found.
  % also plot the result for mvap and mliq
  mvap = p12*mem.kappa/(flsetup.nuapp(T2,(p1+p2)/2)*mem.L);
  pvap = presiduum(mvap,solver);
  % No, really no interval found.
  disp(sprintf('MNUM: No interval found! mguess = %g',mguess));
  %disp('  Call MNUM with a guess for the mass flow,');
  %disp('  mnum(T1,p1,p2,theta,s,mem,f,''m'',mguess)');
  disp('  Mass flow:'); disp(a);
  disp('  p_res :'); disp(fa);
  % reuse fac as figure handle
  fac = figure('Name','Search for Interval','NumberTitle','off');
  plot(a,fa/p1,'ko-',mvap,pvap/p1,'r*'); %,mliq,pliq/p1,'b*');
  xlabel('mflux'); ylabel('p_{1,calc}/p_1 - 1');
  % AND ASK FOR AN INTERVAL
  while size(a,2)~= 2
    b = input(['  Provide an interval [mflux mflux] or try\n' ...
      '  several points [mflux mflux mflux ...]: ']);
    lenb = size(b,2); fb = b;
    for i = 1:lenb,  fb(i) = presiduum(b(i),solver);  end
    % reuse old as figure handle; save a possible different current figure.
    old = get(0,'CurrentFigure'); figure(fac);
    line(b,fb/p1,'Marker','+','LineStyle','none');
    figure(old);
    if lenb == 2 && (fb(1) < 0) == (fb(2) > 0)
      a = b; fa = fb;
    end
  end
  close(fac);
end

%  INTERVAL FOUND OR ASSIGNED
% check for zero at boundaries
b = a(2); a = a(1); fb = fa(2); fa = fa(1);
lowprecision = true;

% if a or b == 0, the while-loop is run once, lowprecision set to false, and the
% result re-calculated.

tol = eps; % we do not care for accuracy in m
fc = fb;
procedure = 'initial';
header2 = ' Func-count    x          f(x)             Procedure';
if trace > 2
    disp(header2)
end
% Main loop, exit from middle of the loop
while fb ~= 0 && a ~= b
    % Insure that b is the best result so far, a is the previous
    % value of b, and c is on the opposite side of the zero from b.
    if (fb > 0) == (fc > 0)
        c = a;  fc = fa;
        d = b - a;  e = d;
    end
    if abs(fc) < abs(fb)
        a = b;    b = c;    c = a;
        fa = fb;  fb = fc;  fc = fa;
    end
    
    % Convergence test and possible exit
    m = 0.5*(c - b);
    toler = 2.0*tol*max(abs(b),1.0);
    if (abs(fb) < p1tol) || (abs(m) <= toler) % || (fb == 0.0) 
        if lowprecision
	  if trace > 1
	    disp(sprintf('p1 = %g, p2 = %g [bar]. We jumped over 100xp1tol.',...
	      p1/1e5,p2/1e5));
	  end
	else
          break
	end
    end
    if trace > 2
        disp(sprintf('%5.0f   %13.6g %13.6g        %s',fcount, b, fb, procedure));
    end
    
    % Choose bisection or interpolation
    if (abs(e) < toler) || (abs(fa) <= abs(fb))
        % Bisection
        d = m;  e = m;
        procedure='bisection';
    else
        % Interpolation
        s = fb/fa;
        if (a == c)
            % Linear interpolation
            p = 2.0*m*s;
            q = 1.0 - s;
        else
            % Inverse quadratic interpolation
            q = fa/fc;
            r = fb/fc;
            p = s*(2.0*m*q*(q - r) - (b - a)*(r - 1.0));
            q = (q - 1.0)*(r - 1.0)*(s - 1.0);
        end;
        if p > 0, q = -q; else p = -p; end;
        % Is interpolated point acceptable
        if (2.0*p < 3.0*m*q - abs(toler*q)) && (p < abs(0.5*e*q))
            e = d;  d = p/q;
            procedure='interpolation';
        else
            d = m;  e = m;
            procedure='bisection';
        end;
    end % Interpolation
    
    % Next point
    a = b;
    fa = fb;
    if abs(d) > toler, b = b + d;
    elseif b > c, b = b - toler;
    else b = b + toler;
    end
    if lowprecision && abs(fb) < 100*p1tol
      solver.writesolution = savewritesolution;
      solver.partialsolution = ~savewritesolution;
      solver = accurate(solver);
      lowprecision = false;
      if trace > 1
	disp(sprintf('p1 = %g, p2 = %g. Switched to high precision.',...
	  p1/1e5,p2/1e5));
      end
    end
    fb = presiduum(b,solver);
    fcount = fcount + 1;
end % Main loop

% Output last chosen b
if trace > 2
    disp(sprintf('%5.0f   %13.6g %13.6g        %s',fcount, b, fb, procedure));
end
end %-------------------------------------------------------------- end findzero

end %%% END MNUM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MNUM %%%


%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%

function flsetup = flowsetup(T2,Tmax,theta,s,mem,f) %------------------- flsetup
%FLOWSETUP  Setup of flow properties.
% FS = FLOWSETUP(T2,Tmax,THETA,S,MEM,F) returns a struct FS that contains flow
% properties.
%  Fields:
%    FS.curv              Curvature.
%    FS.kelv(T,sigma,rho) pk/psat
%    FS.pkelv(T)          pk/psat
%    FS.hgK(T)            Specific enthalpy of the vapor, h(T,pk(T)).
%    FS.intdhdpdpsatdT(T) Int_T2^Tmax dh/dp dpsat/dT dT.
%    FS.nuapp(T,p)        Apparant vapor viscosity (viscous + molecular flow).
%    FS.knudsen(T,p)      Knudsen number.
%    FS.nu2ph(T,pk,a)     Apparent 2ph viscosity, using app. vapor viscosity.
%    FS.kmgas(T)          Heat conductivity of vapor-filled membrane.
%    FS.kmliq(T)          Heat conductivity of liquid-filled membrane.
%    FS.k2ph(T,a)         Heat conductivity of two-phase filled membrane.
%    FS.xdot(T,pk,a)      Vapor mass flow fraction.
%    FS.odemaxstep(range,delta) Minimum number of integration steps.

flsetup = struct('curv',[],'kelv',[],'pkelv',[],'hgK',[],'intdhdpdpsatdT',[],...
  'nuapp',[],'knudsen',[],'nu2ph',[],'kmgas',[],'kmliq',[],'k2ph',[],...
  'xdot',[],'odemaxstep',[]);

flsetup.curv = mem.fcurv(cos(theta*pi/180));
flsetup.kelv = @(T,sigma,rho) exp(-flsetup.curv.*sigma./(s.R.*rho.*T));

% Apparent vapor viscosity
%   nuapp = nug/(1 + beta Kn),  Kn = 3 nug sqrt(pi/(8*R*T)) / dia,
if mem.beta
  kn_nu = @(T) 3*sqrt(pi/(8*s.R)) / (sqrt(T)*mem.dia);
  flsetup.knudsen = @(T,p) kn_nu(T)*s.nug(T,p);
  flsetup.nuapp = @(T,p) 1 / (1/s.nug(T,p) + mem.beta*kn_nu(T));
  flsetup.nu2ph = @nu2phapp;
else % mem.beta = 0
  % Viscous flow
  flsetup.knudsen = @(T,p) 0;
  flsetup.nuapp = @(T,p) s.nug(T,p);
  flsetup.nu2ph = @(T,pk,a) f.nu2ph(a,s.v(T,pk),1/s.rho(T),s.mug(T),s.mul(T));
end

% Flow model and effective two-phase properties.
flsetup.kmgas = @(T) f.kmgas(mem.epsilon,mem.km,s.kg(T));
flsetup.kmliq = @(T) f.kmliq(mem.epsilon,mem.km,s.kl(T));
flsetup.k2ph = @(T,a) f.k2ph(a,mem.epsilon,mem.km,s.kg(T),s.kl(T));
flsetup.xdot = @(T,pk,a) f.xdot(a,s.v(T,pk),1/s.rho(T));
%flsetup.x = @(T,a) f.x(a,s.v(T,pk(T)),1/s.rho(T));

% Step size calculation for ode45, ode23t.
% minsteps (= 2): minimum number of integration steps;
% odemaxstep: number of (power of 2: 2, 4, 8, ...) steps such that T or p does
% not change more than maxTperstep (maxpperstep) within one step.
flsetup.odemaxstep = ... % 1/max(minsteps,...
  @(range,maxstep) 1/max(2,pow2(ceil(log2(abs(range)/maxstep))));

% Calculate hgK, intdhdpsdpsatdT and pkelv. Check if T > Tc, i.e., ps(T) == Inf.
% inttol: relative tolerance for integration
inttol = 1e-6;
stepT = 0.4;

% Return Inf or NaN if the fluid is anywhere supercritical.
% Might be changed with moderate effort to return functions that give exact
% values.
if ~isfinite(s.ps(T2))
  flsetup.pkelv = @(T) Inf;
  flsetup.hgK = @(T) NaN;
  flsetup.intdhdpdpsatdT = @(T) NaN;
  return
end

% do an integer number of steps
%intTrange = ceil(intTrange/stepT)*stepT;
intTrange = ceil((Tmax-T2)/stepT)*stepT;
if intTrange > 4
  intTrange = [T2 T2+4*intTrange]; %DEBUG use 2 for normal use
else
  intTrange = [T2 T2+10];
end
%intTrange = [295 305];

options = odeset('RelTol',inttol,'Refine',1,...
  'InitialStep',stepT,'MaxStep',stepT);
%solhgK = ode45(@inthgK,[Trange],initial condition,options);
solhgK = ode45(@inthgK,intTrange,0,options);
function dy = inthgK(T,y)
  % Calculate the enthalpy of the vapor along (T, pk(T)). h = 0 at T = T2.
  % With dh = (dh/dp) dp + (dh/dT) dT, dp = (dpk/dT)dT, (dh/dT) = cpg, follows
  %   dh/dT = (dh/dp) dpK/dT + cpg
  % One equation: Scaling (making dimensionless) is not necessary.
  [dpk pk] = dpkdT(T);
  [dhdp cpg] = s.dhcpg(T,pk);
  dy = dpk*dhdp + cpg;
end
% Quadrature, with fun(T) = dpk*dhdp + cpg and delhgK = quad(fun,T1,T2) would
% also have been possible. But quadrature would have to be evaluated each time.
% Instead, solhgK returns h(T).

% options stay the same as before
soldhdps = ode45(@intdhdpdpsdT,intTrange,0,options);
function dh = intdhdpdpsdT(T,h)
  % calculate the term
  %   T2_int^T dh/dp dpsat/dT' dT'
  % in the liquid phase, hence dh/dp = v - Tdv/dT = 1/rho + (1/rho^2) drho/dT.
  [drho rho] = s.drho(T); [ps dps] = s.ps(T);
  dh = dps*(1+drho)/rho;
end

flsetup.pkelv = @(T) s.ps(T)*flsetup.kelv(T,s.sigma(T),s.rho(T));
flsetup.hgK = @(T) deval(solhgK,T);
flsetup.intdhdpdpsatdT = @(T) deval(soldhdps,T);

function [dpk pk] = dpkdT(T) %-------------------------------------------- dpkdT
%DPKDT      Derivative of the equilibrium pressure at a curved meniscus, dpk/dT.
%
% [DPK PK] = DPKDT(T) returns the pk and dpk/dT. A copy from HVAPK, for
% convenience.
%
% See also FLOW12/HVAPK.
[psat dps] = s.ps(T);
[dsig sigma] = s.dsig(T);
[drho rho] = s.drho(T);
pk_ps = flsetup.kelv(T,sigma,rho);
pk = pk_ps*psat;
dpk = pk_ps * (dps + psat*flsetup.curv*sigma*(1/T-dsig+drho)/(s.R*rho*T));
end %----------------------------------------------------------------- end dpkdT

function nu2app = nu2phapp(T,pk,a) %----------------------------------- nu2phapp
%NU2PHAPP   Apparent viscosity (Knudsen + viscous flow) of the 2ph-mixture.
%
% NU2PHAPP(T,P,A) returns the apparent kinematic viscosity [m2/s] that results
% from assuming an apparent viscosity for the vapor flow,
%   nuapp = nug / (1 + beta*Kn).
vol = s.v(T,pk);
nu2app = f.nu2ph(a,vol,1/s.rho(T),flsetup.nuapp(T,pk)/vol,s.mul(T));
end %-------------------------------------------------------------- end nu2phapp
end %--------------------------------------------------------------- end flsetup

function fl = flowstruct(theta,s,mem,f) %---------------------------- flowstruct
%FLOWSTRUCT Constructs the structure FL containing the solution.
%  FLOWSTRUCT(THETA,SUBSTANCE,MEMBRANE,FMODEL)
%  FL contains the fields FL.INFO, FL.SOL and FL.FLOW.
%  FL.info contains
%    .theta
%    .substance
%    .membrane
%    .fmodel
%    .m
%    .T1
%    .p1
%    .T2
%    .p2
%    .kappac
%    .Kkappa	kappaK(kappa)
%    .Ccc
%  FL.sol contains
%    .m
%    .T1
%    .p1
%    .q1
%    .T2
%    .p2
%    .len
%    .states
%  FL.flow usually is of size > 1 and contains
%    .z
%    .T
%    .p
%    .Kn (for vapor flow in membrane, only)
%    .q
%    .a
%    .color

fl = struct('info',[],'sol',[],'flow',[]);
fl.info = struct('theta',theta,'substance',s,'membrane',mem,'fmodel',f,...
  'm',[],'T1',[],'p1',[],'T2',[],'p2',[],'kappac',[],'Kkappa',[],'Ccc',[]);
fl.sol = struct('m',[],'T1',[],'p1',[],'q1',[],'T3',[],'T2',[],'p2',[],...
  'len',[],'states',[]);
fl.flow = struct('z',{},'T',{},'p',{},'q',{},'a',{},'color',{});

end %------------------------------------------------------------ end flowstruct

%%% GARBAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GARBAGE %%%

%--------------------------------------------------------------------- garbage 1
%% SHOOT USING FZERO
%% Construct an interval, though; FZERO is very anxiuos. The interval search is
%% copiet from findzero.
%% Assign accurate tolerances and use fzero
%solver = accurate(solver);
%
%% Assign the function to zero in.
%% tolx: Tolerance in the mass flux. p1tol: Tolerance in pressure p1.
%pres = @(m) presiduum(m,solver);
%% See line 463 in fzero: toler = tol*max(abs(m),1.0)
%tolx = min(1,10^(floor(log10(mguess))))*1e-5;
%
%
%a(1) = mguess; fa(1) = presiduum(a(1),solver);
%
%%  SEARCH AN INTERVAL
%% Now look, in which direction we have to search for pres =  0.
%if fa(1) > 0, fac = 0.5; else fac = 1.4; end
%
%for i = 2:6
%  old = i - 1;
%  a(i) = a(old)*fac;
%  fa(i) = presiduum(a(i),solver);
%  if (fa(old) < 0) == (fa(i) > 0) || fa(i) == 0
%    % found an interval.
%    a = [a(i) a(old)];
%    fa = [fa(i) fa(old)];
%    break
%  end
%end
%
%%  PLOT SEARCH, IF NO INTERVAL FOUND
%if size(a,2) < 2
%  % No interval found.
%  % also plot the result for mvap and mliq
%  mliq = p12*mem.kappa/(s.nul(T1)*mem.L);
%  mvap = p12*mem.kappa/(s.nug(T2,(p1+p2)/2)*mem.L);
%  pvap = presiduum(mvap,solver);
%  pliq = presiduum(mliq,solver);
%  % No, really no interval found.
%  disp('MNUM: No interval found! mguess = %g',mguess);
%  disp('  Call MNUM with a guess for the mass flow,');
%  disp('  mnum(T1,p1,p2,theta,s,mem,f,''m'',mguess)');
%  plot(a,fa,'ko-')
%  plot(a,fa/p1,'ko-',mvap,pvap/p1,'r*',mliq,pliq/p1,'b*');
%  xlabel('mflux'); ylabel('p_{1,calc}/p_1 - 1');
%  return
%end
%
%% Setup and invoke fzero
%% 'TolFun' is not used in fzero. % 'TolFun',p12*1e-3
%% 'Display': 'iter', 'notify', 'off', 'final'
%% also, 'TolX' becomes TolX*max(abs(X),1) - not usable
%m = fzero(pres,a,optimset('TolX',tolx,'Display','iter'));
%
%% Control residuum in pressure
%while abs(fl.sol.p1 - p1) > p1tol
%  disp('MNUM: Tolerance in fzero had to be reduced.');
%  disp(sprintf('  TolX = %g, (p1_calc - p1_given)/(p1 - p2) = %g.',tolx,...
%    (fl.sol.p1 - p1)/p12));
%  mguess = [m-tolx m+tolx];
%  tolx = tolx*0.5*p1tol/abs(fl.sol.p1 - p1);
%  disp(sprintf('  New tolerance: TolX = %g. ',tolx));
%  m = fzero(pres,mguess,optimset('TolX',tolx,'Display','iter'));
%end
%fl.sol.m = m;
%----------------------------------------------------------------- end garbage 1
