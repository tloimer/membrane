function [m fl] = mnum(T1,p1,p2,theta,s,mem,f,Tmax,mguess)
%MNUM       Mass flux by numerical integration.
%  MNUM(T1,P1,P2,THETA,S,MEM,F) calculates the mass flux [kg/m2s] for a
%  substance S through a membrane MEM with a two-phase flow model F. The
%  upstream conditions are T1 [K] and P1 [Pa], the downstream pressure is P2
%  [Pa] and the contact angle is THETA, in degrees. Diagnostic output is printed
%  for the GLOBAL variable VERBOSE > 0.
%
%  [M FL] = MNUM(T1,P1,P2,THETA,S,MEM,F) returns the mass flux M and a
%  flowstruct FL containing the solution.
%
%  MNUM(T1,P1,P2,THETA,S,MEM,F,'m',MGUESS) uses MGUESS as the inital guess for
%  the mass flux.
%
%  MNUM(M,T2,P2,THETA,S,MEM,F,TMAX) constructs a struct SOLVER and invokes
%  FLOW12. Returns a FLOWSTRUCT. See FLOW12.
%
%  The flowstruct FL has the fields FL.info, Fl.calc, FL.sol and FL.flow.
%  FL.calc is written by mnum. For gases it is empty. For vapors Fl.calc
%  contains
%    FL.calc.psat1      psat(T1)
%    FL.calc.pK1        pK(T1)
%    FL.calc.n          n = (dT/dp)_h dpK/dT at T1, see Eq. (13) in TL, 2007. 
%    FL.calc.Ccc        (p1-p2)(1-n)/(psat-pK), Eq. (12) in Loimer, 2007.
%    FL.calc.Ccap       n(p1-p2)/pcap, Eq. (14) in Loimer, 2007.
%    FL.calc.kappac     kappa_crit  cf. Eq. (8) in Schneider, 1983.
%    FL.calc.kapK       dpK/dT instead of dpsat/dT in kappac.
%    FL.calc.kapl       (dpK/dT - d pcap/dT) instead of dpsat/dT in kappac.
%    FL.calc.kapKK      kapKK = kapK(kapK), pK and hvapK iteratively evaluated.
%    FL.calc.kapll      kapll = kapl(kapl), iteratively evaluated at kapl.
%    FL.calc.mlinp1sat  mass flux according to linear theory if p1 = psat(T1).
%    FL.calc.mguess     mass flux, estimated by improved linear description,
%                       see (Loimer, eurotherm09)
%    FL.calc.mgas       mass flux for isothermal flow of the gaseous phase
%
%  Calls FLOW12, FLOWSTRUCT, FLOWSETUP.
%
%  See also FLOW12, FMODEL, MEMBRANE, SUBSTANCE.
%
%  Subfunction:   FLCALCVARS.
%  Try, e.g., help mnum>flcalcvars.
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
[psat1 dps1] = s.ps(T1);
if  p1 > psat1
  error(['MNUM: Fluid is in liquid state at upstream location:'...
    ' p1 > psat(T1), p1 = %g bar, psat = %g bar'],p1*1e-5,s.ps(T1)*1e-5);
end
p12 = p1 - p2;
T2 = s.intjt(T1,p1,p2);
fl.info.T1 = T1; fl.sol.T2 = T2;
fl.info.p1 = p1;
flsetup = flowsetup(T2,T1,theta,s,mem,f);
fl.info.flsetup = flsetup;

% Calculate the information already possible.
if isfinite(psat1) % above the critical point, these numbers do not make sense
  fl.calc = flcalcvars(T1,p12,p1,psat1,dps1,flsetup,theta,s,mem,f);
end

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
  % calculate the flow of the fluid in its gaseous state, isothermal flow
  mgas = p12*mem.kappa/(flsetup.nuapp(T1,(p1+p2)/2)*mem.L);
  % psat1 is finite, the struct .calc is only written for vapors.
  if ~isempty(fl.calc) % this test could break easily!
    if theta <= 90
      plo = fl.calc.pK1 - fl.calc.n*p12;
      if p1 > plo
	%if fl.calc.Ccc >= 1
	  % Use linear approximation, see Loimer, eurotherm 2009.
	  % no further distinction between Ccc > 1 or Ccc > 1.
	  mguess = mgas + (fl.calc.mlinp1sat-mgas)*(p1-plo)/(psat1-plo);
	%else
	  % rather err towards the larger values
	  %  fl.caIlc.mlinear = fl.calc.mlinp1sat;
	%end
      else % p1 <= plo
	mguess = mgas;
      end % p1 > plo
    else % theta > 90
      % for non-wetting fluid near saturation use linear theory
      flin = mlinear(p1,p2,T1,theta,s,mem,f);
      mguess = flin.lin.m;
    end % theta
  else % isempty(fl.calc), therefore: psat1 = Inf
    % for unsaturated vapor or gas, use gaseous fluid flow
    mguess = mgas;
  end % ~isempty(fl.calc)
  fl.calc.mgas = mgas;
end
fl.calc.mguess = mguess;

%  SHOOT AT P1
p1tol = 1e-3*p12;
m = findzero(mguess,p1tol);
fl.sol.m = m;


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
global VERBOSE; if VERBOSE > 0, trace = trace + VERBOSE; end
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

%-------------------------------------------------------------------- flcalcvars
function calc = flcalcvars(T1,p12,p1,psat1,dps1,flsetup,theta,s,mem,f)
%FLCALCVARS Calculate div. kappas, Ccc, Ccap, pk1 and other fl.calc variables.

global VERBOSE;

% Initialize the output struct.
calc = struct('psat1',psat1,'pK1',[],'n',[],'Ccc',[],'Ccap',[],...
  'kappac',[],'kapK',[],'kapl',[],'kapKK',[],'kapll',[],'mlinp1sat',[]);

% We need this variables.
[dpk1 calc.pK1] = flsetup.dpkdT(T1);
[dsig1 sigma1] = s.dsig(T1);
[drho1 rho1] = s.drho(T1);
% probably that also
hvap1 = s.hvap(T1);

dsigma_dT1 = sigma1*dsig1;
nulkml = s.nul(T1)*flsetup.kmliq(T1);

calc.n = s.jt(T1,p1)*dpk1;
calc.kappac = nulkml/hvap1; % reuse kappac!
calc.kapK = calc.kappac/dpk1;
calc.kapl = calc.kappac/(dpk1-flsetup.curv*dsigma_dT1);
calc.kappac = calc.kappac/dps1; % kappac now is OK.
%
facdia = mem.dia/sqrt(mem.kappa); %  dia = facdia*sqrt(kappa);
faccurv = flsetup.curv*mem.dia; %   curv = faccurv/dia;
facpk_ps = -sigma1/(s.R*rho1*T1); % pk_ps = exp(curv*facpk_ps)
% pk_ps = exp(-flsetup.curv.*sigma./(s.R.*rho.*T));
% dpk = pk_ps * (dps1 + psat1*curv*sigma1*(1/T1-dsig1+drho1)/(s.R*rho1*T1));
facdpk = psat1*sigma1*(1/T1-dsig1+drho1)/(s.R*rho1*T1);
% dpk = pk_ps * (dps1 + curv*facdpk)
dhldp1 = (1 + T1*drho1)/rho1;
% hvK = hvap1 + (prad-psat)*(s.dhdp(T,prad)-dhldp) + dhldp*pcap;
fachvap = dhldp1*sigma1;

function [curv dpk hvapK] = varskappa(kappa)
  % curv = faccurv/dia; % dia = facdia*sqrt(kappa);
  curv = faccurv/(facdia*sqrt(kappa));
  pk_ps = exp(curv*facpk_ps);
  dpk = pk_ps * (dps1 + curv*facdpk);
  pk = pk_ps*psat1;
  hvapK = hvap1 + ( (pk-psat1)*s.dhdp(T1,pk)-dhldp1 ) + curv*fachvap;
end

function kapKK = kappaK(kappa)
  [curv dpk hvapK] = varskappa(kappa);
  kapKK = nulkml/(dpk*hvapK);
end

function kapll = kappal(kappa)
  [curv dpk hvapK] = varskappa(kappa);
  %dpcap = curv*dsigma_dT1;
  kapll = nulkml/((dpk-curv*dsigma_dT1)*hvapK);
end

options = optimset('fzero');
options = optimset(options,'TolX',1e-6*calc.kappac);
calc.kapKK = fzero(@(kappa) kappa - kappaK(kappa),calc.kapK,options);
if mem.kappa < calc.kapl
  range = calc.kapl*[0.99 1+0.3*calc.kapl/mem.kappa];
else
  range = calc.kapl./[1+0.3*mem.kappa/calc.kapl 1.01];
end

try
  calc.kapll = fzero(@(kappa) kappa - kappal(kappa),range,options);
catch err
  switch err.identifier
    case 'MATLAB:fzero:ValuesAtEndPtsSameSign'
      calc.kapll = fzero(@(kappa) kappa - kappal(kappa),calc.kapl,options);
      if VERBOSE > 0
        warning(['Fzero threw an error when trying to calculate kappa_ll,\n'...
	  'no change of sign between endpoints of the initial, guessed '...
	  'range.\nOnly specify a starting point.\n'...
	  'Initial range [%.3g %.3g], result %.3g'],...
	  range(1),range(2),calc.kapll);
      end
    case 'MATLAB:fzero:ValuesAtEndPtsComplexOrNotFinite'
      warning('Invalid kappa_ll (fl.calc.kapll), a non-wetting system?');
    otherwise
      warning('fl.calc.kapll could not be calculated.');
      %rethrow(err);
  end
end

  if theta ~= 90
    if theta < 90
      calc.Ccc = p12*(1-calc.n)/(psat1-calc.pK1);
    else
      calc.Ccap = calc.n*p12/(flsetup.curv*sigma1);
    end
  end
% Calculate the mass flux for p1 = psat(T1) according to linear theory.
  flin = mlinear(psat1,psat1-p12,T1,theta,s,mem,f);
  calc.mlinp1sat = flin.lin.m;
end %------------------------------------------------------------ end flcalcvars
