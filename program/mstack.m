function [m mguess fl] = mstack(T1,p1,p2,theta,s,mem,f)
%MSTACK     Mass flux for a stack of membranes.
%  MSTACK(T1,P1,P2,THETA,S,MEM,F) returns the mass flux [kg/m2s] for a stack of
%  N independent membranes MEM(1,I:N). MEM is a row vector of membranes
%  MEMBRANE. MEM(1) is the first membrane, on which P1 is applied. A GLOBAL
%  variable VERBOSE controls display of information, if VERBOSE > 0.
%
%  CALLS FLOW12, ...
%
%  Subfunctions: SOLVERSETUP, FINDZERO (functional, but not used)
%
%  Nested functions: PRESIDUUM

if ~isscalar(T1) || ~isscalar(p1) || ~isscalar(p2) || ~isscalar(s)
  error([upper(mfilename) ': Scalar values expected for T1, p1, p2 and s.']);
end

global VERBOSE
fzerooutput = 'notify';

% Make theta and f row vectors of size(mem), (1,1:mem). If theta or f are
% scalar, the following assignment expands them to the desired vectors. If they
% are vectors of the same size, the assignment does nothing. If the sizes do not
% match, an error occurs.
nmem = size(mem,2);
theta(1,1:nmem) = theta;
f(1,1:nmem) = f;

% Calculate an estimate for the mass flux.
nucirca = s.nug(T1,(p1+p2)/2);
T2 = s.intjt(T1,p1,p2);
mguess = (p1-p2)/(nucirca*sum([mem.L]./[mem.kappa]));
if VERBOSE > 0
  fprintf([upper(mfilename) ': mguess = %.3g.\n'],mguess);
  fzerooutput = 'final';
end

% Assign values to flowstruct and flowsetup. Tmin = T2, Tmax = T1 + (T1-T2)/5;
T12 = T1 - T2;
%disp(sprintf('Temperature difference: %.3f, range %.3f -- %.3f',T12,T2,T2));
%Tmin = T2 - max(2,2*T12); Tmax = T1 + max(2,2*T12);
Tmin = T2; Tmax = T1 + max(2,3*T12);
%fprintf('Valid range %.3f -- %.3f.\n',Tmin,Tmax);
for i = nmem:-1:1 % allocate array of correct size in the first loop
  flsetup(i) = flowsetup(Tmin,Tmax,theta(i),s,mem(i),f(i));
  fl(i) = flowstruct(theta(i),s,mem(i),f(i));
end

% Allocate space for p and T between the membranes.
pi = zeros(1,nmem);
Ti = pi;

% Setup solver with crude tolerances, and find an interval.
solver = solversetup('crude');   solver.partialsolution = false;
[minterval pinterval] = findinterval(@(m) presiduum(m,solver),mguess,p2-p1);

% Check, if pinterval has a too high pressure; if so, shoot with one guess
% Gave an error in flow12/flow45
% if isempty(s.Ts(pinterval(2))) || isempty(s.Ts(pinterval(1)))
%   warning([upper(mfilename) ':FINDINTERVAL'],...
%     'Call to findinterval returned p > pcrit');
%   % Set minterval to the interpolated value. See FINDINTERVAL.M for the formula.
%   minterval =  minterval(2) - ...
%     pinterval(2)*(minterval(2)-minterval(1))/(pinterval(2)-pinterval(1));
%   fzerooutput='iter';
% end

%% Shoot with findzero. Would work.
%p1tol = 1e-3*(p1-p2);
%m = findzero(@presiduum,mguess,solver,p1tol);

% Setup solver with fine tolerances.
solver = solversetup('accurate');   solver.partialsolution = false;

% Setup and shoot with fzero.
options = optimset('TolX',mguess/1000,'Display',fzerooutput);
[m,p1res,~,output] = fzero(@(m) presiduum(m,solver),minterval,options);

% Print diagnosis.
if VERBOSE > 0
  fprintf([upper(mfilename) ': Calls to flow12: %u,\n'],output.funcCount);
  fprintf([upper(mfilename)...
    ': Iterations to find an interval: %u, Iterations to zero in: %u.\n'],...
    output.intervaliterations,output.iterations);
  fprintf([upper(mfilename) ': Mass flow fzero %.3g, '...
  'residuum %.3g Pa, or %.3g%%.\n'], m,p1res,100*abs(p1res)/(p1-p2));
end

%% Write the output struct
%sol = struct('mguess',mguess);

%%% NESTED FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NESTED FUNCTIONS %%%
% Nested functions, because we use variables from the main function. Would be
% cleaner, to use subfunctions.

function pres = presiduum(m,solver) %--------------------------------- presiduum
% Input variables:  m, T2, p2, nmem(?), fl, flsetup, solver
% Output variables: pres, pi, Ti, (fl)

if m <=0 % Provide for fzero straying into negative terrain.
  pi = p2; Ti = T2;
  pres = p2 - p1;
return
end

% We are save, compute.
plo=p2; Tlo=T2;
for i = nmem:-1:1
  fl(i) = flow12(m,Tlo,plo,fl(i),flsetup(i),solver);
  plo = fl(i).sol.p1; Tlo = fl(i).sol.T1;
  pi(i) = plo; Ti(i) = Tlo;
end

pres = plo - p1;
end %------------------------------------------------------------- end presiduum

%%% END NESTED FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END NESTED FUNCTIONS %%%

end %%% END MSTACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MSTACK %%%


%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%

function solver = solversetup(solvertype) %------------------------- solversetup
%SOLVERSETUP Set the properties of the ode-solvers in flow12.
%  SOLVERSETUP returns a struct SOLVER set up for crude tolerances
%
%  SOLVERSETUP('accurate','writesolution') sets accurate solver tolerances. In
%  addition, FLOW12 writes the solution to the flow struct FL. Possible
%  arguments are 'crude' and 'accurate'.
%
%  By default, SOLVER.writesolution = false and SOLVER.partialsolution = true.
%
%  See also FLOW12.
solver = struct('rtol',1e-3,'atol',1e-5,'tola',1e-6,'maxTperstep',2,...
  'maxpperstep',1e5,'writesolution',true,'partialsolution',true); %%CHANGED!
%  'maxpperstep',1e5,'writesolution',false,'partialsolution',true); %%CHANGED!
if strcmp(solvertype,'accurate')
  solver.rtol = 5e-2;  solver.atol = 1e-5;
  solver.maxTperstep = .4;  solver.maxpperstep = 4e4;
end
end %----------------------------------------------------------- end solversetup

%%% OBSOLETE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OBSOLETE %%%

function b = findzero(shoot,mguess,requestedsolver,p1tol) %------------ findzero
%FINDZERO  Find a solution to FUN(X,SOLVER) == 0.

%  INITIALIZE
% Information is plotted for trace > 1. Iteration results for trace > 2.
trace = 3;
procedure = ' ';

%  ASSIGN CRUDE SOLVER TOLERANCES
% During iteration, use crude tolerances.
solver = solversetup('crude');
solver.partialsolution = false;

a(1) = mguess; fa(1) = shoot(a(1),solver); fcount = 1;

%  SEARCH AN INTERVAL
% Now look, in which direction we have to search for pres =  0.
if fa(1) > 0, fac = 0.5; else fac = 1.4; end

for i = 2:6
  old = i - 1;
  a(i) = a(old)*fac;
  fa(i) = shoot(a(i),solver);
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
  %mvap = p12*mem.kappa/(flsetup.nuapp(T2,(p1+p2)/2)*mem.L);
  %pvap = shoot(mvap,solver);
  % No, really no interval found.
  fprintf([upper(mfilename) ': No interval found! mguess = %g.\n'],mguess);
  %disp('  Call MNUM with a guess for the mass flow,');
  %disp('  mnum(T1,p1,p2,theta,s,mem,f,''m'',mguess)');
  disp('  Mass flow:'); disp(a);
  disp('  p_res :'); disp(fa);
  % reuse fac as figure handle
  fac = figure('Name','Search for Interval','NumberTitle','off');
  %plot(a,fa/p1,'ko-',mvap,pvap/p1,'r*'); %,mliq,pliq/p1,'b*');
  plot(a,fa,'ko-');
  xlabel('mflux'); ylabel('p_{1,calc}-p_1');
  % AND ASK FOR AN INTERVAL
  while size(a,2)~= 2
    b = input(['  Provide an interval [mflux mflux] or try\n' ...
      '  several points [mflux mflux mflux ...]: ']);
    lenb = size(b,2); fb = b;
    for i = 1:lenb,  fb(i) = shoot(b(i),solver);  end
    % reuse old as figure handle; save a possible different current figure.
    old = get(0,'CurrentFigure'); figure(fac);
    line(b,fb,'Marker','+','LineStyle','none');
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
            fprintf('We jumped over 100 x p1tol.\n');
          end
        else
          break
        end
    end
    if trace > 2
        fprintf('%5.0f   %13.6g %13.6g        %s\n',fcount, b, fb, procedure);
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
      solver = requestedsolver;
      lowprecision = false;
      if trace > 1
        fprintf('Residuum = %g Pa. Switched to requested solver.\n', fb);
      end
    end
    fb = shoot(b,solver);
    fcount = fcount + 1;
end % Main loop

% Output last chosen b
if trace > 2
    fprintf('%5.0f   %13.6g %13.6g        %s\n',fcount, b, fb, procedure);
end
end %-------------------------------------------------------------- end findzero
