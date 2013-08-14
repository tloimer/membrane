function b = findzero(presiduum,mguess,p1tol,solver)
%FINDZERO  Find the zero of a function PRESIDUUM(M,SOLVER)
%  M = FINDZERO(PRESIDUUM,MGUESS,PTOL,SOLVER) finds a zero to the function
%  PRESIDUUM(M,SOLVER). FINDZERO returns a solution where the function value is
%  smaller than PTOL. During iteration, the accuracy of the solver is increased.
%  FINDZERO uses a scalar value MGUESS as an initial guess, an interval MGUESS =
%  [M1 M2] as an interval in which the function changes sign, or an array MGUESS
%  = [M1 M2; PRES(M1) PRES(M2)] as an interval with the function values provided.
%  The provided function values must have been calculated with SOLVER('crude')!
%
%  See also FZERO, SOLVER, MNUMADIABAT.

%  INITIALIZE
% Information is plotted for trace > 1. Iteration results for trace > 2.
trace = 1;
global VERBOSE; if VERBOSE > 0, trace = trace + VERBOSE; end
savesolver = solver;

%  ASSIGN CRUDE SOLVER TOLERANCES
% During iteration, use crude tolerances.
solver = solverstruct('crude');
% defaults: solver.writesolution = false;
%	    solver.fullsolution = false; solver.partialsolution = true;

% INITIALIZE, DEPENDING ON MGUESS
if isscalar(mguess)
  % SEARCH AN INTERVAL
  a(1) = mguess; fa(1) = presiduum(a(1),solver); fcount = 1;

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
      break
    end
  end

  % Plot search, if no interval found
  if size(a,2) ~= 2
    % No interval found.
    fprintf('MNUM: No interval found! mguess = %g\n',mguess);
    fprintf('  Mass flow [kg/m2s]: '); disp(a);
    fprintf('  p_res [Pa]:  '); disp(fa);
    % reuse fac as figure handle
    fac = figure('Name','Search for Interval','NumberTitle','off');
    plot(a,fa/p1,'ko-');
    xlabel('mflux'); ylabel('p_{1,calc}/p_1 - 1');
    % AND ASK FOR AN INTERVAL
    while size(a,2)~= 2
      b = input(['  Provide an interval [mflux mflux] or try\n' ...
		 '  several points [mflux mflux mflux ...]: ']);
      lenb = size(b,2); fb = b;
      for i = 1:lenb,  fb(i) = presiduum(b(i),solver);  end
      % reuse old as figure handle; save a possible different current figure.
      old = get(0,'CurrentFigure');
      figure(fac);
      line(b,fb/p1,'Marker','+','LineStyle','none');
      figure(old);
      if lenb == 2 && (fb(1) < 0) == (fb(2) > 0)
	a = b; fa = fb;
      end
    end
    close(fac);
  end

  % Interval found or assigned
  b = a(2); a = a(1); fb = fa(2); fa = fa(1);
elseif size(mguess) == [1 2]
  % GOT AN INTERVAL
  a = mguess(1);  fa = presiduum(a,solver);
  b = mguess(2);  fb = presiduum(b,solver);
  fcount = 2;
  if (fa > 0) ~= (fb > 0)
    if fa == 0
      b = a;
      if trace > 1, fprintf('Hit a zero!\n'); end
      return
    end
    if fb == 0
      if trace > 1, fprintf('Found an exact zero!\n'); end
      return
    end
    error(['Interval endpoints do not differ in sign.\n mguess = '...
	   '[%.4g %.4g] kg/m2s, p1calc - p1 = [%.4g %.4g] Pa.\n'], a,b,fa,fb);
  end
elseif size(mguess) == [2 2]
  % AN INTERVAL WITH FUNCTION VALUES
  a = mguess(1,1);  fa = mguess(2,1);
  b = mguess(1,2);  fb = mguess(2,2);
  fcount = 0;
else
  error('Wrong size of mguess, [%d %d]. Must be scalar, [1 2] or [2 2].\n',...
	size(mguess));
end
lowprecision = true;

% if a == 0 or b == 0, the while-loop is run once, lowprecision set to false, and the
% result re-calculated.

% Use tolerance larger than computer epsilon, because toler below is computed
% with respect to abs(b), not max(abs(b),1).
tol = 2*eps;
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
    toler = 2.0*tol*abs(b); %max(abs(b),1.0);
    if (abs(fb) < p1tol) || (abs(m) <= toler) % || (fb == 0.0)
        if lowprecision
	  if trace > 1
	    fprintf(['p1calc - p1 = %g Pa, m = %g kg/m2s. We jumped over '...
		     '100xp1tol.\n'], fb, b);
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
		% avoid an infinite loop, see abs(m) above
    if lowprecision && ( abs(fb) < 100*p1tol || abs(m) <= toler )
      % Continue with the demanded solver tolerance
      solver = savesolver;
      lowprecision = false;
      if trace > 1
	fprintf(['p1calc - p1 = %g Pa, m = %g kg/m2s. Switched to '...
		 'required precision.\n'], fb, a);
      end
    end
    fb = presiduum(b,solver);
    fcount = fcount + 1;
end % Main loop

% Output last chosen b
if trace > 2
    fprintf('%5.0f   %13.6g %13.6g        %s\n',fcount, b, fb, procedure);
end

%%% END FINDZERO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END FINDZERO %%%
